library(modeldata)
library(pls)
library(tidymodels)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(RankProd)
library(dplyr)
library(maditr)

#### load("protsData.RData")
########Since we are working with variances and covariances, we need to standardize the data.###########
###### the ab values are log2 already.######
norm_rec <- 
  recipe(log2(H1N1) + log2(H3N2) + log2(B) ~ ., data = protsT1_PLS) %>%
  step_normalize(everything()) 

####Before we can finalize the PLS model, the number of PLS components to retain must be determined######
set.seed(57343)
folds <- vfold_cv(protsT1_PLS, repeats = 10)

######The folds can be created using the rsample ######
#######package and the recipe can be estimated for each resample using the prepper() function ####

folds <- 
  folds %>%
  mutate(recipes = map(splits, prepper, recipe = norm_rec))

#####PLS using tidymodel ######

get_var_explained <- function(recipe, ...) {
  
  # Extract the predictors and outcomes into their own matrices
  y_mat <- bake(recipe, new_data = NULL, composition = "matrix", all_outcomes())
  x_mat <- bake(recipe, new_data = NULL, composition = "matrix", all_predictors())
  
  # The pls package prefers the data in a data frame where the outcome
  # and predictors are in _matrices_. To make sure this is formatted
  # properly, use the `I()` function to inhibit `data.frame()` from making
  # all the individual columns. `pls_format` should have two columns.
  pls_format <- data.frame(
    endpoints = I(y_mat),
    measurements = I(x_mat)
  )
  # Fit the model
  mod <- plsr(endpoints ~ measurements, data = pls_format,na.action = na.exclude)
  
  # Get the proportion of the predictor variance that is explained
  # by the model for different number of components. 
  xve <- explvar(mod)/100 
  #coefficients <- mod$coefficients
  modLoadings = loadings(mod)
  # To do the same for the outcome, it is more complex. This code 
  # was extracted from pls:::summary.mvr. 
  explained <- 
    drop(pls::R2(mod, estimate = "train", intercept = FALSE)$val) %>% 
    # transpose so that components are in rows
    t() %>% 
    as_tibble() %>%
    # Add the predictor proportions
    mutate(predictors = cumsum(xve) %>% as.vector(),
           components = seq_along(xve)) %>%
    # Put into a tidy format that is tall
    pivot_longer(
      cols = c(-components),
      names_to = "source",
      values_to = "proportion"
    )
  
  return(list(
    'explained' = explained, 
    'loadings' = modLoadings
  ))
}

####We compute this data frame for each resample and save the results in the different columns.
lapply(X = folds$recipes, FUN = get_var_explained) -> info

# Variances explained per component
vars <- lapply(info, function(x) {x$explained} )
for (i in 1:length(vars)) {vars[[i]]$fold <- i}
vars <- do.call(rbind, vars)

vars %<>% 
  group_by(source, components) %>% 
  summarise(meanProp = mean(proportion), sd = sd(proportion)) %>% 
  filter(components <= 15)


# Coefficients fitted per protein per iteration of the model
#coefs <- lapply(info, function(x) {x$coefficients})
Allloadings <- lapply(info,function(x) as.data.frame(unclass({x$loadings})))

# At this point, coefs is a list of 3D arrays containing: protein x predictor x component (e.g. 100 iterations of 311 proteins x 3 predictors x 100 components)
# We want the mean and sd of each coefficient per protein per strain over all iterations, for one specific component
#COMP=1:3 # For plotting only
#lapply(coefs, function(iteration, n) {iteration[,,n] %>% as.data.frame()}, n=COMP) -> coefs
#for (i in 1:length(coefs)) {coefs[[i]]$iteration = i}

COMP=1:3 # For plotting only
lapply(Allloadings, function(iteration, n) {iteration[,n] %>% as.data.frame()}, n=COMP) -> Allloadings
for (i in 1:length(Allloadings)) {Allloadings[[i]]$iteration = i}

lapply(Allloadings, function(x) {x %>% tibble::rownames_to_column('protein') }) %>% 
  do.call(rbind, .) %>% as.data.frame() %>% 
  melt(c('protein', 'iteration')) -> Allloadings


######## RankProduct analysis ##########
Comp1 = subset(Allloadings,Allloadings$variable == "Comp 1")
Comp1Dcast = dcast(Comp1, protein + variable ~ iteration, value.var = "value")
Comp1Dcast = data.frame(Comp1Dcast)
rownames(Comp1Dcast) = Comp1Dcast$protein
Comp1Dcast$protein = NULL
Comp1Dcast$variable = NULL
cl <- rep(1,100)
Comp1Dcast_RP.out = RankProducts(abs(Comp1Dcast),cl,rand = 123,gene.names = rownames(Comp1Dcast))
Comp1Dcast_RP_merged = merge(Comp1Dcast_RP.out$RPrank,Comp1Dcast_RP.out$pfp,by=0)
colnames(Comp1Dcast_RP_merged) = c("ProtName","blah1","Rank","blah2","PFP")
Comp1Dcast_RP_merged$PFP = ifelse(Comp1Dcast_RP_merged$PFP == 0.0 , 2.225070e-308, Comp1Dcast_RP_merged$PFP)
Comp1Dcast_RP_merged$negLog10 = -(log10(Comp1Dcast_RP_merged$PFP))


