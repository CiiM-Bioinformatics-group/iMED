try(dev.off())
rm(list = ls())

library(dplyr)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(rsample)
library(purrr)
library(recipes)
library(reshape2)
library(ggsci)
library(pls)

load('data_flt.Rdata')

# Response matrix: log2 normalised Ab fold changes (3 columns x 200 rows)
# Predictor matrix: log2 normalised metabolic features for one specific timepoint (836 columns x 200 rows)

meta %>% filter(time == 'T3') %>%
  select(ab_H1N1, ab_H3N2, ab_B, ProbandID) %>%
  tibble::column_to_rownames('ProbandID') -> response.matrix
response.matrix$ab_H1N1 <- log2(response.matrix$ab_H1N1)
response.matrix$ab_H3N2 <- log2(response.matrix$ab_H3N2)
response.matrix$ab_B <- log2(response.matrix$ab_B)
response.matrix <- as.matrix(response.matrix)

t1.meta <- meta %>% filter(time == 'T1')
predictor.matrix <- t(df)
predictor.matrix <- predictor.matrix[t1.meta$name, ]
# colnames(predictor.matrix) <- make.unique(annot$finalname)
colnames(predictor.matrix) <- make.unique(annot$hmdb)
rownames(predictor.matrix) <- t1.meta$ProbandID

# Equal order
predictor.matrix <- predictor.matrix[rownames(response.matrix), ]
stopifnot(all(rownames(predictor.matrix) == rownames(response.matrix)))
#
# # 1. PLSR
# # https://stackoverflow.com/questions/71728756/how-to-use-r-package-caret-to-run-plsplsr-with-multiple-responses
# data <- data.frame(
#   'responses' = I(response.matrix),
#   'predictors' = I(predictor.matrix)
# )
#
# # plsr takes a data frame of matrices as input
# fit <- plsr(formula = responses ~ predictors, data = data, ncomp=50, method = 'kernelpls')
# validationplot(fit, estimate = "all", legendpos = "top")
# validationplot(fit, val.type = "R2")
# plot(fit, "validation", val.type = "MSEP", legendpos = "top")



# 2. Tidymodels
data <- cbind(response.matrix, predictor.matrix)
norm_rec <- recipe(ab_B + ab_H3N2 + ab_H1N1 ~ ., data = data) %>% step_normalize(everything())
folds <- vfold_cv(data, repeats = 10)
folds <- folds %>% mutate(recipes = map(splits, prepper, recipe = norm_rec))


get_info <- function(recipe, ...) {

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
  mod <- plsr(endpoints ~ measurements, data = pls_format, ncomp=100)

  # Get the coefficients
  coefficients <- mod$coefficients

  # Get the proportion of the predictor variance that is explained
  # by the model for different number of components.
  xve <- explvar(mod)/100

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
    tidyr::pivot_longer(
      cols = c(-components),
      names_to = "source",
      values_to = "proportion"
    )

  return(list(
    'explained' = explained,
    'coefficients' = coefficients
  ))
}

lapply(X = folds$recipes, FUN = get_info) -> info

# Variances explained per component
vars <- lapply(info, function(x) {x$explained} )
for (i in 1:length(vars)) {vars[[i]]$fold <- i}
vars <- do.call(rbind, vars)

vars %<>%
  group_by(source, components) %>%
  summarise(mean = mean(proportion), sd = sd(proportion)) %>%
  filter(components <= 30)

# png('output/var_explained_PLS.png', width = 5, height = 5, units = 'in', res = 300)
pdf('output/var_explained_plsr.pdf', width = 2, height = 2)
ggplot(vars, aes(x=components, y=mean, color = source)) +
  geom_point() +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd)) +
  geom_line() +
  theme_bw() +
  ggsci::scale_color_nejm() +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = .5), aspect.ratio = 1) +
  labs(x = 'Component', y = 'Mean', title = 'Explained variance PLS', color = 'Source')
dev.off()

# Coefficients fitted per metabolite per iteration of the model
coefs <- lapply(info, function(x) {x$coefficients})

# At this point, coefs is a list of 3D arrays containing: metabolite x predictor x component (e.g. 100 iterations x 836 metabolites x 3 predictors x 100 components)
# We want the mean and sd of each coefficient per metabolite per strain over all iterations, for one specific component
COMP=1 # For plotting only
lapply(coefs, function(iteration, n) {iteration[, , n] %>% as.data.frame()}, n=COMP) -> coefs
for (i in 1:length(coefs)) {coefs[[i]]$iteration = i}

lapply(coefs, function(x) {x %>% tibble::rownames_to_column('metabolite') }) %>%
  do.call(rbind, .) %>% as.data.frame() %>%
  melt(c('metabolite', 'iteration')) -> coefs

coefs %>% group_by(metabolite, variable) %>% summarise(mean = mean(value), sd = sd(value)) -> vals

# Sort by one of the strains
vals %>% filter(variable == 'ab_B') %>% arrange(mean) %>% pull(metabolite) %>% as.factor(.) -> ord
vals$metabolite <- factor(vals$metabolite, levels = ord)


pdf('output/importance_metabolites_PLS_component1.pdf', width = 8, height = 4)
ggplot(vals) +
  geom_point(aes(x = metabolite, y = mean), size = .5) +
  # geom_linerange(aes(x = metabolite, ymin = mean - sd, ymax = mean + sd)) +
  facet_grid('variable') +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(), plot.title = element_text(hjust = .5)) +
  labs(x = 'Ranked metabolites', y = 'Coefficient', title = 'PLS Component 1')
dev.off()
write.csv(x = vals, file = 'output/PLS_component1_importancescores.csv')


# Write the first five components and put the correct names
names <- annot %>% select(hmdb, finalname)

exp <- list()
for (comp in 1:5) {
  print(comp)
  coefs <- lapply(info, function(x) {x$coefficients})

  lapply(coefs, function(iteration, n) {iteration[, , n] %>% as.data.frame()}, n=COMP) -> coefs
  for (i in 1:length(coefs)) {coefs[[i]]$iteration = i}
  lapply(coefs, function(x) {x %>% tibble::rownames_to_column('metabolite') }) %>%
    do.call(rbind, .) %>% as.data.frame() %>%
    melt(c('metabolite', 'iteration')) -> coefs

  coefs %>% group_by(metabolite, variable) %>% summarise(mean = mean(value), sd = sd(value)) -> vals
  vals$component <- comp
  vals <- merge(vals, names, by.x = 'metabolite', by.y = 'hmdb')
  exp[[comp]] <- vals
}

do.call(rbind, exp) -> exp

write.csv(x = exp, 'output/PLS_5components_scores.csv')


# Rank Product stuff
Comp1 = subset(Allloadings,Allloadings$variable == "Comp 1")
Comp1Dcast = dcast(Comp1, metabolite + variable ~ iteration, value.var = "value")
Comp1Dcast = data.frame(Comp1Dcast)
rownames(Comp1Dcast) = Comp1Dcast$metabolite
Comp1Dcast$metabolite = NULL
Comp1Dcast$variable = NULL
cl <- rep(1,100)
Comp1Dcast_RP.out = RankProducts(abs(Comp1Dcast),cl,rand = 123,gene.names = rownames(Comp1Dcast))
Comp1Dcast_RP_merged = merge(Comp1Dcast_RP.out$RPrank,Comp1Dcast_RP.out$pfp,by=0)
colnames(Comp1Dcast_RP_merged) = c("ProtName","blah1","Rank","blah2","PFP")
Comp1Dcast_RP_merged$PFP = ifelse(Comp1Dcast_RP_merged$PFP == 0.0 , 2.225070e-308, Comp1Dcast_RP_merged$PFP)
Comp1Dcast_RP_merged$negLog10 = -(log10(Comp1Dcast_RP_merged$PFP))
