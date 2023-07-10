library(MOFA2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)

load("MofaData.RData")

finalDF = objList$longDF

MOFAobject <- create_mofa(finalDF)
print(MOFAobject)
plot_data_overview(MOFAobject)
data_opts <- get_default_data_options(MOFAobject)
head(data_opts)
data_opts$scale_views = TRUE

model_opts <- get_default_model_options(MOFAobject)
head(model_opts)
train_opts <- get_default_training_options(MOFAobject)
head(train_opts)


MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path(getwd(),"MOFAobject.trained.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile)

metaData = objList$metaData
colnames(metaData) = c("sample","patientID","Gender","Category","Time","Age") ##### rename the column names to avoid issues
samples_metadata(MOFAobject.trained) = metaData

head(MOFAobject.trained@cache$variance_explained$r2_total[[1]])
variance_perFactor = MOFAobject.trained@cache$variance_explained$r2_per_factor
variance_perFactorMelted = melt(variance_perFactor)

