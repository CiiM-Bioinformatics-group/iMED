try(dev.off())
rm(list = ls())

library(dplyr)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggsci)
library(OlinkAnalyze)

# Cohort 1
df <- read_NPX(filename = '20210222_Li_NPX_2021-05-18.csv')
meta <- read.csv('meta_cohort1.csv', row.names = 1)
annot <- df %>% distinct(Assay, .keep_all = T) %>% select(OlinkID, UniProt, Assay)

df2 <- df %>%
  filter(SampleID %in% meta$name) %>%
  filter(Assay_Warning == 'PASS') %>%
  filter(QC_Warning == 'PASS') %>%
  filter(MissingFreq < 0.30) %>%
  reshape2::dcast(data = ., SampleID ~ OlinkID, value.var = 'NPX') %>%
  tibble::column_to_rownames(var = 'SampleID')

annot %<>%
  filter(OlinkID %in% colnames(df)) %>%
  arrange(match(OlinkID, colnames(df)))


stopifnot(all(rownames(df) == meta$name))
stopifnot(all(colnames(df) == annot$OlinkID))

save.image('output/cohort1.RData')



# Cohort 2
df <- read_NPX(filename = '20210222_Li_NPX_2021-05-18.csv')
meta <- read.csv('meta_cohort2.csv', row.names = 1)
annot <- df %>% distinct(Assay, .keep_all = T) %>% select(OlinkID, UniProt, Assay)

df %<>%
  filter(SampleID %in% meta$name) %>%
  filter(Assay_Warning == 'PASS') %>%
  filter(QC_Warning == 'PASS') %>%
  filter(MissingFreq < 0.30) %>%
  reshape2::dcast(data = ., SampleID ~ OlinkID, value.var = 'NPX') %>%
  tibble::column_to_rownames(var = 'SampleID')

annot %<>%
  filter(OlinkID %in% colnames(df)) %>%
  arrange(match(OlinkID, colnames(df)))


stopifnot(all(rownames(df) == meta$name))
stopifnot(all(colnames(df) == annot$OlinkID))

save.image('output/cohort2.RData')
