library(dplyr)
library(reshape2)
library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)

load("protsData.RData")

prots= t(df)
metaFile = subset(meta, meta$time == "T1")

protsT1 = subset(prots, rownames(prots) %in% metaFile$name)
proteins = colnames(protsT1)

proteins = gsub("HLA-E","HLA_E",proteins)
colnames(protsT1) = gsub("HLA-E","HLA_E",colnames(protsT1))

protsT1 = merge(protsT1,metaFile,by.x=0,by.y="name")


####strain B####
protAssoc_B = data.frame(Protein = character(),IndependantVariable = character(),estimate = double(),std.error = double(),statistic = double(),p.value = double())

for(i in proteins) {
  print(i)
  vals = unlist(protsT1[i])
  t = (tidy(lm( vals ~ gender + age + log2(ab_B) ,data = protsT1))) %>%
    select(term,estimate,std.error,statistic,p.value) %>%
    as.data.frame()
  protAssoc_B[nrow(protAssoc_B)+1,] <- c(i,t[4,])
  
}

###### strain H1N1 #####
protAssoc_H1N1 = data.frame(Protein = character(),IndependantVariable = character(),estimate = double(),std.error = double(),statistic = double(),p.value = double())

for(i in proteins) {
  print(i)
  vals = unlist(protsT1[i])
  t = (tidy(lm( vals ~ gender + age + log2(ab_H1N1) ,data = protsT1))) %>%
    select(term,estimate,std.error,statistic,p.value) %>%
    as.data.frame()
  protAssoc_H1N1[nrow(protAssoc_H1N1)+1,] <- c(i,t[4,])
  
}

##### strain H3N2 #####
protAssoc_H3N2 = data.frame(Protein = character(),IndependantVariable = character(),estimate = double(),std.error = double(),statistic = double(),p.value = double())

for(i in proteins) {
  print(i)
  vals = unlist(protsT1[i])
  t = (tidy(lm( vals ~ gender + age + log2(ab_H3N2) ,data = protsT1))) %>%
    select(term,estimate,std.error,statistic,p.value) %>%
    as.data.frame()
  protAssoc_H3N2[nrow(protAssoc_H3N2)+1,] <- c(i,t[4,])
  
}