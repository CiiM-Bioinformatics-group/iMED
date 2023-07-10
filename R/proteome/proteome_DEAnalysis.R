
library(edgeR)
library(limma)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)

metaFile = read.csv("cohort2_metadata.csv",header = T)
load("protsData.RData")

metaFile = metaFile[match(rownames(df),metaFile$name),]#### to make sure the order of metadataFile is same as the matrix file
prots= t(df)

######## model for differential abundance analysis#######
mod = model.matrix(~ gender + age + responder.x*time, data=metaFile)

corfit = duplicateCorrelation(prots,mod,block = metaFile$ProbandID)
corfit$consensus.correlation
fit = lmFit(prots,mod,block=metaFile$ProbandID,correlation = corfit$consensus.correlation)
fit = eBayes(fit)
summary(decideTests(fit))


###### analysis over time per category#######
metaFile$responderTime = paste0(metaFile$responder.x,"_",metaFile$time)
modTime = model.matrix(~ 0 + responderTime+gender + age,data=metaFile)

corfit = duplicateCorrelation(prots,mod,block = metaFile$ProbandID)
corfit$consensus.correlation
fitTime = lmFit(prots,modTime,block=metaFile$ProbandID,correlation = corfit$consensus.correlation)
fitTime = eBayes(fitTime)
summary(decideTests(fitTime))
cont.matrix = makeContrasts(TR_T2vsT1 = responderTimeTR_T2 - responderTimeTR_T1,
                            TR_T3vsT1 = responderTimeTR_T3 - responderTimeTR_T1 ,
                            NR_T2vsT1 = responderTimeNR_T2 -responderTimeNR_T1,
                            NR_T3vsT1 = responderTimeNR_T3 -responderTimeNR_T1, levels = modTime)

fit2 = contrasts.fit(fitTime,cont.matrix)
fit2 = eBayes(fit2)
summary(decideTests(fit2))
