library(lme4)
library(fgsea)
library(lmerTest)
library(broom.mixed)


#### TR T3 vs T1 BTM score calculation #######
TR_T3vsT1_BTM_signi_split = split(TR_T3vsT1_BTM_signi,TR_T3vsT1_BTM_signi$pathway)
TR_T3vsT1_BTM_signi_genes = lapply(TR_T3vsT1_BTM_signi_split,"[",,8)


TR_T3vsT1_BTM_signi_meanPerSample = list()
for (i in 1:length(TR_T3vsT1_BTM_signi_genes)){
  temp = unname(unlist(TR_T3vsT1_BTM_signi_genes[[i]]))
  #print(names(TR_T3vsT1_BTM_signi_genes[i]))
  expr = subset(TR_T3T1_expr,rownames(TR_T3T1_expr) %in% temp)
  subt_expr = expr - rowMeans(expr)
  mean_expr = colMeans(subt_expr)
  TR_T3vsT1_BTM_signi_meanPerSample[[i]] = mean_expr
}
names(TR_T3vsT1_BTM_signi_meanPerSample) = names(TR_T3vsT1_BTM_signi_genes)

#### creating matrix of the list ######
TR_T3vsT1_meanMat_BTM = as.matrix(do.call(rbind,TR_T3vsT1_BTM_signi_meanPerSample))
TR_T3vsT1_BTM_MeanDF = t(TR_T3vsT1_meanMat_BTM)

temp = colnames(TR_T3vsT1_BTM_MeanDF)
temp2 <- stringr::str_extract(string = temp, pattern = "(?<=\\()[SM]\\d.*(?=\\))")
temp3 = gsub("\\s*\\([^\\)]+\\)\\s*$","",temp)
TR_T3vsT1_BTM_list = list()
for (i in 1:length(temp2)){ TR_T3vsT1_BTM_list[temp2[i]]<- temp3[i]}
colnames(TR_T3vsT1_BTM_MeanDF) = temp2


######NRs T3 vs T1 BTM score calculation#####
NR_T3vsT1_BTM_signi_split = split(NR_T3vsT1_BTM_signi,NR_T3vsT1_BTM_signi$pathway)
NR_T3vsT1_BTM_signi_genes = lapply(NR_T3vsT1_BTM_signi_split,"[",,8)


NR_T3vsT1_BTM_signi_meanPerSample = list()
for (i in 1:length(NR_T3vsT1_BTM_signi_genes)){
  temp = unname(unlist(NR_T3vsT1_BTM_signi_genes[[i]]))
  #print(names(NR_T3vsT1_BTM_signi_genes[i]))
  expr = subset(NR_T3T1_expr,rownames(NR_T3T1_expr) %in% temp)
  subt_expr = expr - rowMeans(expr)
  mean_expr = colMeans(subt_expr)
  NR_T3vsT1_BTM_signi_meanPerSample[[i]] = mean_expr
}
names(NR_T3vsT1_BTM_signi_meanPerSample) = names(NR_T3vsT1_BTM_signi_genes)

#### creating matrix of the list ######
NR_T3vsT1_meanMat_BTM = as.matrix(do.call(rbind,NR_T3vsT1_BTM_signi_meanPerSample))
NR_T3vsT1_BTM_MeanDF = t(NR_T3vsT1_meanMat_BTM)

temp = colnames(NR_T3vsT1_BTM_MeanDF)
temp2 <- stringr::str_extract(string = temp, pattern = "(?<=\\()[SM]\\d.*(?=\\))")
temp3 = gsub("\\s*\\([^\\)]+\\)\\s*$","",temp)
NR_T3vsT1_BTM_list = list()
for (i in 1:length(temp2)){ NR_T3vsT1_BTM_list[temp2[i]]<- temp3[i]}
colnames(NR_T3vsT1_BTM_MeanDF) = temp2

temp = TR_T2vsT1_BTM_signi$pathway
temp2 <- stringr::str_extract(string = temp, pattern = "(?<=\\()[SM]\\d.*(?=\\))")


temp = NR_T2vsT1_BTM_signi$pathway
temp2 <- stringr::str_extract(string = temp, pattern = "(?<=\\()[SM]\\d.*(?=\\))")

#######subset metaData file to only include T3 and T1 for each of TRs and NRs file######
TR_metaDataT3T1 = subset(metaFile, metaFile$Category == "TR" & metaFIle$SampleTime %in% c("T1","T3"))
NR_metaDataT3T1 = subset(metaFile, metaFile$Category == "NR" & metaFIle$SampleTime %in% c("T1","T3"))


############ cibersort subsets ############
ciberSortResults = read.csv("../newCiberSortResults.csv",header = TRUE)

TR_T3T1_NewcellProps = subset(ciberSortResults, rownames(ciberSortResults) %in% TR_metaDataT3T1$SampleName)
NR_T3T1_NewcellProps = subset(ciberSortResults, rownames(ciberSortResults) %in% NR_metaDataT3T1$SampleName)

TR_T3vsT1_newFinalDF = cbind(TR_T3T1_NewcellProps,TR_metaDataT3T1,TR_T3vsT1_BTM_MeanDF)
NR_T3vsT1_newFinalDF = cbind(NR_T3T1_NewcellProps,NR_metaDataT3T1,NR_T3vsT1_BTM_MeanDF)


######### running associations #########
assocs = data.frame(LMCompared = character(),CellType = character(),IndependantVariable = character(),estimate = double(),std.error = double(),statistic = double(),df = double(),p.value = double())
for (i in newCellTypesToAnalyze){
  print(i)
  for (j in 24:length(colnames(TR_T3vsT1_newFinalDF))) {
    pathwayName = names(TR_T3vsT1_newFinalDF[j])
    print(pathwayName)
    pred_vars = c(pathwayName,"Gender", "Age","(1 | patientID)")
    modToFit = reformulate(pred_vars, response = i)
    t = tidy(lmer(modToFit,data=TR_T3vsT1_newFinalDF)) %>%
      #filter(p.value < 0.05) %>%
      select(term,estimate,std.error,statistic,df,p.value) %>%
      as.data.frame()
    assocs[nrow(assocs)+1,] <- c(paste0(i,":",pathwayName),i,t[2,])
  }
}
TR_T3vsT1_BTM_NewcellPropAssoc = assocs

assocs = data.frame(LMCompared = character(),CellType = character(),IndependantVariable = character(),estimate = double(),std.error = double(),statistic = double(),df = double(),p.value = double())
for (i in newCellTypesToAnalyze){
  print(i)
  for (j in 24:length(colnames(NR_T3vsT1_newFinalDF))) {
    pathwayName = names(NR_T3vsT1_newFinalDF[j])
    print(pathwayName)
    pred_vars = c(pathwayName,"Gender", "Age","(1 | patientID)")
    modToFit = reformulate(pred_vars, response = i)
    t = tidy(lmer(modToFit,data=NR_T3vsT1_newFinalDF)) %>%
      #filter(p.value < 0.05) %>%
      select(term,estimate,std.error,statistic,df,p.value) %>%
      as.data.frame()
    assocs[nrow(assocs)+1,] <- c(paste0(i,":",pathwayName),i,t[2,])
  }
}
NR_T3vsT1_BTM_NewcellPropAssoc = assocs
