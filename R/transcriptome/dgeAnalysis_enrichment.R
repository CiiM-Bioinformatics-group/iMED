library(tximport)
#library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(edgeR)
library(limma)
library(fgsea)


########transcript to gene file readily available for either ENSEMBL or NCBI
new_tx2gene = read.table("your file", header = T,na.strings = TRUE,fill = TRUE) ### human ###
###### location of all the salmon count files per transcript ####
dataFiles="salmonFilesLocation"

####metadataFile###########
f1 = read.csv("metaFile.csv",header = T,sep = ";")
f1$fileName = paste0(f1$SampleID,"-quant.sf")

Names_count = data.frame(patientID = as.factor(f1$iMedID),SampleName = f1$SampleID,Gender = as.factor(f1$Sex), Time = as.factor(f1$SampleCollect), Category = as.factor(f1$Category), Age = f1$Age)
#Names_count

txi.salmon.limma <- tximport(paste0(dataFiles,f1$fileName), type = "salmon", tx2gene = new_tx2gene, countsFromAbundance = "lengthScaledTPM")
colnames(txi.salmon.limma$counts) = Names_count$SampleName

rawCounts = txi.salmon.limma$counts
write.table(rawCounts, file = "rawCounts_all_entriesWithGeneNamesOnly.txt",quote = F)
A <- rowSums(rawCounts)
isexpr <- A > 10
table(isexpr)
genesTouse = rawCounts[isexpr,]

y = DGEList(genesTouse)
y <- calcNormFactors(y)
cpmCounts = cpm(y)



Names_count$CategoryTime = paste0(Names_count$Category,"_",Names_count$SampleTime)

Names_count$Category = factor(Names_count$Category,levels = c("NR","TR"))

###########TR vs NR for all time points #################

mod = model.matrix(~0+Category+Gender+SampleTime+Age,data=Names_count)

vun = voom(y,mod, plot=TRUE,normalize = "quantile")
corfit = duplicateCorrelation(vun,mod,block=Names_count$patientID)
corfit$consensus.correlation
vun = voom(y,mod,normalize.method = "quantile",plot = TRUE,block=Names_count$patientID,correlation=corfit$consensus.correlation)
fit = lmFit(vun,mod,block = Names_count$patientID,correlation = corfit$consensus.correlation)
fit = eBayes(fit)
summary(decideTests(fit))
cont.matrix = cbind("TRvsNR" =c(-1,1,0,0,0,0,0,0,0,0,0,0))
fit2 = contrasts.fit(fit,cont.matrix)
fit2 = eBayes(fit2)
summary(decideTests(fit2))

TRvsNR = topTable(fit2,coef = 1, number = Inf)

#########TR vs NR for each timepoint and time dynamics for each category #######
mod_catTime = model.matrix(~0+ CategoryTime + Gender + Age,data=Names_count)
vunCatTime = voom(y,mod_catTime, plot=TRUE,normalize = "quantile")

corfit2 = duplicateCorrelation(vunCatTime,mod_catTime,block=Names_count$patientID)
corfit2$consensus.correlation
vunCatTime = voom(y,mod_catTime,normalize.method = "quantile",plot = TRUE,block=Names_count$patientID,correlation=corfit2$consensus.correlation)
fitTime = lmFit(vunCatTime,mod_catTime,block = Names_count$patientID,correlation = corfit2$consensus.correlation)
fitTime = eBayes(fitTime)
summary(decideTests(fitTime))

cont.matrix1 = makeContrasts(T1_TRvsNR = (CategoryTimeTR_T1 - CategoryTimeNR_T1),
                             T2_TRvsNR = (CategoryTimeTR_T2 - CategoryTimeNR_T2),
                             T3_TRvsNR = (CategoryTimeTR_T3 - CategoryTimeNR_T3),
                             T4_TRvsNR = (CategoryTimeTR_T4 - CategoryTimeNR_T4),
                             T5_TRvsNR = (CategoryTimeTR_T5 - CategoryTimeNR_T5),levels = mod_catTime)

cont.matrix2 = makeContrasts(TR_T2vsT1 = CategoryTimeTR_T2 - CategoryTimeTR_T1,
                             TR_T3vsT1 = CategoryTimeTR_T3 -CategoryTimeTR_T1,
                             TR_T4vsT1 = CategoryTimeTR_T4 - CategoryTimeTR_T1, 
                             TR_T5vsT1 = CategoryTimeTR_T5 - CategoryTimeTR_T1,
                             NR_T2vsT1 = CategoryTimeNR_T2 - CategoryTimeNR_T1,
                             NR_T3vsT1 = CategoryTimeNR_T3 - CategoryTimeNR_T1,
                             NR_T4vsT1 = CategoryTimeNR_T4 - CategoryTimeNR_T1,
                             NR_T5vsT1 = CategoryTimeNR_T5 - CategoryTimeNR_T1 ,levels = mod_catTime)

fitTime2 = contrasts.fit(fitTime,cont.matrix1)
fitTime2 = eBayes(fitTime2)
summary(decideTests(fitTime2))

timeT1_TRvsNR = topTable(fitTime2,coef = 1,number = Inf)
timeT2_TRvsNR = topTable(fitTime2,coef = 2,number = Inf)
timeT3_TRvsNR = topTable(fitTime2,coef = 3,number = Inf)
timeT4_TRvsNR = topTable(fitTime2,coef = 4,number = Inf)
timeT5_TRvsNR = topTable(fitTime2,coef = 5,number = Inf)


fitTime3 = contrasts.fit(fitTime,cont.matrix2)
fitTime3 = eBayes(fitTime3)
summary(decideTests(fitTime3))

TRs_T2vsT1 = toptable(fitTime3,coef = 1, number = Inf)
TRs_T3vsT1 = toptable(fitTime3,coef = 2, number = Inf)
TRs_T4vsT1 = toptable(fitTime3,coef = 3, number = Inf)
TRs_T5vsT1 = toptable(fitTime3,coef = 4, number = Inf)
NRs_T2vsT1 = toptable(fitTime3,coef = 5, number = Inf)
NRs_T3vsT1 = toptable(fitTime3,coef = 6, number = Inf)
NRs_T4vsT1 = toptable(fitTime3,coef = 7, number = Inf)
NRs_T5vsT1 = toptable(fitTime3,coef = 8, number = Inf)

#########enrichment analysis ################
ranked_list <- TRvsNR$t
names(ranked_list) <- rownames(TRvsNR)

fgsea::gmtPathways("path_To_BTMs_gmt") -> genesets_BTM

ranked_list <- ranked_list[!duplicated(names(ranked_list))]
sort(ranked_list, decreasing = T) -> ranked_list

fgseaRes <- fgsea(pathways = genesets_BTM, 
                  stats    = ranked_list, nPermSimple = 5000)


ranked_list_TRsT2vsT1<- TRs_T2vsT1$t
names(ranked_list_TRsT2vsT1) <- rownames(TRs_T2vsT1)


fgseaRes_TRsT2vsT1 <- fgsea(pathways = genesets_BTM, 
                  stats    = ranked_list_TRsT2vsT1, nPermSimple = 5000)


ranked_list_TRsT3vsT1<- TRs_T3vsT1$t
names(ranked_list_TRsT3vsT1) <- rownames(TRs_T3vsT1)


fgseaRes_TRsT3vsT1 <- fgsea(pathways = genesets_BTM, 
                            stats    = ranked_list_TRsT3vsT1, nPermSimple = 5000)

ranked_list_NRsT2vsT1<- NRs_T2vsT1$t
names(ranked_list_NRsT2vsT1) <- rownames(NRs_T2vsT1)


fgseaRes_NRsT2vsT1 <- fgsea(pathways = genesets_BTM, 
                            stats    = ranked_list_NRsT2vsT1, nPermSimple = 5000)


ranked_list_NRsT3vsT1<- NRs_T3vsT1$t
names(ranked_list_NRsT3vsT1) <- rownames(NRs_T3vsT1)


fgseaRes_NRsT3vsT1 <- fgsea(pathways = genesets_BTM, 
                            stats    = ranked_list_NRsT3vsT1, nPermSimple = 5000)

