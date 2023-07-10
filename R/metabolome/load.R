try(dev.off())
rm(list = ls())

library(dplyr)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggsci)

# Meta metabolites
ions <- read.xlsx(xlsxFile = 'DATA_CURATED_reformatted.xlsx', sheet = 'ions')
annot <- data.frame(ids = ions$Top.annotation.name)

# Split into three cols: KEGG (CXXXX), CHEBI IDs (CHEBI: XXXX) and HMDB: HMDBXXXXX.
# Everything is separated by ;
spl <- strsplit(x = ions$Top.annotation.ids, split = ';')
getID <- function(x, pattern) { ind <- grep(pattern = pattern, x = x)[1]; return(x[ind]) }

# Chebi
chebi <- lapply(spl, FUN = getID, pattern = 'CHEBI')
chebi <- lapply(chebi, FUN = function(x) { strsplit(x = x, split = ':')[[1]][2] -> num; return(as.numeric(num))} )
annot$chebi <- unlist(chebi)

# HMBD
hmdb <- lapply(spl, FUN = getID, pattern = 'HMDB')
hmdb <- gsub(pattern = ' ', replacement = '', x = hmdb)
annot$hmdb <- hmdb

# KEGG
kegg <- lapply(spl, FUN = getID, pattern = 'C')
kegg <- lapply(kegg, FUN = function(x) {x <- ifelse(grepl(':', x), NA, x); return(x)} )
annot$kegg <- unlist(kegg)
annot$idx <- rownames(annot)

# Annotation of the metabolites
annot2 <- read.xlsx(xlsxFile = 'DATA_CURATED_reformatted.xlsx',
                    sheet = 'annotation') %>%
  filter(!is.na(ionIdx))

# Filters based on Bryan's email
# mzDelta
annot2 %<>% filter(abs(mzDelta) < 0.001)
# Remove compounds without HMDB annotation
annot %<>% filter(!is.na(hmdb))

keep <- intersect(annot2$ionIdx, annot$idx)


annot %<>% filter(idx %in% keep) %>% arrange(match(idx, keep))
annot2 %<>% filter(ionIdx %in% keep) %>% arrange(match(ionIdx, keep))

df <- read.xlsx(xlsxFile = 'DATA_CURATED_reformatted.xlsx', sheet = 'ions')
df %<>% filter(ionIdx %in% keep) %>% arrange(match(ionIdx, keep)) %>% tibble::column_to_rownames('ionIdx')
df %<>% select(-all_of(c("ionMz", "ion.Average.Intensity", "Top.annotation.name", "Top.annotation.ids", "Top.annotation.formula", "Top.annotation.ion", "Top.annotation.score")))

# Meta samples
meta <- read.csv('meta.csv', row.names=1, header=T)
df <- df[, which(colnames(df) %in% meta$name)]
meta %<>% arrange(match(name, colnames(df)))
stopifnot(all(meta$name == colnames(df)))


# Log-2 normalisation
df <- log2(df)

# Add the novel HMDB names
hmdbnames <- data.table::fread('metabolite_HMDB_names_small.csv') %>%
  as.data.frame() %>%
  set_colnames(c('hmdb', 'hmdb_name'))

table(hmdbnames$hmdb %in% annot$hmdb)
merge(annot, hmdbnames, by = 'hmdb', sort=F, all.x=T, all.y=F) -> annot
annot %<>% arrange(match(idx, rownames(df)))
stopifnot(all(annot$idx == rownames(df)))

annot$finalname <- ifelse(is.na(annot$hmdb_name), annot$ids, annot$hmdb_name)
rm(chebi, kegg, hmdb, keep, getID, ions, spl)

save.image('output/data_flt.Rdata')
