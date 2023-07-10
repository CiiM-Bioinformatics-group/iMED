try(dev.off())
rm(list = ls())

library(dplyr)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(reshape2)
library(limma)
library(ggsci)
library(patchwork)

load('output/data_flt.Rdata')
rm(annot, annot2, df, hmdbnames)
annotation <- read.csv('metabolite_HMDB_taxonomy.csv', sep = '\t') %>% tibble::rownames_to_column('HMDB_ID')

df <- data.table::fread('endgenousMetabolites.txt')
load('20220902_metabolites_noDrug_NN.RData')

annot <- merge(df, annot_metabolites_noDrug, by.x = 'CHEMICAL_FORMULA', by.y = 'Formula')
annot <- merge(annot, annotation, by = 'HMDB_ID', all.y=F, all.x = T)

dat_metabolites_noDrug %<>% tibble::rownames_to_column('ionIdx') %>% filter(ionIdx %in% annot$ionIdx) %>% tibble::column_to_rownames('ionIdx')

df <- dat_metabolites_noDrug

annot %>% count(ionIdx, super_class, class, sub_class) %>% group_by(ionIdx) %>% arrange(n) %>% summarise(super_class = super_class[1], class = class[1], sub_class = sub_class[1]) -> annotfinal

grps <- read.csv('meta.csv') %>% select(ProbandID, Category) %>% distinct(ProbandID, .keep_all=T)
merge(meta, grps, by = 'ProbandID', sort=F) %>% arrange(name, colnames(df)) -> meta
meta$Category2 <- ifelse(meta$Category == 'TR', 'TR', ifelse(meta$Category == 'NR', 'NR', 'Other'))

all(colnames(df) == meta$name)
all(rownames(df) == annotfinal$ionIdx)

cond_time <- factor(paste0(meta$Category2, '_', meta$time))
gender <- factor(meta$gender, levels = c('male', 'female'))
age <- meta$age
probands <- factor(meta$ProbandID)

design <- model.matrix(~ 0 + cond_time + age + gender)
dupcor <- duplicateCorrelation(object = df, design = design, block = probands)
fit <- lmFit(as.matrix(df), design = design, block = probands, correlation = dupcor$consensus.correlation)

fit <- eBayes(fit)
DE_time <- list()



for (t in c('T2', 'T3')) {
  for (ct in unique(meta$Category2)) {

    # Each timepoint vs baseline per category
    cnt <- paste0('cond_time', ct, '_', t, '- cond_time', ct, '_T1')
    outname <- paste0(t, 'vsT1_', ct)

    contrasts.fit(fit, makeContrasts(cnt, levels = design)) %>% eBayes() %>%
      topTable(sort.by = 'none', adjust.method = 'BH', number = Inf) %>%
      mutate(ionIdx = rownames(.),
             significance = ifelse(adj.P.Val < 0.05, T, F),
             direction = ifelse(logFC < 0.0, 'Downregulated', 'Upregulated')) %>%
      arrange(adj.P.Val) %>%
      set_rownames(NULL) %>%
      left_join(x = ., y = annotfinal %>% mutate(ionIdx = as.character(ionIdx)), by = 'ionIdx')  -> DE_time[[outname]]

  }
}

write.xlsx(x = DE_time, file = 'output/DE_classes_endogeneous_metabolites.xlsx', overwrite=T)



DE_classes <- list()
for (t in c('T1', 'T2', 'T3')) {
  for (ct in c('NR', 'Other')) {

    # Each timepoint vs baseline per category
    cnt <- paste0('cond_time', ct, '_', t, '- cond_timeTR_', t)
    outname <- paste0(ct, 'vsTR_', t)
    paste0()

    contrasts.fit(fit, makeContrasts(cnt, levels = design)) %>% eBayes() %>%
      topTable(sort.by = 'none', adjust.method = 'BH', number = Inf) %>%
      mutate(metabolite = annot$ids,
             hmdb = annot$hmdb,
             kegg = annot$kegg,
             chebi = annot$chebi,
             significance = ifelse(adj.P.Val < 0.05, T, F),
             direction = ifelse(logFC < 0.0, 'Downregulated', 'Upregulated')) %>%
      arrange(adj.P.Val) %>%
      set_rownames(NULL) -> DE_classes[[outname]]



  }
}
DE_classes
lapply(DE_classes, function(x) {print(table(x$significance))})



# Nr of DEG
# Over time, split for up and down
rbind(
  lapply(DE_time, function(x) {x %>% filter(significance) %>% filter(direction == 'Upregulated') %>% nrow()}) %>% do.call(rbind, .) %>%
    as.data.frame() %>% set_colnames('deg') %>% cbind(., colsplit(rownames(.), pattern = '_', names = c('time', 'resp'))) %>% mutate('direction' = 'up'),

  lapply(DE_time, function(x) {x %>% filter(significance) %>% filter(direction == 'Downregulated') %>% nrow()}) %>% do.call(rbind, .) %>%
    as.data.frame() %>% set_colnames('deg') %>% cbind(., colsplit(rownames(.), pattern = '_', names = c('time', 'resp'))) %>% mutate('direction' = 'down')
) -> time

time <- rbind(time, data.frame(deg = 0, time = rep('T1', 3), resp = c('TR', 'Other', 'NR', 'TR', 'Other', 'NR'), direction = c(rep('up', 3), rep('down', 3))))

time$y <- (time$deg / nrow(df)) * 100

time$y <- ifelse(time$direction == 'down', -1 * time$y, time$y)
time$grp <- paste0(time$resp, '_', time$direction)
time$resp <- factor(time$resp, levels = rev(c('NR', 'Other', 'TR')))

time$time[which(time$time == 'T2vsT1')] <- 'T3'
time$time[which(time$time == 'T3vsT1')] <- 'T4'

time %<>% filter(resp != 'Other')


time$time <- dplyr::recode(time$time, 'T1' = 'Day 0', 'T3' = 'Day 7', 'T4' = 'Day 35')
time$time <- factor(time$time, levels = c('Day 0', 'Day 7', 'Day 35'))


pdf('output/DE_subsetted_metabolites_overtime.pdf', width = 2, height = 2)
ggplot(time) +
  geom_hline(yintercept = 0, lty=2, color = 'darkgray') +
  geom_point(aes(x = time, y = y, color = resp)) +
  geom_line(aes(x = time, y= y, group = grp, color = resp)) +
  labs(y = 'Percentage DAM', color = 'Category', title = 'Metabolites', x = 'Time') +
  # scale_color_manual(values = c('indianred2', 'brown', 'cyan4','dodgerblue')) +
  scale_color_manual(values = c('TR' = 'indianred2', 'other' = 'brown',"NR" = 'dodgerblue')) +
  theme_bw() +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = .5, size=8),
        panel.grid = element_blank(),
        axis.ticks = element_line(colour = 'darkgrey'),
        panel.border = element_rect(colour = 'darkgrey'),
        axis.title = element_text(size=8),
        axis.text = element_text(size= 6)) +
  scale_y_continuous(labels=function(x) {paste0(abs(x), '%')} ) +
  guides(color = 'none')
  # annotate(geom = 'text', x = 0.6, y = 2, label = 'Up', size=3, fontface = 'bold') +
  # annotate(geom = 'text', x = 0.7, y = -2, label = 'Down', size=3, fontface = 'bold')
dev.off()





# Barplot categories
table(DE_time$T2vsT1_TR$super_class, useNA = 'ifany') %>% as.data.frame() %>%
  set_colnames(c('Class', 'Freq')) %>% mutate(Proportion = Freq / sum(Freq), time = 'Baseline') -> baseline

deg <- DE_time

lapply(deg, function(x) {x %>% filter(significance) %>% filter(direction == 'Upregulated') %>%
    select(super_class) %>% table(useNA = 'ifany') %>% as.data.frame()}) -> up
up$T2vsT1_NR <- NULL
for (i in names(up)) {up[[i]]$cohort_time <- i; up[[i]]$prop <- up[[i]]$Freq / sum(up[[i]]$Freq)}
do.call(rbind.data.frame, up) %>% set_rownames(NULL) %>% set_colnames(c('Class', 'Freq', 'time_class', 'Proportion')) %>%  cbind(., colsplit(.$time_class, pattern = '_', names = c('time', 'class'))) %>% select(-time_class) -> up

lapply(deg, function(x) {x %>% filter(significance) %>% filter(direction == 'Downregulated') %>%
    select(super_class) %>% table(useNA = 'ifany') %>% as.data.frame()}) -> down
for (i in names(down)) {down[[i]]$cohort_time <- i; down[[i]]$prop <- down[[i]]$Freq / sum(down[[i]]$Freq)}
do.call(rbind.data.frame, down) %>% set_rownames(NULL) %>% set_colnames(c('Class', 'Freq', 'time_class', 'Proportion')) %>%  cbind(., colsplit(.$time_class, pattern = '_', names = c('time', 'class'))) %>% select(-time_class) -> down


########
aspect.ratio = 3.5
# Measured plot only
baseline$label <- ifelse(baseline$Proportion > 0.05, paste0(round(baseline$Proportion * 100), '% (', baseline$Freq, ')'), NA)
baseline$x <- 'Measured'

m <- ggplot(baseline, aes(x = 'Measured', y = Proportion)) +
  geom_col(aes(fill = Class)) +
  # theme_classic() +
  theme_bw() +
  # theme(plot.title = element_text(hjust = .5, size=8),
  #       panel.grid = element_blank(), legend.text = element_text(size = 6), legend.title = element_text(size = 8),
  #       strip.text = element_text(size = 8),
  #       axis.title = element_text(size=8), axis.text = element_text(size= 6))  +
  theme(plot.title = element_text(hjust = .5, size=8),
        panel.grid = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        strip.text = element_text(size = 8, face = 'bold'),
        axis.ticks = element_line(colour = 'darkgrey'),
        panel.border = element_rect(colour = 'darkgrey'),
        plot.margin = margin(0, 0, 0, 0, 'pt'),

        # strip.background.y = element_rect(fill = 'white'),
        strip.background = element_blank(),
        # panel.border = element_blank(),
        # panel.spacing = unit(0.1, "lines"),
        # panel.border = element_rect(fill = NA, colour = "black"),
        # axis.line.y = element_line(colour = 'darkgrey'),
        # axis.line.x = element_line(colour = 'darkgrey'),
        # plot.margin = unit(c(0,0,0,0), 'cm'),
        # panel.border = element_blank(), axis.line.x = element_line(),
        axis.title = element_text(size=8),
        axis.title.x = element_blank(),
        legend.position = 'bottom',
        aspect.ratio = 3.5,
        axis.text = element_text(size= 6)) +
  scale_fill_manual(values = ggpubr::get_palette('npg', 9)) +
  geom_text(aes(label = label, group=Class), position = position_stack(vjust = .5), size = 0.36 * 4) +
  scale_y_continuous(label = function(x) {paste0(x * 100, '%')}) +
  guides(fill = guide_legend(ncol = 2)) +
  labs(y = 'Percentage')

m

pdf('output/legend.pdf', width = 4, height = 4)
cowplot::get_legend(m) %>% plot()
dev.off()


#### Up and down
df <- rbind(
  # baseline %>% mutate(class = 'Baseline', direction = 'Upregulated'),
  # baseline %>% mutate(class = 'Baseline', direction = 'Downregulated'),
  up %>% mutate(direction = 'Upregulated'),
  down %>% mutate(direction = 'Downregulated')
)

df$label <- ifelse(df$Proportion > 0.1, paste0(round(df$Proportion * 100), '% (', df$Freq, ')'), NA)
# df$Proportion <- ifelse(df$direction == 'Upregulated', df$Proportion, -1 * df$Proportion)

# df$time <- dplyr::recode(df$time, 'Baseline' = 'Measured', 'T2vsT1' = 'Day 7', 'T3vsT1' = 'Day 35')
df$time <- dplyr::recode(df$time, 'T2vsT1' = 'Day 7', 'T3vsT1' = 'Day 35')
# df$class <- dplyr::recode(df$class, 'Baseline' = 'Measured')
# df$class <- factor(df$class, levels = c('Measured', 'NR', 'TR'))
df$class <- factor(df$class, levels = c('NR', 'TR'))
df$time <- factor(df$time, levels = c('Day 7', 'Day 35'))
# df$time <- factor(df$time, levels = c('Measured', 'Day 7', 'Day 35'))
df %<>% filter(class != 'Other')
df$direction <- factor(df$direction, levels = c('Upregulated', 'Downregulated'))



# pdf('output/barplot_proportions_subsetted_metabolites.pdf', width = 8, height = 4)
trnr <- ggplot(df, aes(x = time, y = Proportion, fill = Class)) +
  geom_bar(position = 'fill', stat = 'identity') +
  facet_grid(direction ~ class, scales = 'free', drop=T) +
  # geom_hline(yintercept = 0, lty=2, color = 'black') +
  geom_text(aes(label = label), position = position_fill(vjust = .5), size = 0.36 * 4) +
  theme_bw() +
  # theme_classic() +
  scale_fill_manual(values = ggpubr::get_palette('npg', 9)) +
  labs(y = 'Proportion', x = 'Time') +
  theme(plot.title = element_text(hjust = .5, size=8),
        panel.grid = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        axis.ticks = element_line(colour = 'darkgrey'),
        # panel.border = element_rect(colour = 'darkgrey', fill=NA),
        strip.background = element_rect(fill = 'white', colour = 'darkgrey', linewidth = .5),
        # strip.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        # strip.background = element_rect(colour = 'darkgrey'),
        panel.border = element_blank(),
        axis.line.y = element_line(colour = 'darkgrey', linewidth = .25),
        axis.line.x = element_line(colour = 'darkgrey', linewidth = .25),
        legend.position = 'none',
        aspect.ratio = 1.5,
        plot.margin = margin(0, 0, 0, 0, 'pt'),
        axis.title = element_text(size=8),
        axis.text = element_text(size= 6),
        axis.title.y = element_blank()) +
  # scale_y_continuous(labels = function(x) {abs(x)}) +
  scale_y_continuous(label = function(x) {paste0(x * 100, '%')}) +
  guides(color='none')

trnr

pdf('output/newbarplot.pdf', width = 4, height = 3)
(m + plot_layout(guides = "collect", widths = c(1, 3)) & theme(legend.position = "none")) + (trnr + theme(legend.position = "none"))
dev.off()
