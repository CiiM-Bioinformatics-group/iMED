try(dev.off())
rm(list = ls())

library(dplyr)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(ggVennDiagram)
library(ComplexHeatmap)
library(reshape2)
library(ggsci)

load('output/data_flt.Rdata')
grps <- read.csv('meta.csv') %>% select(ProbandID, Category) %>% distinct(ProbandID, .keep_all=T)
merge(meta, grps, by = 'ProbandID', sort=F) %>% arrange(name, colnames(df)) -> meta

all(colnames(df) == meta$name)

meta$ab_B <- log2(meta$ab_B)
meta$ab_H1N1 <- log2(meta$ab_H1N1)
meta$ab_H3N2 <- log2(meta$ab_H3N2)
rownames(df) <- make.unique(annot$hmdb)

meta %<>% arrange(ProbandID, time)
df <- df[, meta$name]
rownames(df) <- annot$idx
# rownames(df) <- make.unique(annot$hmdb)


# Fit linear model for each time, strain separately.
# Model: metabolite - LFC + age + gender + error
res <- data.frame()

for (t in c('T1', 'T2', 'T3')) {
  for (strain in c('ab_B', 'ab_H3N2', 'ab_H1N1')) {

      meta.sub <- meta %>% filter(time == t)


      for (metabolite in rownames(df)) {

        age <- meta.sub$age
        gender <- meta.sub$gender
        met = df[meta.sub$name] %>% t() %>% as.data.frame() %>% pull(metabolite)
        lfc = meta %>% filter(time == 'T3') %>% pull(strain)

        lm(met ~ lfc + age + gender) %>%
          broom.mixed::tidy() %>%
          select(term, estimate, std.error, statistic, p.value) %>%
          mutate(strain = strain, time = t, metabolite = metabolite) %>%
          as.data.frame() -> r

        res <- rbind(res, r)

      }
  }
}


# Plots
# Rankplot per time and strain separately, colored by significance. Do we find sig. associations to each of the strains? at which timepoint?
# What are the metabolites that are associated to the logFC? HMDB classes
# Density plots of estimates per strain, colored by time: What is the best time to associate to the logFC?
# Rankproduct to produce the following vectors: H1N1, H3N2, B, H1N1+H3N2, H1N1+B, H3N2+B and triple. Viridis heatmap to show the ranks for each metabolite in each comparison
res %<>%
  filter(term == 'lfc') %>%
  mutate(sig = ifelse(p.value < 0.05, T, F)) %>%
  group_by(strain, time) %>%
  arrange(p.value) %>%
  mutate(rank = seq_along(p.value))

write.csv('output/assoc_metabolites_lfc.csv', x = res)

rm(list = ls())
res <- read.csv('output/assoc_metabolites_lfc.csv')

res <- merge(res, annot, by.x = 'metabolite', by.y = 'hmdb', all.x=T, all.y=F, sort=F)
# Rankplot
pdf('output/rankplot_assoc_metabolome_lfc.pdf', width = 10, height = 5)
ggplot(res) +
  geom_point(aes(x = rank, y =  -log10(p.value), color = sig), size=.8) +
  facet_grid(time ~ strain) +
  scale_color_manual(values = c('gray', 'red')) +
  theme_classic() +
  theme(strip.text = element_text(face = 'bold'), aspect.ratio = 1) +
  labs(y = '-log10(Pvalue)', x = 'Rank', color = 'Significance') +
  geom_hline(yintercept = 0, lty=2, col = 'black') +
  scale_x_continuous(breaks = seq(0, 800, 400)) +
  scale_y_continuous(breaks = c(0, 2, 4)) +
  ggrepel::geom_text_repel(data = res %>% group_by(time, strain) %>% slice_min(order_by = p.value, n = 3),
                           aes(x = rank, y =  -log10(p.value), label = ids), size=2, min.segment.length = 0)
dev.off()

# Overlap between the sig ones?
# Between times within strain
doVenn <- function(x) {

  venn <- Venn(x)
  data <- process_data(venn)
  # data@region$perc <- ((data@region$count / sum(data@region$count)) * 100) %>% round() %>% paste0(., '%')
  # data@region$label <- paste0(data@region$count, '\n', data@region$perc)

  return(ggplot() +
          geom_sf(aes(fill = count), data = venn_region(data), show.legend = F) +
          geom_sf(data = venn_setedge(data)) +
          geom_sf_text(aes(label = name), data = venn_setlabel(data), nudge_y = .05) +
          geom_sf_label(aes(label = count), data = venn_region(data), label.size = 0, alpha = 0) +
          theme_void() +
          scale_fill_gradientn(colours = c('white', 'red')) +
          # labs(title = title) +
          theme(plot.title = element_text(hjust = 0.5))
         )
}

# between strains within time
x <- list("B" = res %>% filter(sig) %>% filter(strain == 'ab_B') %>% filter(time == 'T1') %>% pull(metabolite),
          "H1N1" = res %>% filter(sig) %>% filter(strain == 'ab_H1N1') %>% filter(time == 'T1') %>% pull(metabolite),
          "H3N2" = res %>% filter(sig) %>% filter(strain == 'ab_H3N2') %>% filter(time == 'T1') %>% pull(metabolite))
pdf('output/overlap_logFC_assoc_met_T1.pdf', width = 3, height = 3)
doVenn(x)
dev.off()

x <- list("B" = res %>% filter(sig) %>% filter(strain == 'ab_B') %>% filter(time == 'T2') %>% pull(metabolite),
          "H1N1" = res %>% filter(sig) %>% filter(strain == 'ab_H1N1') %>% filter(time == 'T2') %>% pull(metabolite),
          "H3N2" = res %>% filter(sig) %>% filter(strain == 'ab_H3N2') %>% filter(time == 'T2') %>% pull(metabolite))
pdf('output/overlap_logFC_assoc_met_T2.pdf', width = 3, height = 3)
doVenn(x)
dev.off()

x <- list("B" = res %>% filter(sig) %>% filter(strain == 'ab_B') %>% filter(time == 'T3') %>% pull(metabolite),
          "H1N1" = res %>% filter(sig) %>% filter(strain == 'ab_H1N1') %>% filter(time == 'T3') %>% pull(metabolite),
          "H3N2" = res %>% filter(sig) %>% filter(strain == 'ab_H3N2') %>% filter(time == 'T3') %>% pull(metabolite))
pdf('output/overlap_logFC_assoc_met_T3.pdf', width = 3, height = 3)
doVenn(x)
dev.off()

res[which(res$time == 'T3'), ]$time <- 'T4'
res[which(res$time == 'T2'), ]$time <- 'T3'


pdf('output/densities_estimate_assoclfc.pdf', width = 8, height = 2)
ggplot(res) +
  geom_density(aes(x = estimate, color = time)) +
  facet_wrap('strain') +
  theme_classic() +
  scale_color_lancet()
dev.off()

# Heatmap of test statistics
res %>% filter(sig) %>% pull(metabolite) %>% unique() -> mets

res %>%
  filter(metabolite %in% mets) %>%
  filter(time == 'T2') %>%
  dcast(data = ., formula = strain ~ metabolite, value.var = 'estimate') %>% tibble::column_to_rownames('strain')-> mat


png('output/heatmap_assoc_lfc_bestimates.png', width = 10, height = 2, res = 1200, units = 'in')
Heatmap(mat, cluster_rows = F, show_column_names = F, row_names_side = 'left', name = 'B estimate')
dev.off()


# show the ones that are neg. correlated at T3
res %>% filter(time == 'T2') %>% arrange(estimate) %>% head(15)
head(res)
meta.sub <- meta %>% filter(time == 'T2')
m <- 'HMDB0135284'
met = df[meta.sub$name] %>% t() %>% as.data.frame() %>% pull(m)

pdf('output/corr_met_lfc.pdf', width = 8, height = 3)
ggpubr::ggarrange(
  ggplot() +
    geom_point(aes(x = meta %>% filter(time == 'T3') %>% pull('ab_B'), y = met)) +
    labs(x = 'LogFC', y = 'Metabolite', title = 'B') +
    theme_classic() +
    theme(plot.title = element_text(hjust = .5)) +
    geom_smooth(method = 'lm', aes(meta %>% filter(time == 'T3') %>% pull('ab_B'), y = met)) +
    annotate(geom = 'text', x = 3, y = 6, label = 'B: -0.09\nP: 0.003'),

  ggplot() +
    geom_point(aes(x = meta %>% filter(time == 'T3') %>% pull('ab_H1N1'), y = met)) +
    labs(x = 'LogFC', y = 'Metabolite', title = 'H1N1') +
    theme_classic() +
    theme(plot.title = element_text(hjust = .5)) +
    geom_smooth(method = 'lm', aes(meta %>% filter(time == 'T3') %>% pull('ab_H1N1'), y = met)) +
    annotate(geom = 'text', x = 3, y = 6, label = 'B: -0.05\nP: 0.06'),

  ggplot() +
    geom_point(aes(x = meta %>% filter(time == 'T3') %>% pull('ab_H3N2'), y = met)) +
    labs(x = 'LogFC', y = 'Metabolite', title = 'H3N2') +
    theme_classic() +
    theme(plot.title = element_text(hjust = .5)) +
    geom_smooth(method = 'lm', aes(meta %>% filter(time == 'T3') %>% pull('ab_H3N2'), y = met)) +
    annotate(geom = 'text', x = 3, y = 6, label = 'B: -0.07\nP: 0.015'), nrow = 1

)
dev.off()
cor.test(meta %>% filter(time == 'T3') %>% pull('ab_B'), y = met)
cor.test(meta %>% filter(time == 'T3') %>% pull('ab_H3N2'), y = met)
cor.test(meta %>% filter(time == 'T3') %>% pull('ab_H1N1'), y = met)

for (metabolite in rownames(df)) {

  age <- meta.sub$age
  gender <- meta.sub$gender
  met = df[meta.sub$name] %>% t() %>% as.data.frame() %>% pull(metabolite)
  lfc = meta %>% filter(time == 'T3') %>% pull(strain)

  lm(met ~ lfc + age + gender) %>%
    broom.mixed::tidy() %>%
    select(term, estimate, std.error, statistic, p.value) %>%
    mutate(strain = strain, time = t, metabolite = metabolite) %>%
    as.data.frame() -> r
