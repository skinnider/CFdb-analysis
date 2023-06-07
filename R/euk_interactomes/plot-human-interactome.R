setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(flavin)
library(PrInCE)
library(VennDiagram)
source("R/theme.R")

## PR curve, v1 vs. v2
# read networks
v1 = read.delim("data/euk_interactomes/CF-MS-interactome.tsv.gz") %>% 
  mutate(version = 'Version 1') ## from Zenodo
v2 = readRDS("data/euk_interactomes/networks/network-species=Homo_sapiens-n_replicates=166-version=V2-n_features=1.rds") %>% 
  mutate(version = 'Version 2')
# combine
pr = bind_rows(v1, v2) %>% 
  group_by(version) %>% 
  mutate(idx = row_number()) %>% 
  ungroup()

# label
pct = (nrow(v2) - nrow(v1)) / nrow(v1)
lab = paste0('+', round(100 * pct, digits = 0),
             '% increase\nin human interactions')

# plot
version_pal = c('grey88', pals::stepped()[3])
p1 = pr %>%
  ggplot(aes(x = idx, y = precision, color = version)) +
  geom_path(linewidth = 0.3) +
  geom_label(data = head(pr, 1) %>% mutate(label = lab),
             aes(label = lab, x = Inf, y = Inf), hjust = 1, vjust = 1, 
             size = 1.75, color = 'black', label.size = NA, lineheight = 0.95,
             label.padding = unit(0.7, 'lines'), fill = NA) +
  scale_y_continuous('Precision') +
  scale_x_continuous('Interactions (thousands)', labels = function(x) x / 1e3) +
  scale_color_manual('', values = version_pal) +
  guides(color = guide_legend(override.aes = list(linewidth = 0.5))) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.position = 'top',
        legend.key.size = unit(0.4, 'lines'))
p1
ggsave("fig/euk_interactomes/human-PR-v1-v2.pdf", p1, 
       width = 4.5, height = 4.5, units = "cm", useDingbats = FALSE)

## PR curve, v2 vs. Y2H/AP-MS 
# re-read networks
v1 = read.delim("data/euk_interactomes/CF-MS-interactome.tsv.gz") %>% 
  mutate(version = 'Version 1')
v2 = readRDS("data/euk_interactomes/human-interactome-v2-PR.rds") %>% 
  mutate(version = 'Version 2')
# re-label and recalculate precision
complexes = read.delim("data/complex/CORUM/complexes_human.txt") %>%
  as_annotation_list('gene_name', 'complex')
labels = complexes %>%
  adjacency_matrix_from_list() %>%
  make_labels(v2)
precision = calculate_precision(labels)
v2$label = labels
v2$precision = precision
# combine
pr = bind_rows(v1, v2) %>% 
  group_by(version) %>% 
  mutate(idx = row_number()) %>% 
  ungroup()

# plot RF against human high-throughput interactomes
hts_files = list.files("data/high_throughput", full.names = TRUE, 
                       pattern = '\\.gz$')
hts_ppis = map(hts_files, read.delim) %>%
  setNames(gsub("\\..*$", "", basename(hts_files)))
# calculate precision and recall
complexes = read.delim("data/complex/CORUM/complexes_human.txt") %>%
  flavin::as_annotation_list('gene_name', 'complex')
adj = adjacency_matrix_from_list(complexes)
precision = map_dbl(hts_ppis, ~ {
  make_labels(adj, .) %>% 
    mean(na.rm = TRUE)
})
recall = map_int(hts_ppis, nrow)
n_proteins = map_int(hts_ppis, ~ with(., n_distinct(c(gene_A, gene_B))))
hts = data.frame(interactome = names(hts_ppis),
                 precision = precision,
                 recall = recall,
                 n_proteins = n_proteins) %>%
  mutate(interactome = fct_recode(interactome,
                                  'BioPlex' = 'Huttlin2015',
                                  'BioPlex 2' = 'Huttlin2017',
                                  'BioPlex 3\nHCT116' = 'Huttlin2021-HCT116',
                                  'BioPlex 3\n293T' = 'Huttlin2021-293T',
                                  'BioPlex 2' = 'Huttlin2017',
                                  'QUBIC' = 'Hein2015',
                                  'HI-II-14' = 'Rolland2014',
                                  'HuRI' = 'Luck2019'))

pal = pals::stepped3()[seq(1, 19, 4)] %>% c('pink', 'black') %>% setNames(hts$interactome)
p2 = v2 %>%
  mutate(idx = row_number()) %>% 
  filter(idx %% 14 == 0) %>%
  ggplot(aes(x = idx, y = precision)) +
  geom_path(linetype = 'dotted', linewidth = 0.3) +
  # geom_path(linetype = 'dotted', size = 0.3) +
  geom_point(data = hts, aes(x = recall, color = interactome), size = 0.9) +
  geom_label_repel(data = hts, 
                   aes(x = recall, color = interactome, label = interactome),
                   size = 2.25, label.padding = 0.25, label.size = NA,
                   fill = NA, min.segment.length = 0, segment.size = 0.4,
                   lineheight = 0.85) +
  scale_y_continuous('Precision') +
  scale_x_continuous('Interactions (thousands)', labels = function(x) x / 1e3) +
  scale_color_manual(values = pal) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.position = 'none')
p2
ggsave("fig/euk_interactomes/human-PR-vs-HTS.pdf", p2, 
       width = 4.5, height = 4.5, units = "cm", useDingbats = FALSE)

## Venn diagram, v1 vs. v2
# get overlaps
overlap1 = v1 %>%
  left_join(v2, by = c('protein_A', 'protein_B')) %>% 
  mutate(flag = !is.na(precision.y))
overlap2 = v2 %>%
  left_join(v1, by = c('protein_A', 'protein_B')) %>% 
  mutate(flag = !is.na(precision.y))

## Venn diagram
venn_file = "fig/euk_interactomes/human-venn-v1-v2.pdf"
pdf(venn_file, width = 1.15, height = 1.15)
draw.pairwise.venn(area1 = nrow(v1),
                   area2 = nrow(v2),
                   cross.area = sum(overlap1$flag),
                   category = c('Version 1', 'Version 2'),
                   scaled = TRUE,
                   lwd = c(0.5, 0.5),
                   col = c('grey20', 'grey20'),
                   fill = version_pal,
                   fontfamily = rep('sans', 3),
                   cat.fontfamily = rep('sans', 2),
                   # alpha = rep(1, 2),
                   cex = rep(0.514, 3),
                   cat.cex = rep(0.514, 2))
dev.off()

## number of publications per interaction
# read all known interactions 
db_files = list.files("~/git/network-validation/data/database", 
                      pattern = '*-human\\.txt\\.gz$', full.names = TRUE)
dbs = map(db_files, read.delim)
db = bind_rows(dbs)
## alphabetize
sorted = t(apply(db[, 1:2], 1, sort))
db$gene_A = sorted[, 1]
db$gene_B = sorted[, 2]
db %<>% distinct(gene_A, gene_B, pmid)
# count PMIDs per interaction
known = db %>%
  count(gene_A, gene_B)
# subtract number of NAs
known_na = db %>%
  mutate(is_na = is.na(pmid)) %>%
  filter(is_na) %>%
  dplyr::select(-pmid) %>%
  distinct()
known %<>% 
  left_join(known_na) %>%
  replace_na(list(is_na = FALSE)) %>%
  mutate(n2 = ifelse(is_na, n - 1, n))

# flag number of PMIDs
pubs = v2 %>%
  dplyr::rename(gene_A = protein_A, gene_B = protein_B) %>% 
  left_join(known) %>%
  replace_na(list(n = 0, n2 = 0))

## plot as a bar chart instead
pal = brewer.pal(6, 'Blues')
p3 = pubs %>%
  mutate(fill = ifelse(n >= 5, 5, n) %>% 
           as.character() %>%
           fct_recode('0' = '0',
                      '5+' = '5')) %>%
  ggplot(aes(x = '1', fill = fill)) +
  geom_segment(aes(x = -Inf, xend = -Inf, y = 0, yend = 1),
               color = 'grey50') +
  geom_bar(color = 'grey15', size = 0.15, position = 'fill') +
  scale_fill_manual('Publications', values = pal) +
  scale_y_reverse('% of interactions', # expand = c(0, 0),
                  breaks = seq(0, 1, 0.25), 
                  labels = function(x) 100 * (1 - x)) +
  coord_flip() +
  clean_theme() +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(0.4, 'lines'),
        legend.position = 'top')
p3
ggsave("fig/euk_interactomes/human-novel-interactions-bar.pdf",
       p3, width = 4.25, height = 2.5, units = 'cm', useDingbats = FALSE)

## number of publications per interaction, v1 only vs. v2 only 
union = inner_join(v1, v2, by = c('protein_A', 'protein_B'))
v1_only = anti_join(v1, v2, by = c('protein_A', 'protein_B'))
v2_only = anti_join(v2, v1, by = c('protein_A', 'protein_B'))
pubs = bind_rows(
  union %>% mutate(xval = 'Both'),
  v1_only %>% mutate(xval = 'Version 1 only'),
  v2_only %>% mutate(xval = 'Version 2 only')
) %>% 
  mutate(xval = fct_relevel(xval, 'Both', 'Version 2 only', 'Version 1 only')) %>% 
  dplyr::rename(gene_A = protein_A, gene_B = protein_B) %>% 
  left_join(known) %>%
  replace_na(list(n = 0, n2 = 0))
pal = brewer.pal(6, 'Blues')
p4 = pubs %>%
  mutate(fill = ifelse(n >= 5, 5, n) %>% 
           as.character() %>%
           fct_recode('0' = '0',
                      '5+' = '5')) %>%
  ggplot(aes(x = xval, fill = fill)) +
  geom_segment(aes(x = -Inf, xend = -Inf, y = 0, yend = 1),
               color = 'grey50') +
  geom_bar(color = 'grey15', size = 0.15, position = 'fill', width = 0.8) +
  scale_fill_manual('Publications', values = pal) +
  scale_y_reverse('% of interactions', # expand = c(0, 0),
                  breaks = seq(0, 1, 0.25), 
                  labels = function(x) 100 * (1 - x)) +
  coord_flip() +
  clean_theme() +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        # axis.text.y = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, 'lines'),
        legend.key.size = unit(0.4, 'lines'),
        legend.position = 'top')
p4
ggsave("fig/euk_interactomes/human-novel-interactions-v1-v2-bar.pdf",
       p4, width = 4.25, height = 2.75, units = 'cm', useDingbats = FALSE)

library(broom)
pubs %>% 
  filter(xval != 'Both') %>% 
  do(tidy(chisq.test(.$xval, .$n == '0')))
