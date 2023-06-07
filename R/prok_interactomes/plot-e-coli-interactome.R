setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(flavin)
library(PrInCE)
source("R/theme.R")

## PR curve, v2 vs. Y2H/AP-MS
# plot v2 against E. coli high-throughput interactomes
interactomes = c('Butland2005',
                 'Arifuzzaman2006',
                 'Hu2009',
                 'Rajagopala2014',
                 'Babu2017')
hts_files = paste0("data/resources/interactomes/", interactomes, ".csv")
hts_ppis = map(hts_files, read.csv) %>%
  setNames(gsub("\\..*$", "", basename(hts_files)))

# calculate precision and recall
complexes = read.csv("data/prok_interactomes/EcoCyc/complexes_UP000000625.csv") %>% 
  as_annotation_list('gene_name', 'complex')
adj = adjacency_matrix_from_list(complexes)
precision = map_dbl(hts_ppis, ~ {
  make_labels(adj, .) %>% 
    mean(na.rm = TRUE)
})
recall = map_int(hts_ppis, nrow)
n_proteins = map_int(hts_ppis, ~ with(., n_distinct(c(protein_A, protein_B))))
hts = data.frame(interactome = names(hts_ppis),
                 precision = precision,
                 recall = recall,
                 n_proteins = n_proteins)

# read PR curve
pr = readRDS("data/prok_interactomes/e-coli-PR.rds") %>% 
  mutate(idx = row_number())

pal = pals::stepped3()[seq(1, 19, 4)] %>% c('black') %>% setNames(hts$interactome)
p1 = pr %>%
  filter(idx %% 10 == 0) %>% 
  ggplot(aes(x = idx, y = precision)) +
  geom_line(linetype = 'dotted', linewidth = 0.3) +
  # geom_path(linetype = 'dotted', size = 0.3) +
  geom_point(data = hts, aes(x = recall, color = interactome), size = 0.9) +
  geom_label_repel(data = hts, 
                   aes(x = recall, color = interactome, label = interactome),
                   size = 2.25, label.padding = 0.25, label.size = NA,
                   fill = NA, min.segment.length = 0, segment.size = 0.4) +
  scale_y_continuous('Precision') +
  scale_x_continuous('Interactions (thousands)', labels = function(x) x / 1e3) +
  scale_color_manual(values = pal) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.position = 'none')
p1
ggsave("fig/prok_interactomes/e-coli-PR-vs-HTS.pdf", p1, 
       width = 4.5, height = 4.5, units = "cm", useDingbats = FALSE)
