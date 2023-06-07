setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(flavin)
library(PrInCE)
library(VennDiagram)
source("R/theme.R")

# read CF-MS network
v2 = readRDS("data/euk_interactomes/networks/network-species=Saccharomyces_cerevisiae-n_replicates=8-version=V2-n_features=2.rds") %>% 
  mutate(version = 'Version 2')

# plot against yeast high-throughput interactomes
hts_networks = c('Schwikowski2000', 
                 'Gavin2006',
                 'Krogan2006-intersection',
                 'Krogan2006-merged',
                 'Yu2008',
                 'Tarassov2008',
                 'Babu2012')
hts_files = paste0("data/resources/interactomes/", hts_networks, ".csv") %>% 
  setNames(hts_networks)
hts_ppis = map(hts_files, ~ read.csv(.x))

# calculate precision and recall
complexes = read.csv("data/euk_interactomes/CORUM/complexes_UP000002311.txt") %>%
  flavin::as_annotation_list('target', 'complex')
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
                 n_proteins = n_proteins) %>% 
  # rename Krogan
  filter(!grepl('intersection', interactome)) %>% 
  mutate(interactome = gsub("-.*$", "", interactome))

pal = pals::stepped3()[seq(1, 19, 4)] %>% c('black') %>% setNames(hts$interactome)
p1 = v2 %>%
  mutate(idx = row_number()) %>% 
  filter(idx %% 3 == 0) %>% 
  ggplot(aes(x = idx, y = precision)) +
  geom_path(linetype = '1f', linewidth = 0.3) +
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
ggsave("fig/euk_interactomes/yeast-PR-vs-HTS.pdf", p1, 
       width = 4.5, height = 4.5, units = "cm", useDingbats = FALSE)
