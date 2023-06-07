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
v1 = readRDS("data/euk_interactomes/networks/network-species=Mus_musculus-n_replicates=21-version=V1-n_features=1.rds") %>% 
  mutate(version = 'Version 1')
v2 = readRDS("data/euk_interactomes/networks/network-species=Mus_musculus-n_replicates=40-version=V2-n_features=1.rds") %>% 
  mutate(version = 'Version 2')
# combine
pr = bind_rows(v1, v2) %>% 
  group_by(version) %>% 
  mutate(idx = row_number()) %>% 
  ungroup()

# label
pct = (nrow(v2) - nrow(v1)) / nrow(v1)
lab = paste0('+', round(100 * pct, digits = 0),
             '% increase\nin mouse interactions')

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
ggsave("fig/euk_interactomes/mouse-PR-v1-v2.pdf", p1, 
       width = 4.5, height = 4.5, units = "cm", useDingbats = FALSE)

## PR curve, v2 vs. tissue interactomes
# plot RF against tissue interactomes
tissue_files = list.files("data/resources/interactomes", full.names = TRUE, 
                          pattern = 'Skinnider2021') %>% 
  extract(grepl('csv', .))
tissue_ppis = map(tissue_files, read.csv) %>%
  setNames(gsub("^.*-|\\.csv$", "", basename(tissue_files)))
# calculate precision and recall
species_map = read.csv("data/GO/species-map.csv")
proteome = filter(species_map, species == 'Mus musculus') %>% pull(proteome)
complexes = paste0("data/euk_interactomes/CORUM/complexes_", 
                   proteome, ".txt") %>% 
  read.csv() %>%
  as_annotation_list('target', 'complex')
adj = adjacency_matrix_from_list(complexes)
precision = map_dbl(tissue_ppis, ~ {
  make_labels(adj, .) %>% 
    mean(na.rm = TRUE)
})
recall = map_int(tissue_ppis, nrow)
n_proteins = map_int(tissue_ppis, ~ with(., n_distinct(c(protein_A, protein_A))))
tissue = data.frame(interactome = names(tissue_ppis),
                    precision = precision,
                    recall = recall,
                    n_proteins = n_proteins)

# pal = pals::stepped3()[seq(1, 19, 4)] %>% c('black')
p2 = v2 %>%
  mutate(idx = row_number()) %>% 
  filter(idx %% 14 == 0) %>% 
  ggplot(aes(x = idx, y = precision)) +
  geom_path(linetype = '1f', linewidth = 0.3) +
  # geom_path(linetype = 'dotted', size = 0.3) +
  geom_point(data = tissue, aes(x = recall), size = 0.9, 
             color = pals::stepped2()[13]) +
  geom_label_repel(data = tissue, aes(x = recall, label = interactome),
                   size = 2.25, label.padding = 0.25, label.size = NA,
                   fill = NA, min.segment.length = 0, segment.size = 0.4, 
                   color = pals::stepped2()[13], max.overlaps = 999) +
  scale_y_continuous('Precision') +
  scale_x_continuous('Interactions (thousands)', labels = function(x) x / 1e3) +
  # scale_color_manual(values = pal) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.position = 'none')
p2
ggsave("fig/euk_interactomes/mouse-PR-vs-tissue-interactomes.pdf", p2, 
       width = 4.5, height = 4.5, units = "cm", useDingbats = FALSE)
