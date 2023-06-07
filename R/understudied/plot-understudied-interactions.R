# Plot the number of understudied proteins with interactions in CF-MS networks
# (2021 vs. CFdb)
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(igraph)
source("R/theme.R")

# read gene2pubmed
g2p = read.delim("data/resources/uncharacterized/gene2pubmed.gz") %>% 
  filter(X.tax_id == 9606) %>% 
  dplyr::count(GeneID)
# map to gene names
map = read_tsv("data/identifier/HUMAN_9606_idmapping.dat.gz",
               col_names = c("uniprot", "db", "id"))
gn = filter(map, db == 'Gene_Name') %>% dplyr::select(-db)
id = filter(map, db == 'GeneID') %>% dplyr::select(-db)
gn2id = left_join(gn, id, by = 'uniprot') %>%
  dplyr::select(starts_with('id')) %>%
  set_colnames(c('gene', 'id'))
# execute mapping
g2p %<>%
  dplyr::rename(id = GeneID, n_pubs = n) %>% 
  mutate(id = as.character(id)) %>% 
  left_join(gn2id, by = 'id') %>% 
  drop_na(gene) %>% 
  distinct(gene, n_pubs)
# get understudied genes
understudied_genes = g2p %>% filter(n_pubs <= 10) %>% pull(gene)

# read networks
v1 = read.delim("data/euk_interactomes/CF-MS-interactome.tsv.gz") %>% 
  mutate(version = 'Version 1')
v2 = readRDS("data/euk_interactomes/networks/network-species=Homo_sapiens-n_replicates=166-version=V2-n_features=1.rds") %>% 
  mutate(version = 'Version 2')

# bar chart: percentage of understudied proteins with an interaction in v2
cdf = readRDS("data/QC/genes-CDF.rds") %>% filter(species == 'Homo sapiens')
understudied = g2p %>%
  filter(n_pubs <= 10) %>% 
  distinct(gene) %>% 
  mutate(has_interaction = gene %in% names(degrees),
         ever_detected = gene %in% cdf$gene,
         color = ifelse(has_interaction, 'Interaction detected',
                        ifelse(ever_detected, 'Protein detected', 
                               'Not detected')) %>% 
           fct_relevel('Interaction detected',
                       'Protein detected', 
                       'Not detected'))
pal = c('Not detected' = 'grey88',
        'Protein detected' = pals::kovesi.linear_gow_65_90_c35(100)[75], # 'grey60',
        'Interaction detected' = pals::kovesi.linear_gow_65_90_c35(100)[40]
        )
counts = dplyr::count(understudied, color)
labels = with(counts, paste0(color, ' (n = ', format(n, big.mark = ',') %>% 
                               trimws(), ')')) %>% 
  setNames(counts$color)
p = understudied %>% 
  ggplot(aes(x = '1', fill = color)) +
  geom_bar(size = 0.15, color = 'grey20') + 
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 2.5e3),
               color = 'grey50', size = 0.4) +
  scale_y_continuous('Understudied proteins', expand = c(0, 200),
                     limits = c(0, 2800),
                     breaks = seq(0, 2500, 500)) +
  scale_fill_manual('', values = pal, breaks = names(pal), labels = labels) +
  guides(fill = guide_legend(nrow = 3)) +
  coord_flip() +
  boxed_theme() +
  theme(panel.border = element_blank(),
        aspect.ratio = 0.2,
        legend.key.size = unit(0.3, 'lines'),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
p
ggsave("fig/understudied/interaction-vs-protein-detected-bar-chart.pdf", p,
       width = 3.75, height = 3, units = "cm", useDingbats = FALSE)
