# Plot the breakdown of phospho-STY sites, in aggregate and by species.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read localization scores
dat = readRDS("data/QC/phosphosites-localization-scores.rds") %>% 
  filter(n_fractions > 0) %>% 
  # cut fractions at a max of 10
  mutate(n_fractions = winsorize(n_fractions, c(NA, 10))) %>% 
  # flag STY
  mutate(sty = strsplit(gene, '-') %>% map_chr(~ extract(., length(.) - 1))) %>% 
  # recode one-off species
  mutate(species = fct_recode(species,
                              'Other\nplant' = 'Solanum lycopersicum',
                              'Other\nplant' = 'Selaginella moellendorffii',
                              'Other\nplant' = 'Ceratopteris richardii',
                              'Other\nplant' = 'Glycine max',
                              'Other\nplant' = 'Brassica oleracea var. italica',
                              'Other\nplant' = 'Zea mays',
                              'S. cerevisiae' = 'Saccharomyces cerevisiae',
                              'Other\nfungus' = 'Chaetomium thermophilum',
                              'Bacteria' = 'Cyanothece ATCC 51142',
                              'Plasmodium\nspp.' = 'Plasmodium berghei',
                              'Plasmodium\nspp.' = 'Plasmodium falciparum',
                              'Plasmodium\nspp.' = 'Plasmodium knowlesi',
                              'H. sapiens' = 'Homo sapiens',
                              'H. sapiens' = 'Homo sapiens/Salmonella enterica',
                              'A. thaliana' = 'Arabidopsis thaliana',
                              'T. brucei' = 'Trypanosoma brucei',
                              'M. musculus' = 'Mus musculus',
                              'C. elegans' = 'Caenorhabditis elegans',
                              'X. laevis' = 'Xenopus laevis',
                              'T. aestivum' = 'Triticum aestivum',
                              'Other\neukaryote' = 'Strongylocentrotus purpuratus',
                              'Other\neukaryote' = 'Dictyostelium discoideum',
                              'Other\nplant' = 'Chlamydomonas reinhardtii',
                              'Other\neukaryote' = 'Drosophila melanogaster',
                              'Other\neukaryote' = 'Nematostella vectensis',
                              'O. sativa' = 'Oryza sativa',
                              'Bacteria' = 'Synechocystis sp. PCC 6803', 
                              'Bacteria' = 'Kuenenia stuttgartiensis',
                              'Bacteria' = 'Escherichia coli',
                              'Bacteria' = 'Anabaena sp. 7120',
                              'Bacteria' = 'Salmonella typhimurium SL1344',
                              'Other\nfungus' = 'Podospora anserina',
                              'Other\nplant' = 'Gossypium hirsutum',
                              'Other\neukaryote' = 'Rattus norvegicus',
                              'Other\neukaryote' = 'Plasmodium\nspp.')) %>% 
  # let's just do it manually to get it right
  mutate(species = fct_relevel(species, 'H. sapiens', 'A. thaliana',
                               'M. musculus', 'Bacteria', 
                               'C. elegans', 'T. brucei', 'S. cerevisiae',
                               'X. laevis', 'O. sativa', 'T. aestivum', 
                               'Other\nplant', 'Other\nfungus', 'Other\neukaryote',
                               'Prokaryote'))

## pie chart, aggregate
pal = pals::stepped3()[c(1, 5) + 2] %>% c(pals::kelly()[1])
p1 = dat %>% 
  distinct(species, gene, sty) %>% 
  ggplot(aes(x = '1', fill = sty)) +
  geom_bar(size = 0.15, color = 'grey20') + 
  coord_polar(theta = 'y') +
  scale_fill_manual('', labels = ~ paste0('p', .), values = pal) +
  boxed_theme() +
  theme(panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.key.size = unit(0.35, 'lines'),
        aspect.ratio = 1)
p1
ggsave("fig/phosphorylation/STY-overall-pie.pdf", p1,
       width = 4, height = 3, units = "cm", useDingbats = FALSE)

## bar chart ####
p2 = dat %>% 
  distinct(species, gene, sty) %>% 
  mutate(sty = factor(sty, levels = rev(c('S', 'T', 'Y')))) %>% 
  ggplot(aes(x = species, fill = sty)) +
  geom_bar(position = 'fill', color = 'grey20', size = 0.15, width = 0.8) +
  scale_fill_manual('', labels = ~ paste0('p', .), values = pal,
                    limits = c('S', 'T', 'Y')) +
  scale_x_discrete(labels = ~ gsub("\n", " ", .)) +
  scale_y_continuous('% of phosphosites', labels = ~ . * 100, 
                     expand = c(0, 0)) +
  boxed_theme() +
  theme(aspect.ratio = 0.8,
        legend.key.size = unit(0.35, 'lines'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())
p2
ggsave("fig/phosphorylation/STY-by-species-bar.pdf", p2,
       width = 18, height = 5, units = "cm", useDingbats = FALSE)
