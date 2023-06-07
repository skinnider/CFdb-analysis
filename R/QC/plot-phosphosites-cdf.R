# Plot a CDF of phosphosites vs. # of fractions in each experiment and 
# on average.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read experiments
expts = read.csv("~/git/CFdb-searches/data/experiments.csv") %>% 
  mutate(Replicate = fct_recode(Replicate,
                                'BN3to12_RAW' = 'BN3to12',
                                'BN4to16_RAW' = 'BN4to16',
                                'Ghosts_DDM_BiobasicSEC-2' = 'Ghosts_6061_DDM_BiobasicSEC-2',
                                'Hemolysate_old_9330_SEC' = 'Hemolysate_9330_old_SEC',
                                'whole_rbc_IP_lysis_6840' = 'whole_rbc_IP_lysis',
                                'Hemolysate_8994_IEX-3' = 'Hemolysate_8994_IEX') %>% 
           as.character())

# read protein quantitation CDF
cdf = readRDS("data/QC/phosphosites-CDF.rds") %>% 
  # drop duplicates
  filter(!(grepl("-[2-6]$", experiment) & 
             accession %in% c('PXD027704', 'PXD023211', 'MSV000082468')))

# keep intensity only
best = cdf %>% 
  filter(grepl('intensity', quant_mode, ignore.case = TRUE))

# calculate the mean curve
mean = best %>%
  group_by(n_fractions) %>%
  summarise(mean = mean(n_phosphosites, na.rm = T),
            sd = sd(n_phosphosites, na.rm = T),
            median = median(n_phosphosites, na.rm = T),
            q1 = quantile(n_phosphosites, na.rm = T, probs = 0.25),
            q3 = quantile(n_phosphosites, na.rm = T, probs = 0.75),
            cov = mean(coverage, na.rm = T),
            n = n()) %>%
  ungroup()

# calculate the mean curve for proteins
mean_proteins = best %>%
  group_by(n_fractions) %>%
  summarise(mean = mean(n_proteins, na.rm = T),
            sd = sd(n_proteins, na.rm = T),
            median = median(n_proteins, na.rm = T),
            q1 = quantile(n_proteins, na.rm = T, probs = 0.25),
            q3 = quantile(n_proteins, na.rm = T, probs = 0.75),
            cov = mean(coverage, na.rm = T),
            n = n()) %>%
  ungroup()

# plot each curve, plus the mean curve
## label n=10 fractions
lab = filter(mean, n_fractions == 1) %>%
  mutate(text = paste0(format(round(mean), big.mark = ','), 
                       ' phosphosites'))
col = BuenColors::jdb_palette("solar_extra") %>% extract(1)
p1 = best %>%
  unite(group, accession, experiment, remove = FALSE) %>%
  ggplot(aes(x = n_fractions, y = n_phosphosites)) +
  geom_path(aes(group = group), color = 'grey90', size = 0.1) + 
  geom_path(data = mean, aes(y = mean), color = col, size = 0.75) +
  geom_label_repel(data = lab, aes(label = text, y = mean), size = 1.75,
                   segment.size = 0.25, label.padding = 0.1, label.size = NA,
                   nudge_y = +4000, nudge_x = +50, hjust = 0, fill = NA,
                   lineheight = 1) +
  geom_point(data = lab, aes(y = mean), size = 0.8, color = col) +
  scale_x_continuous('Fractions', expand = c(0, 0), limits = c(1, 50),
                     breaks = c(1, seq(10, 50, 10)),
                     labels = c(1, seq(10, 40, 10), '>50')) +
  scale_y_continuous('Phosphosites', # limits = c(0, 6400), 
                     breaks = seq(0, 1600, 400),
                     expand = c(0, 0)) +
  coord_cartesian(xlim = c(1, 50), ylim = c(0, 1.4e3)) +
  boxed_theme()
p1
ggsave("fig/QC/phosphosites-cdf.pdf", p1, width = 4.6, height = 4.2, 
       units = "cm", useDingbats = F)

# plot number of phosphoproteins
## label n=10 fractions
lab = filter(mean_proteins, n_fractions == 1) %>%
  mutate(text = paste0(format(round(mean), big.mark = ','), 
                       ' phosphoproteins'))
p2 = best %>%
  unite(group, accession, experiment, remove = FALSE) %>%
  ggplot(aes(x = n_fractions, y = n_proteins)) +
  geom_path(aes(group = group), color = 'grey90', size = 0.1) + 
  geom_path(data = mean_proteins, aes(y = mean), color = col, size = 0.75) +
  geom_label_repel(data = lab, aes(label = text, y = mean), size = 1.75,
                   segment.size = 0.25, label.padding = 0.1, label.size = NA,
                   nudge_y = +4000, nudge_x = +50, hjust = 0, fill = NA,
                   lineheight = 1) +
  geom_point(data = lab, aes(y = mean), size = 0.8, color = col) +
  scale_x_continuous('Fractions', expand = c(0, 0), limits = c(1, 50),
                     breaks = c(1, seq(10, 50, 10)),
                     labels = c(1, seq(10, 40, 10), '>50')) +
  scale_y_continuous('Phosphoproteins', # limits = c(0, 6400), 
                     expand = c(0, 0)) +
  coord_cartesian(xlim = c(1, 50), ylim = c(0, 630)) +
  boxed_theme()
p2
ggsave("fig/QC/phosphoproteins-cdf.pdf", p2, width = 4.6, height = 4.2, 
       units = "cm", useDingbats = F)

# plot number of phosphosites in 1+ fractions per species
species = best %>% 
  filter(n_fractions == 1) %>% 
  dplyr::rename(Accession = accession, Replicate = experiment) %>% 
  left_join(expts, by = c('Accession', 'Replicate')) %>% 
  # recode one-off species
  mutate(species = fct_recode(Species,
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
species_pal = pals::kelly()
species_color = species_pal 
species_color[1] = 'grey30'
p3 = species %>% 
  ggplot(aes(x = species, y = n_phosphosites, color = species, fill = species)) +
  geom_boxplot(alpha = 0.4, width = 0.6, size = 0.3, outlier.shape = NA) + 
  scale_x_discrete(labels = ~ gsub("\n", " ", .)) +
  scale_y_continuous('# of phosphosites', breaks = seq(0, 1200, 300)) +
  scale_fill_manual(values = species_pal) +
  scale_color_manual(values = species_color) +
  coord_cartesian(ylim = c(0, 1400)) +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = 'none',
        aspect.ratio = 0.8)
p3
ggsave("fig/QC/phosphosites-by-species.pdf", p3, width = 7, height = 5,
       units = "cm", useDingbats = FALSE)

# repeat, but for phosphoproteins
p4 = species %>% 
  ggplot(aes(x = species, y = n_proteins, color = species, fill = species)) +
  geom_boxplot(alpha = 0.4, width = 0.6, size = 0.3, outlier.shape = NA) + 
  scale_x_discrete(labels = ~ gsub("\n", " ", .)) +
  scale_y_continuous('# of phosphoproteins') +
  scale_fill_manual(values = species_pal) +
  scale_color_manual(values = species_color) +
  coord_cartesian(ylim = c(0, 670)) +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = 'none',
        aspect.ratio = 0.8)
p4
ggsave("fig/QC/phosphoproteins-by-species.pdf", p4, width = 7, height = 5,
       units = "cm", useDingbats = FALSE)

# histogram: number of phosphosites
best0 = filter(best, n_fractions == 1)
median_phosphosites = median(best0$n_phosphosites)
mean_phosphosites = mean(best0$n_phosphosites)
p5 = best0 %>% 
  mutate(n_phosphosites = winsorize(n_phosphosites, limits = c(NA, 2000))) %>% 
  ggplot(aes(x = n_phosphosites)) +
  geom_histogram(aes(fill = '1', color = '1'), 
                 bins = 30, alpha = 0.4, size = 0.3) + 
  geom_vline(aes(xintercept = mean_phosphosites), color = 'black', size = 0.3,
             linetype = 'dotted') +
  geom_label(data = data.frame(),
             aes(x = mean_phosphosites, y = Inf,
                 label = paste(round(mean_phosphosites), 'phosphosites')),
             hjust = 0, vjust = 1, size = 1.75, fill = NA, label.size = NA,
             label.padding = unit(0.45, 'lines')) + 
  scale_color_manual(values = 'grey75') +
  scale_fill_manual(values = 'grey75') +
  scale_x_continuous('# of phosphosites', 
                     labels = ~ replace(., . == 2000, '>2000')) +
  scale_y_continuous('Experiments', expand = c(0, 0), limits = c(0, 100)) +
  boxed_theme() +
  theme(legend.position = 'none',
        aspect.ratio = 0.7)
p5

# histogram: number of phosphoproteins
median_phosphoproteins = median(best0$n_proteins)
mean_phosphoproteins = mean(best0$n_proteins)
p6 = best0 %>% 
  mutate(n_proteins = winsorize(n_proteins, limits = c(NA, 2000))) %>% 
  ggplot(aes(x = n_proteins)) +
  geom_histogram(aes(fill = '1', color = '1'), 
                 bins = 30, alpha = 0.4, size = 0.3) + 
  geom_vline(aes(xintercept = mean_phosphoproteins), color = 'black', size = 0.3,
             linetype = 'dotted') +
  geom_label(data = data.frame(),
             aes(x = mean_phosphoproteins, y = Inf,
                 label = paste(round(mean_phosphoproteins), 'phosphoproteins')),
             hjust = 0, vjust = 1, size = 1.75, fill = NA, label.size = NA,
             label.padding = unit(0.45, 'lines')) + 
  scale_color_manual(values = 'grey75') +
  scale_fill_manual(values = 'grey75') +
  scale_x_continuous('# of proteins', 
                     labels = ~ replace(., . == 2000, '>2000')) +
  scale_y_continuous('Experiments', expand = c(0, 0), limits = c(0, 100)) +
  boxed_theme() +
  theme(legend.position = 'none',
        aspect.ratio = 0.7)
p6

# combine 
p = p5 | p6
ggsave("fig/QC/phosphosites-histograms.pdf", p,
       width = 10, height = 4, units = "cm", useDingbats = FALSE)
