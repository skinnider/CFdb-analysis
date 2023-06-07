# Plot a summary figure of all experiments analyzed here.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(data.table)
source("R/theme.R")

# list all summary files
files = list.files('~/git/CFdb-searches/data/maxquant', pattern = 'summary.txt', 
                   full.names = T, recursive = T)

# count number of input files and fractions per experiment
summary = map_dfr(files, ~ {
  tab = fread(.) %>% filter(Experiment != "")
  n_files = n_distinct(tab$`Raw file`)
  n_fractions = n_distinct(tab$Experiment)
  data.frame(file = ., raw_files = n_files, fractions = n_fractions)
})

# flag accession, replicate, and search
summary %<>%
  mutate(file = gsub("^.*maxquant/", "", file)) %>%
  separate(file, c("accession", "replicate", "search", "x"), '/') %>%
  dplyr::select(-x)

# filter out phosphorylation
summary %<>% filter(search != 'phosphorylation') %>%
  dplyr::select(-search)

# merge in species
experiments = read.csv('~/git/CFdb-searches/data/experiments.csv') %>%
  # drop channels
  mutate(Replicate = gsub("-heavy|-medium|-[1-6]$", "", Replicate) %>% 
           fct_recode('BN3to12_RAW' = 'BN3to12',
                      'BN4to16_RAW' = 'BN4to16',
                      'Ghosts_DDM_BiobasicSEC-2' = 'Ghosts_6061_DDM_BiobasicSEC',
                      'Hemolysate_old_9330_SEC' = 'Hemolysate_9330_old_SEC',
                      'whole_rbc_IP_lysis_6840' = 'whole_rbc_IP_lysis',
                      'Hemolysate_8994_IEX-3' = 'Hemolysate_8994_IEX') %>% 
           as.character())
colnames(experiments) %<>% tolower()
summary %<>% left_join(experiments, ., by = c('accession', 'replicate'))

# create data frame for plotting
df = summary %>%
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
  # group by accession
  group_by(species, accession) %>%
  arrange(fractions) %>%
  mutate(yval = factor(as.character(row_number()),
                       levels = as.character(seq_len(99))))

# facet by species and plot by accession
df %<>%
  # arrange facets by species
  group_by(species) %>%
  mutate(n_expts = n(),
         n_xvals = n_distinct(accession)) %>%
  ungroup() %>%
  arrange(desc(n_xvals), desc(n_expts)) %>%
  mutate(species = factor(species, levels = unique(species))) %>%
  # let's just do it manually to get it right
  mutate(species = fct_relevel(species, 'H. sapiens', 'A. thaliana',
                               'M. musculus', 'Bacteria', 
                               'C. elegans', 'T. brucei', 'S. cerevisiae',
                               'X. laevis', 'O. sativa', 'T. aestivum', 
                               'Other\nplant', 'Other\nfungus', 'Other\neukaryote',
                               'Prokaryote')) %>%
  # arrange accessions within facets by experiments
  group_by(species, accession) %>%
  mutate(n_yvals = n()) %>%
  ungroup() %>%
  arrange(desc(n_yvals)) %>%
  mutate(accession = factor(accession, levels = unique(.$accession))) %>%
  # bin fractions
  mutate(bin = cut(fractions, breaks = c(0, 20, 30, 40, 50, 75, 100, 999)),
         bin = fct_recode(bin, 
                          '<20' = '(0,20]',
                          '<25' = '(0,25]',
                          '21-30' = '(20,30]',
                          '31-40' = '(30,40]',
                          '41-50' = '(40,50]',
                          '26-50' = '(25,50]',
                          '51-75' = '(50,75]',
                          '76-100' = '(75,100]',
                          '>100' = '(100,999]'))
labels = df %>%
  group_by(species) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(label = paste0('n = ', n))

pal = viridis::cividis(n = 20) %>% extract(-20)
p1 = df %>%
  mutate(fractions = winsorize(fractions, c(NA, 100))) %>%
  ggplot(aes(x = accession, y = yval)) + 
  facet_grid(~ species, scales = 'free_x', space = 'free') +
  geom_tile(color = 'white', aes(fill = fractions)) + 
  geom_label(data = labels, aes(x = +Inf, y = +Inf, label = label),
             hjust = 1, vjust = 1, label.size = NA, size = 1.75,
             label.padding = unit(0.3, 'lines')) +
  scale_fill_gradientn('Fractions   ', colors = pal, na.value = 'grey90',
                       breaks = c(20, 100), labels = c('<20', '>100')) +
  scale_x_reordered('Accession') +
  scale_y_discrete('Experiments', breaks = as.character(seq(0, 30, 2))) +
  guides(fill = guide_colorbar(frame.colour = 'black', ticks = F)) +
  boxed_theme(size_sm = 5, size_lg = 6) +
  theme(axis.text.x = element_blank(), # element_text(angle = 45, hjust = 1),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.6, 'lines'),
        legend.key.height = unit(0.2, 'lines'),
        legend.key.width = unit(0.22, 'lines')
  )
p1
ggsave("fig/QC/dataset-overview.pdf", width = 18.3, height = 5.0, 
       units = "cm", useDingbats = F)

## color instead by version
version_pal = c('grey88', pals::stepped()[3]) %>% 
  setNames(c('V1', 'V2'))
p2 = df %>%
  ggplot(aes(x = accession, y = yval)) + 
  facet_grid(~ species, scales = 'free_x', space = 'free',
             # labeller = as_labeller(function(x) gsub("\n", " ", x))
  ) +
  geom_tile(color = 'white', aes(fill = version)) + 
  geom_label(data = labels, aes(x = +Inf, y = +Inf, label = label),
             hjust = 1, vjust = 1, label.size = NA, size = 1.75,
             label.padding = unit(0.3, 'lines')) +
  scale_x_reordered('Accession') +
  scale_y_discrete('Experiments', breaks = as.character(seq(0, 30, 2))) +
  scale_fill_manual('', values = version_pal, 
                    labels = ~ gsub("V", "Version ", .)) +
  guides(fill = guide_legend(override.aes = list(color = 'black', size = 0.2))) +
  boxed_theme(size_sm = 5, size_lg = 6) +
  theme(axis.text.x = element_blank(), # element_text(angle = 45, hjust = 1),
        # axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.ticks.length.x = unit(0, 'lines'),
        panel.spacing = unit(0.6, 'lines'),
        # legend.key.size = unit(0.5, "lines"),
        legend.key.height = unit(0.35, 'lines'),
        legend.key.width = unit(0.35, 'lines')
  )
p2
ggsave("fig/QC/dataset-overview-by-version.pdf", p2, 
       width = 18.3, height = 4.9, units = "cm", useDingbats = F)

# plot bar chart: total number of experiments per species
species_pal = pals::kelly()
species_color = species_pal
species_color[1] = 'grey30'
counts = df %>% count(species)
p3 = df %>% 
  ggplot(aes(x = species, color = species, fill = species)) +
  geom_bar(size = 0.2, alpha = 0.4, width = 0.75) +
  geom_segment(aes(y = 0, yend = 150, x = -Inf, xend = -Inf), color = 'grey50',
               size = 0.4) +
  geom_text(data = counts, aes(y = n, label = n), size = 1.75, vjust = 0, 
            nudge_y = +5) +
  scale_x_discrete(labels = ~ chartr('\n', ' ', .)) +
  scale_y_continuous('# of experiments', breaks = seq(0, 150, 50),
                     expand = c(0.01, 0), limits = c(0, 185)) +
  scale_fill_manual(values = species_pal) +
  scale_color_manual(values = species_color) +
  boxed_theme() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length.y = unit(0.15, 'lines'),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        aspect.ratio = 0.8,
        panel.border = element_blank())
p3
ggsave("fig/QC/dataset-overview-species.pdf", p3, width = 4.5, height = 6, 
       units = "cm")

# re-plot, colored by version
p4 = df %>% 
  mutate(version = factor(version, levels = c('V2', 'V1'))) %>% 
  ggplot(aes(x = species, fill = version)) +
  geom_bar(size = 0.2, alpha = 1, width = 0.75, color = 'grey20') +
  geom_segment(aes(y = 0, yend = 150, x = -Inf, xend = -Inf), color = 'grey50',
               size = 0.4) +
  scale_x_discrete(labels = ~ chartr('\n', ' ', .)) +
  scale_y_continuous('# of experiments', breaks = seq(0, 150, 50),
                     expand = c(0.01, 0), limits = c(0, 185)) +
  scale_fill_manual('', values = version_pal, 
                    labels = ~ gsub("V", "Version ", .)) +
  boxed_theme() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length.y = unit(0.15, 'lines'),
        axis.ticks.x = element_blank(),
        legend.position = 'top',
        aspect.ratio = 0.8,
        legend.key.width = unit(0.35, 'lines'),
        legend.key.height = unit(0.35, 'lines'),
        panel.border = element_blank())
p4
ggsave("fig/QC/dataset-overview-species-by-version.pdf", p4, 
       width = 4.5, height = 6, units = "cm")
