# Plot a CDF of the number of _experiments_ containing at least n 
# phosphosites.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read experiments
expts = read.csv("~/git/CFdb-searches/data/experiments.csv")
versions = distinct(expts, Accession, Replicate, Version) %>% 
  dplyr::rename(version = Version, 
                accession = Accession, 
                experiment = Replicate) %>% 
  mutate(experiment = fct_recode(experiment, 
                                 'BN3to12_RAW' = 'BN3to12',
                                 'BN4to16_RAW' = 'BN4to16',
                                 'Ghosts_DDM_BiobasicSEC-2' = 'Ghosts_6061_DDM_BiobasicSEC-2',
                                 'Hemolysate_old_9330_SEC' = 'Hemolysate_9330_old_SEC',
                                 'whole_rbc_IP_lysis_6840' = 'whole_rbc_IP_lysis',
                                 'Hemolysate_8994_IEX-3' = 'Hemolysate_8994_IEX') %>% 
           as.character())

# read CDF
cdf = readRDS('data/QC/phosphosites-CDF.rds') %>% 
  filter(quant_mode == 'intensity') %>% 
  # drop duplicates
  filter(!(grepl("-[2-6]$", experiment) & accession %in% c('PXD027704',
                                                           'PXD023211',
                                                           'MSV000082468')))

# reshape to # datasets/# of phosphosites
dataset_cdf = tidyr::crossing(n_phosphosites = seq(0, 2.5e3, 10),
                              n_fractions = c(seq_len(10)))
dataset_cdf$n_datasets = map2_int(
  dataset_cdf$n_phosphosites, dataset_cdf$n_fractions,
  ~ cdf %>%
    filter(n_fractions == .y) %>%
    filter(n_phosphosites > .x) %>%
    distinct(accession, experiment) %>%
    nrow()
)

# plot
pal = BuenColors::jdb_palette("wolfgang_extra") %>% colorRampPalette() 
p1 = dataset_cdf %>%
  mutate(n_fractions = as.factor(n_fractions)) %>% 
  arrange(desc(n_fractions)) %>% 
  filter(n_datasets > 0) %>% 
  ggplot(aes(x = n_phosphosites, y = n_datasets, color = factor(n_fractions))) +
  geom_path(size = 0.35) + 
  scale_color_manual("min. # of fractions", values = pal(12) %>% rev) +
  scale_x_continuous("# of phosphosites") +
  scale_y_continuous("# of datasets") +
  guides(color = guide_legend(title.position = 'top', title.hjust = 0.5,
                              nrow = 1)) +
  boxed_theme()
p1
ggsave("fig/QC/phosphosites-cdf-datasets.pdf", p1,
       width = 6, height = 5.5, units = "cm", useDingbats = F)

# repeat for phosphoproteins
dataset_cdf = tidyr::crossing(n_proteins = seq(0, 2.5e3, 10),
                              n_fractions = c(seq_len(10)))
dataset_cdf$n_datasets = map2_int(
  dataset_cdf$n_proteins, dataset_cdf$n_fractions,
  ~ cdf %>%
    filter(n_fractions == .y) %>%
    filter(n_proteins > .x) %>%
    distinct(accession, experiment) %>%
    nrow()
)
p2 = dataset_cdf %>%
  mutate(n_fractions = as.factor(n_fractions)) %>% 
  arrange(desc(n_fractions)) %>% 
  filter(n_datasets > 0) %>% 
  ggplot(aes(x = n_proteins, y = n_datasets, color = factor(n_fractions))) +
  geom_path(size = 0.35) + 
  scale_color_manual("min. # of fractions", values = pal(12) %>% rev) +
  scale_x_continuous("# of phosphoproteins") +
  scale_y_continuous("# of datasets") +
  coord_cartesian(xlim = c(0, 1000)) +
  guides(color = guide_legend(title.position = 'top', title.hjust = 0.5,
                              nrow = 1)) +
  boxed_theme()
p2
ggsave("fig/QC/phosphoproteins-cdf-datasets.pdf", p2,
       width = 6, height = 5.5, units = "cm", useDingbats = F)
