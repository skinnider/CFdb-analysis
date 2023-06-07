# Plot panels showing the increase in the dataset from 2021 to CFdb.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(data.table)
source("R/theme.R")

###############################################################################-
## plot old vs new summary statistics ####
###############################################################################-

# list all summary files
files = list.files('~/git/CFdb-searches/data/maxquant', pattern = 'summary.txt', 
                   full.names = T, recursive = T) %>%
  # ignore phosphorylation
  extract(!grepl("phosphorylation", .))

# loop over files
summary = map_dfr(files, ~ {
  tab = fread(.) %>% filter(Experiment != "")
  
  ## number of input raw files
  n_files = n_distinct(tab$`Raw file`)
  
  ## number of fractions
  n_fractions = n_distinct(tab$Experiment)
  
  ## number of MS/MS spectra
  msms = sum(tab$`MS/MS`)
  
  ## number of sequenced peptides
  peptides = sum(tab$`MS/MS Identified`)
  
  data.frame(file = ., raw_files = n_files, fractions = n_fractions,
             msms = msms, peptides = peptides)
})

# flag accession, replicate, and search
summary %<>%
  mutate(file = gsub("^.*maxquant/", "", file)) %>%
  separate(file, c("accession", "replicate", "search", "x"), '/') %>%
  dplyr::select(-x)

# merge in versions 1 vs 2
expts = read.csv("~/git/CFdb-searches/data/experiments.csv") %>% 
  # drop channels
  mutate(Replicate = gsub("-heavy|-medium|-[1-6]$", "", Replicate) %>% 
           fct_recode('BN3to12_RAW' = 'BN3to12',
                      'BN4to16_RAW' = 'BN4to16',
                      'Ghosts_DDM_BiobasicSEC-2' = 'Ghosts_6061_DDM_BiobasicSEC',
                      'Hemolysate_old_9330_SEC' = 'Hemolysate_9330_old_SEC',
                      'whole_rbc_IP_lysis_6840' = 'whole_rbc_IP_lysis',
                      'Hemolysate_8994_IEX-3' = 'Hemolysate_8994_IEX') %>% 
           as.character()) %>% 
  # fix one species
  mutate(Species = fct_recode(Species, 'Salmonella enterica' = 
                                'Salmonella typhimurium SL1344') %>% 
           as.character())
colnames(expts) %<>% tolower()
expts %<>% left_join(summary, by = c('accession', 'replicate'))
# confirm no missing data
filter(expts, is.na(msms))

# o plot number of species ####
version_pal = c('grey88', pals::stepped()[3])
species = expts %>% 
  distinct(species, version) %>% 
  group_by(species) %>% 
  arrange(version) %>% 
  dplyr::slice(1) %>% 
  ungroup()
p0 = species %>% 
  mutate(version = gsub("V", "Version ", version) %>% 
           fct_relevel('Version 2', 'Version 1')) %>% 
  ggplot(aes(x = '1', fill = version)) +
  geom_bar(color = 'grey20', size = 0.15) +
  geom_segment(aes(y = 0, yend = 30, x = -Inf, xend = -Inf), color = 'grey50',
               size = 0.4) +
  scale_y_continuous('# of species') +
  scale_fill_manual('', values = version_pal,
                    breaks = c('Version 1', 'Version 2')) +
  boxed_theme() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15, 'lines'),
        aspect.ratio = 4,
        legend.position = 'right',
        legend.key.width = unit(0.35, 'lines'),
        legend.key.height = unit(0.35, 'lines'),
        panel.border = element_blank())
p0

# o plot number of experiments searched ####
p1 = expts %>% 
  mutate(version = gsub("V", "Version ", version) %>% 
           fct_relevel('Version 2', 'Version 1')) %>% 
  ggplot(aes(x = '1', fill = version)) +
  geom_bar(color = 'grey20', size = 0.15) +
  geom_segment(aes(y = 0, yend = 400, x = -Inf, xend = -Inf), color = 'grey50',
               size = 0.4) +
  scale_y_continuous('# of experiments') +
  scale_fill_manual('', values = version_pal,
                    breaks = c('Version 1', 'Version 2')) +
  boxed_theme() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15, 'lines'),
        aspect.ratio = 4,
        legend.position = 'right',
        legend.key.width = unit(0.35, 'lines'),
        legend.key.height = unit(0.35, 'lines'),
        panel.border = element_blank())
p1

# o plot number of fractions searched ####
n_fractions = expts %>% 
  group_by(version) %>% 
  summarise(n = sum(fractions, na.rm = TRUE)) %>% 
  ungroup()
p2 = n_fractions %>%
  mutate(version = gsub("V", "Version ", version) %>% 
           fct_relevel('Version 2', 'Version 1')) %>% 
  ggplot(aes(x = '1', y = n, fill = version)) +
  geom_col(color = 'grey20', size = 0.15) +
  geom_segment(aes(y = 0, yend = 2.5e4, x = -Inf, xend = -Inf), color = 'grey50',
               size = 0.4) +
  scale_y_continuous(expression('# of fractions'~(10^3)), 
                     labels = function(x) x / 1e3) +
  scale_fill_manual('', values = version_pal,
                    breaks = c('Version 1', 'Version 2')) +
  boxed_theme() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15, 'lines'),
        aspect.ratio = 4,
        legend.position = 'right',
        legend.key.width = unit(0.35, 'lines'),
        legend.key.height = unit(0.35, 'lines'),
        panel.border = element_blank())
p2

# o plot number of papers included ####
n_papers = expts %>% 
  distinct(author, year, version) %>% 
  mutate(year = as.character(year)) %>% 
  bind_rows(data.frame(author = 'Lee', year = '2021b', version = 'V2')) %>% 
  dplyr::count(version)
p3 = n_papers %>%
  mutate(version = gsub("V", "Version ", version) %>% 
           fct_relevel('Version 2', 'Version 1')) %>% 
  ggplot(aes(x = '1', y = n, fill = version)) +
  geom_col(color = 'grey20', size = 0.15) +
  geom_segment(aes(y = 0, yend = 60, x = -Inf, xend = -Inf), color = 'grey50',
               size = 0.4) +
  scale_y_continuous('# of studies') +
  scale_fill_manual('', values = version_pal,
                    breaks = c('Version 1', 'Version 2')) +
  boxed_theme() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15, 'lines'),
        aspect.ratio = 4,
        legend.position = 'right',
        legend.key.width = unit(0.35, 'lines'),
        legend.key.height = unit(0.35, 'lines'),
        panel.border = element_blank())
p3

# o plot number of sequenced peptides ####
n_peptides = expts %>% 
  group_by(version) %>% 
  summarise(n = sum(peptides, na.rm = TRUE)) %>% 
  ungroup()
p4 = n_peptides %>%
  mutate(version = gsub("V", "Version ", version) %>% 
           fct_relevel('Version 2', 'Version 1')) %>% 
  ggplot(aes(x = '1', y = n, fill = version)) +
  geom_col(color = 'grey20', size = 0.15) +
  geom_segment(aes(y = 0, yend = 150e6, x = -Inf, xend = -Inf), color = 'grey50',
               size = 0.4) +
  scale_y_continuous(expression('# of peptides'~(10^6)),
                     labels = function(x) x / 1e6,
                     breaks = seq(0, 150, 50) * 1e6) +
  scale_fill_manual('', values = version_pal,
                    breaks = c('Version 1', 'Version 2')) +
  boxed_theme() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15, 'lines'),
        aspect.ratio = 4,
        legend.position = 'right',
        legend.key.width = unit(0.35, 'lines'),
        legend.key.height = unit(0.35, 'lines'),
        panel.border = element_blank())
p4

# o plot number of MS/MS spectra ####
n_msms = expts %>% 
  group_by(version) %>% 
  summarise(n = sum(msms, na.rm = TRUE)) %>% 
  ungroup()
p5 = n_msms %>%
  mutate(version = gsub("V", "Version ", version) %>% 
           fct_relevel('Version 2', 'Version 1')) %>% 
  ggplot(aes(x = '1', y = n, fill = version)) +
  geom_col(color = 'grey20', size = 0.15) +
  geom_segment(aes(y = 0, yend = 6e8, x = -Inf, xend = -Inf), color = 'grey50',
               size = 0.4) +
  scale_y_continuous(expression('# of MS/MS'~(10^6)),
                     labels = function(x) x / 1e6,
                     limits = c(0, 700e6),
                     breaks = seq(0, 600, 200) * 1e6) +
  scale_fill_manual('', values = version_pal,
                    breaks = c('Version 1', 'Version 2')) +
  boxed_theme() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15, 'lines'),
        aspect.ratio = 4,
        legend.position = 'right',
        legend.key.width = unit(0.35, 'lines'),
        legend.key.height = unit(0.35, 'lines'),
        panel.border = element_blank())
p5

# o plot instrument time ####
instrument_time = read.csv("data/QC/instrument-time.csv")
colnames(instrument_time) %<>% tolower()
time_totals = instrument_time %>% 
  group_by(author, year, accession, replicate) %>% 
  summarise(time = sum(time)) %>% 
  ungroup()
times = expts %>% 
  left_join(time_totals) %>% 
  group_by(version) %>% 
  summarise(time = sum(time, na.rm = TRUE) / 60) %>% 
  ungroup()
p6 = times %>%
  mutate(version = gsub("V", "Version ", version) %>% 
           fct_relevel('Version 2', 'Version 1')) %>% 
  ggplot(aes(x = '1', y = time, fill = version)) +
  geom_col(color = 'grey20', size = 0.15) +
  geom_segment(aes(y = 0, yend = 30e3, x = -Inf, xend = -Inf), color = 'grey50',
               size = 0.4) +
  scale_y_continuous(expression('Hours'~(10^3)),
                     labels = function(x) x / 1e3,
                     limits = c(0, 30e3),
                     breaks = seq(0, 30, 10) * 1e3) +
  scale_fill_manual('', values = version_pal,
                    breaks = c('Version 1', 'Version 2')) +
  boxed_theme() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15, 'lines'),
        aspect.ratio = 4,
        legend.position = 'right',
        legend.key.width = unit(0.35, 'lines'),
        legend.key.height = unit(0.35, 'lines'),
        panel.border = element_blank())
p6

# o plot number of observations ####
cdf = readRDS("data/QC/protein-groups-CDF.rds") 
expts = read.csv("~/git/CFdb-searches/data/experiments.csv")
colnames(expts) %<>% tolower()
expts %<>%
  dplyr::rename(experiment = replicate) %>% 
  mutate(experiment = fct_recode(experiment, 'BN3to12_RAW' = 'BN3to12',
                                 'BN4to16_RAW' = 'BN4to16',
                                 'Ghosts_DDM_BiobasicSEC-2' = 'Ghosts_6061_DDM_BiobasicSEC-2',
                                 'Hemolysate_old_9330_SEC' = 'Hemolysate_9330_old_SEC',
                                 'whole_rbc_IP_lysis_6840' = 'whole_rbc_IP_lysis',
                                 'Hemolysate_8994_IEX-3' = 'Hemolysate_8994_IEX') %>% 
           as.character())
cdf %<>% left_join(expts, by = c('accession', 'experiment'))
filter(cdf, is.na(author))
# calculate total number of observations 
n_obs = cdf %>%
  filter(grepl('intensity', quant_mode, ignore.case = TRUE)) %>% 
  filter(search == 'default') %>% 
  group_by(version) %>%
  summarise(n_quants = sum(n_proteins)) %>%
  ungroup() %>%
  arrange(desc(n_quants))
p7 = n_obs %>%
  mutate(version = gsub("V", "Version ", version) %>% 
           fct_relevel('Version 2', 'Version 1')) %>% 
  ggplot(aes(x = '1', y = n_quants, fill = version)) +
  geom_col(color = 'grey20', size = 0.15) +
  geom_segment(aes(y = 0, yend = 20e6, x = -Inf, xend = -Inf), color = 'grey50',
               size = 0.4) +
  scale_y_continuous(expression('Quantifications'~(10^6)),
                     labels = function(x) x / 1e6,
                     limits = c(0, 21e6),
                     breaks = seq(0, 20, 5) * 1e6) +
  scale_fill_manual('', values = version_pal,
                    breaks = c('Version 1', 'Version 2')) +
  boxed_theme() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15, 'lines'),
        aspect.ratio = 4,
        legend.position = 'right',
        legend.key.width = unit(0.35, 'lines'),
        legend.key.height = unit(0.35, 'lines'),
        panel.border = element_blank())
p7

# o plot number of phosphosite quantifications ####
cdf = readRDS("data/QC/phosphosites-CDF.rds") 
expts = read.csv("~/git/CFdb-searches/data/experiments.csv")
colnames(expts) %<>% tolower()
expts %<>%
  dplyr::rename(experiment = replicate) %>% 
  mutate(experiment = fct_recode(experiment, 'BN3to12_RAW' = 'BN3to12',
                                 'BN4to16_RAW' = 'BN4to16',
                                 'Ghosts_DDM_BiobasicSEC-2' = 'Ghosts_6061_DDM_BiobasicSEC-2',
                                 'Hemolysate_old_9330_SEC' = 'Hemolysate_9330_old_SEC',
                                 'whole_rbc_IP_lysis_6840' = 'whole_rbc_IP_lysis',
                                 'Hemolysate_8994_IEX-3' = 'Hemolysate_8994_IEX') %>% 
           as.character())
cdf %<>% left_join(expts, by = c('accession', 'experiment'))
filter(cdf, is.na(author))
# calculate total number of observations 
n_obs = cdf %>%
  filter(grepl('intensity', quant_mode, ignore.case = TRUE)) %>% 
  # drop duplicates
  filter(!(grepl("-[2-6]$", experiment) & accession %in% c('PXD027704',
                                                           'PXD023211',
                                                           'MSV000082468'))) %>% 
  summarise(n_quants = sum(n_phosphosites)) %>%
  arrange(desc(n_quants))
p8 = n_obs %>%
  mutate(version = 'Version 2') %>% 
  ggplot(aes(x = '1', y = n_quants, fill = version)) +
  geom_col(color = 'grey20', size = 0.15) +
  geom_segment(aes(y = 0, yend = 7.5e5, x = -Inf, xend = -Inf), color = 'grey50',
               size = 0.4) +
  scale_y_continuous(expression('Phosphosites'~(10^3)),
                     labels = function(x) x / 1e3,
                     breaks = seq(0, 750, 250) * 1e3,
                     # limits = c(0, 21e6),
                     # breaks = seq(0, 20, 5) * 1e6
                     ) +
  scale_fill_manual('', values = version_pal,
                    limits = c('Version 1', 'Version 2'),
                    breaks = c('Version 1', 'Version 2')) +
  boxed_theme() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15, 'lines'),
        aspect.ratio = 4,
        legend.position = 'right',
        legend.key.width = unit(0.35, 'lines'),
        legend.key.height = unit(0.35, 'lines'),
        panel.border = element_blank())
p8

###############################################################################-
## combine all plots ####
###############################################################################-

p_comb =   (p3 + theme(legend.position = 'none')) + 
  (p0 + theme(legend.position = 'none')) + 
  (p1 + theme(legend.position = 'none')) + 
  (p2 + theme(legend.position = 'none')) + 
  (p4 + theme(legend.position = 'none')) + 
  # (p5 + theme(legend.position = 'none')) + 
  # (p6 + theme(legend.position = 'none')) + 
  (p7 + theme(legend.position = 'none')) + 
  p8 + plot_layout(nrow = 1)
p_comb
ggsave("fig/QC/overview.pdf", p_comb, width = 13.5, height = 5,
       units = "cm", useDingbats = FALSE)
