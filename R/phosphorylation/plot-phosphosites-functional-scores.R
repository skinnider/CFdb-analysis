# Plot the distribution of 'functional scores' as computed by Ochoa et al.
# for the phosphosites detected by CF-MS.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(readxl)
source("R/theme.R")

# read Ochoa dataset
func = read_excel("data/resources/Ochoa/41587_2019_344_MOESM5_ESM.xlsx") %>% 
  mutate(position = as.character(position))

# read localization scores
dat = readRDS("data/QC/phosphosites-localization-scores.rds") %>% 
  filter(n_fractions > 0) %>% 
  # human only
  filter(species == 'Homo sapiens') %>% 
  # separate to merge
  mutate(uniprot = strsplit(gene, '-') %>% 
           map_chr(~ head(., -2) %>% paste0(collapse = '-')),
         position = strsplit(gene, '-') %>% 
           map_chr(~ tail(., 1)))
# merge
func %<>% left_join(dat, by = c('uniprot', 'position'))

# plot: ever quantified vs. never quantified
ever_quant = func %>% 
  distinct(uniprot, position, gene, functional_score) %>% 
  mutate(ever_quant = !is.na(gene))
pal = pals::kelly()[c(1, 6)]
col = c('grey30', pal[2])
p1 = ever_quant %>% 
  ggplot(aes(x = ever_quant, y = functional_score,
             color = ever_quant, fill = ever_quant)) + 
  geom_violin(color = NA, width = 0.6, alpha = 0.7, size = 0.1) +
  geom_boxplot(alpha = 0.5, width = 0.6, size = 0.3, outlier.shape = NA,
               color = 'grey30') + 
  scale_x_discrete(labels = c('Never\ndetected', 'Ever\ndetected')) +
  scale_y_continuous('Functional score', breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = col) +
  boxed_theme() +
  theme(aspect.ratio = 1.8,
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())
p1
ggsave("fig/phosphorylation/functional-scores-ever-detected.pdf", p1,
       width = 5, height = 4, units = "cm", useDingbats = FALSE)

## number of fractions vs. functional score
func0 = func %>% 
  replace_na(list(n_fractions = 0)) %>% 
  mutate(n_fractions = winsorize(n_fractions, c(NA, 5)))
p2 = func0 %>% 
  ggplot(aes(x = factor(n_fractions), y = functional_score)) +
  geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.6, size = 0.3,
               fill = pals::kelly()[1], color = 'grey30') + 
  scale_x_discrete('# of fractions', labels = ~ replace(., . == '5', '5+')) +
  scale_y_continuous('Functional score') +
  boxed_theme() +
  theme(aspect.ratio = 1.6)
p2

## number of experiments vs. functional score
func1 = func %>% 
  group_by(gene, uniprot, position, functional_score) %>% 
  summarise(n_expts = sum(!is.na(gene))) %>% 
  ungroup() %>% 
  mutate(n_expts = winsorize(n_expts, c(NA, 5)))
p3 = func1 %>% 
  ggplot(aes(x = factor(n_expts), y = functional_score)) +
  # geom_violin(color = NA, fill = 'grey75', width = 1) +
  geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.6, size = 0.3,
               fill = pals::kelly()[1], color = 'grey30') + 
  scale_x_discrete('# of experiments', labels = ~ replace(., . == '5', '5+')) +
  scale_y_continuous('Functional score') +
  boxed_theme() +
  theme(aspect.ratio = 1.6)
p3

# combine and save
p = p1 | p2 | p3
p
ggsave("fig/phosphorylation/functional-scores.pdf", p,
       width = 18, height = 5, units = "cm", useDingbats = FALSE)

## plot number of phosphosites in non-human species found in 5+ fractions
dat = readRDS("data/QC/phosphosites-localization-scores.rds") %>% 
  filter(n_fractions > 0) %>% 
  # recode one-off species
  mutate(species_bk = species,
         species = fct_recode(species,
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
                              ## version 2
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
pal = pals::stepped3()[16] %>% c(pals::kelly()[1])
ff = dat %>% 
  filter(species != 'H. sapiens') %>% 
  group_by(species, gene) %>% 
  summarise(n_fractions = sum(n_fractions)) %>% 
  ungroup() %>% 
  mutate(five_fractions = n_fractions >= 5) %>% 
  mutate(five_fractions = as.character(five_fractions) %>% 
           fct_recode('1-4' = 'FALSE', '5+' = 'TRUE') %>% 
           fct_relevel('5+', '1-4'))
labs = ff %>% 
  dplyr::count(species, five_fractions) %>% 
  spread(five_fractions, n)
p4 = ff %>% 
  ggplot(aes(x = species)) +
  geom_bar(aes(fill = five_fractions, color = five_fractions),
           width = 0.75, size = 0.2, color = 'grey20') + 
  geom_text(data = labs, aes(y = `1-4` + `5+`, label = `5+`), size = 1.75,
            vjust = 0, nudge_y = 700) +
  scale_x_discrete(labels = ~ gsub("\n", " ", .), expand = c(0, 0.75)) +
  scale_y_continuous(expression('# of phosphosites'~(10^3)),
                     limits = c(0, 21e3),
                     expand = c(0, 0),
                     labels = ~ . / 1e3) +
  scale_fill_manual('# of fractions', values = rev(pal),
                    limits = c('1-4', '5+'),
                    breaks = c('1-4', '5+')) +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = 'top',
        legend.key.size = unit(0.35, 'lines'),
        aspect.ratio = 0.75)
p4
ggsave("fig/phosphorylation/functional-scores-5-fractions.pdf", p4,
       width = 7, height = 5.25, units = "cm", useDingbats = FALSE)

# total count
regulatory = filter(ff, five_fractions == '5+')
regulatory = dat %>% 
  filter(species != 'H. sapiens') %>% 
  group_by(species_bk, gene) %>% 
  summarise(n_fractions = sum(n_fractions)) %>% 
  ungroup() %>% 
  mutate(five_fractions = n_fractions >= 5) %>% 
  mutate(five_fractions = as.character(five_fractions) %>% 
           fct_recode('1-4' = 'FALSE', '5+' = 'TRUE') %>% 
           fct_relevel('5+', '1-4')) %>% 
  filter(five_fractions == '5+')
n_distinct(regulatory$gene)
n_distinct(regulatory$species_bk)
table(regulatory$species_bk) %>% mean
