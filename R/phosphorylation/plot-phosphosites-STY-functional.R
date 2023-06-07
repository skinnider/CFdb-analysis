# Plot the functional scores of human phospho-STY sites.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(readxl)
source("R/theme.R")

# read localization scores
dat = readRDS("data/QC/phosphosites-localization-scores.rds") %>% 
  filter(n_fractions > 0) %>% 
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

# filter to human
hs = filter(dat, species == 'H. sapiens') %>% 
  # separate to merge
  mutate(uniprot = strsplit(gene, '-') %>% 
           map_chr(~ head(., -2) %>% paste0(collapse = '-')),
         position = strsplit(gene, '-') %>% 
           map_chr(~ tail(., 1)))

# read Ochoa dataset
func = read_excel("data/resources/Ochoa/41587_2019_344_MOESM5_ESM.xlsx") %>% 
  mutate(position = as.character(position))

# merge
func %<>% left_join(hs, by = c('uniprot', 'position'))

## boxplot: functional scores of pY, pS, and pT detected by CF-MS
pal = pals::stepped3()[c(1, 5) + 2] %>% c(pals::kelly()[1])
col = rep(pals::kelly()[2], 3)
p = func %>% 
  filter(n_fractions > 0) %>% 
  ggplot(aes(x = reorder(sty, functional_score, stats::median), 
             y = functional_score,
             fill = sty, color = sty)) +
  geom_violin(alpha = 0.25, color = NA) +
  geom_boxplot(alpha = 0.4, size = 0.25, width = 0.6, outlier.shape = NA) +
  scale_x_discrete('Residue', labels = ~ paste0('p', .)) +
  scale_y_continuous('Functional score') +
  scale_fill_manual('', values = pal) +
  scale_color_manual('', values = col) +
  boxed_theme() +
  theme(aspect.ratio = 1.4,
        legend.position = 'none')
p
ggsave("fig/phosphorylation/STY-functional-scores.pdf", p,
       width = 4.5, height = 4.5, units = "cm", useDingbats = FALSE)

## Sankey plot: ever detected vs. never detected

# first, merge in acceptor residues
phos = read_excel("data/resources/Ochoa/41587_2019_344_MOESM4_ESM.xlsx")
residues = dplyr::select(phos, uniprot, position, residue) %>% 
  mutate(position = as.character(position))
func %<>% left_join(residues, by = c('uniprot', 'position'))
table(func$residue == func$sty) ## confirmed

# plot Sankey: ever vs. never detected
props_overall = func %>% 
  dplyr::count(residue)
props_cf = dat %>% 
  distinct(gene, sty) %>% 
  dplyr::rename(residue = sty) %>% 
  dplyr::count(residue)
props_5 = dat %>% 
  filter(n_fractions >= 5) %>% 
  distinct(gene, sty) %>% 
  dplyr::rename(residue = sty) %>% 
  dplyr::count(residue)
props = bind_rows(
  props_overall %>% mutate(label = 'Overall'),
  props_cf %>% mutate(label = 'Detected by CF-MS'),
  props_5 %>% mutate(label = 'Detected in 5+ fractions')
) %>% 
  group_by(label) %>% 
  mutate(prop = n / sum(n)) %>% 
  ungroup() %>% 
  mutate(label = factor(label, levels = c('Overall', 'Detected by CF-MS',
                                          'Detected in 5+ fractions')),
         xval = as.integer(label))

# plot
wide = props %>% 
  dplyr::select(label, residue, prop) %>% 
  spread(label, prop) %>% 
  mutate_all(~ replace(., is.na(.), 0))
colnames = colnames(wide)[-1] %>% head(-1)
polygons = map_dfr(seq_along(colnames), ~ {
  colname1 = colnames(wide)[-1][.x]
  colname2 = colnames(wide)[-1][.x + 1]
  print(colname1)
  print(colname2)
  bind_rows(
    mutate(wide, id = residue,
           x = .x + 0.25, y = c(0, cumsum(wide[[colname1]]) %>% head(-1))),
    mutate(wide, id = residue,
           x = .x + 0.25, y = c(cumsum(wide[[colname1]]))),
    mutate(wide, id = residue,
           x = .x + 0.75, y = c(cumsum(wide[[colname2]]))),
    mutate(wide, id = residue,
           x = .x + 0.75, y = c(0, cumsum(wide[[colname2]]) %>% head(-1)))
  ) %>% 
    mutate(group = as.character(.x))
}) %>% 
  mutate(y = 1 - y)
p = polygons %>% 
  ggplot(aes(x = x, y = y, value = residue, fill = residue)) +
  geom_polygon(alpha = 0.4, color = 'white', size = 0.17, 
               aes(group = paste(residue, group))) +
  geom_col(data = props,
           aes(x = xval, y = prop),
           color = 'grey15', size = 0.15, width = 0.5) +
  geom_segment(aes(y = 0, yend = 1, x = -Inf, xend = -Inf), size = 0.4,
               color = 'grey50') +
  scale_x_continuous(breaks = seq_len(3), labels = colnames(wide)[-1],
                     limits = c(0.5, 7.5)) +
  scale_y_continuous('% of phosphosites', labels = function(x) x * 100,
                     limits = c(0, 1)) +
  scale_fill_manual(values = pal, name = '') +
  clean_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(-0.2, 'lines'),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        aspect.ratio = 1.2,
        legend.key.size = unit(0.3, 'lines'),
        legend.position = 'right')
p
ggsave("fig/phosphorylation/STY-sankey-ever-detected.pdf", p,
       width = 7.5, height = 4.5, units = "cm", useDingbats = FALSE)

# test
x1 = func$residue == 'Y'
y1 = rep('Proteome', nrow(func))
x2 = dat %>% 
  distinct(gene, sty) %>% 
  pull(sty) %>% 
  `==`('Y')
y2 = rep('CF-MS', length(x2))
chisq.test(c(x1, x2), c(y1, y2))
