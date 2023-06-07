# Plot the abundance of proteins detected in v1 vs. v2. 
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read PaxDB files
paxdb_files = list.files("data/resources/PaxDB", 
                         pattern = '*WHOLE_ORGANISM*', full.names = T)
paxdbs = map(paxdb_files, read.delim, comment.char = '#',
             header = F, col.names = c('paxid', 'id', 'abundance')) %>%
  setNames(c("Mus musculus", "Arabidopsis thaliana", "Homo sapiens"))

# read UniProt ID maps
id_files = paste0("data/identifier/", 
                  c("MOUSE_10090", "ARATH_3702", "HUMAN_9606"), 
                  "_idmapping.dat.gz")
maps = map(id_files, read_tsv, col_names = c("uniprot", "db", "id")) %>%
  setNames(c("Mus musculus", "Arabidopsis thaliana", "Homo sapiens"))

# map PaxDB to UniProt
species = names(maps)
abundance = map_dfr(species, ~ {
  paxdb = paxdbs[[.x]]
  map = maps[[.x]]
  # Ensembl to UniProt
  paxdb %<>%
    mutate(id = gsub("^(3702|10090|9606)\\.", "", id)) %>%
    left_join(map, by = 'id') %>%
    drop_na(uniprot, abundance) %>%
    # keep only one protein ID per gene
    group_by(paxid) %>%
    arrange(desc(abundance)) %>%
    dplyr::slice(1) %>%
    ungroup() %>% 
    distinct(uniprot, abundance)
  # UniProt to gene name
  gn = filter(map, db == 'Gene_Name') %>% 
    distinct(uniprot, id) %>% 
    dplyr::rename(gene = id)
  paxdb %<>%
    left_join(gn, by = 'uniprot') %>% 
    drop_na(gene) %>% 
    distinct(gene, abundance) %>% 
    mutate(species = .x)
})

# read gene CDF
cdf = readRDS("data/QC/genes-CDF.rds")

# extract v1 vs v2 genes
genes_version = cdf %>% 
  group_by(species, gene) %>% 
  arrange(version) %>% 
  dplyr::slice(1) %>% 
  ungroup()
# merge in abundance
abundance %<>%
  left_join(genes_version, by = c('species', 'gene')) %>% 
  replace_na(list(version = 'Never detected'))

# plot
version_pal = c('grey75', pals::stepped()[3]) %>% c('grey94')
p1 = abundance %>% 
  mutate(version = fct_recode(version, 
                              'Version 1' = 'V1', 
                              'Version 2 only' = 'V2') %>% 
           fct_relevel('Version 1', 'Version 2 only', 'Undetected'),
         species = fct_relevel(species, 'Homo sapiens', 'Mus musculus', 
                               'Arabidopsis thaliana')) %>% 
  ggplot(aes(x = version, y = abundance, color = version, fill = version)) +
  facet_grid(~ species) +
  geom_violin(alpha = 0.3, color = NA, width = 0.65) +
  geom_boxplot(outlier.shape = NA, size = 0.3, width = 0.65, alpha = 0.4) + 
  scale_y_log10('Protein abundance (ppm)', labels = fancy_scientific,
                expand = c(0.02, 0)) +
  annotation_logticks(sides = 'l', short = unit(0.05, 'lines'),
                      mid = unit(0.1, 'lines'), long = unit(0.15, 'lines'),
                      size = 0.2) +
  scale_color_manual('Detected in:', values = version_pal %>% 
                       replace(. == 'grey94', 'grey65')) +
  scale_fill_manual('Detected in:', values = version_pal) +
  boxed_theme() +
  theme(aspect.ratio = 1.4,
        # legend.position = 'none',
        legend.key.width = unit(0.35, 'lines'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())
p1
ggsave("fig/abundance/PaxDb-v1-v2-boxplot.pdf", p1,
       width = 8, height = 4.85, units = "cm", useDingbats = FALSE)

# plot number of experiments vs. abundance
abundance = map_dfr(species, ~ {
  paxdb = paxdbs[[.x]]
  map = maps[[.x]]
  # Ensembl to UniProt
  paxdb %<>%
    mutate(id = gsub("^(3702|10090|9606)\\.", "", id)) %>%
    left_join(map, by = 'id') %>%
    drop_na(uniprot, abundance) %>%
    # keep only one protein ID per gene
    group_by(paxid) %>%
    arrange(desc(abundance)) %>%
    dplyr::slice(1) %>%
    ungroup() %>% 
    distinct(uniprot, abundance)
  # UniProt to gene name
  gn = filter(map, db == 'Gene_Name') %>% 
    distinct(uniprot, id) %>% 
    dplyr::rename(gene = id)
  paxdb %<>%
    left_join(gn, by = 'uniprot') %>% 
    drop_na(gene) %>% 
    distinct(gene, abundance) %>% 
    mutate(species = .x)
})
expts_abu = cdf %>% 
  dplyr::count(species, gene) %>% 
  left_join(abundance, by = c('species', 'gene')) %>% 
  drop_na(abundance) %>% 
  filter(abundance > 0) %>% 
  group_by(species) %>% 
  mutate(density = get_density(n, log(abundance))) %>% 
  ungroup() %>% 
  arrange(density) 
pal = pals::ocean.deep(100)
p2 = map(c('Homo sapiens', 'Mus musculus', 'Arabidopsis thaliana'), ~ {
  expts_abu0 = filter(expts_abu, species == .x)
  range = range(expts_abu0$density)
  brks = c(range[1] + 0.1 * diff(range),
           range[2] - 0.1 * diff(range))
  p = expts_abu0 %>% 
    ggplot(aes(x = n, y = abundance, fill = density, color = density)) +
    facet_wrap(~ species, scales = 'free') +
    geom_point(shape = 21, size = 0.25, stroke = 0) +
    scale_x_continuous('# of experiments', expand = c(0.02, 0)) +
    scale_y_log10('Protein abundance (ppm)', labels = fancy_scientific,
                  expand = c(0.02, 0)) +
    annotation_logticks(sides = 'l', short = unit(0.05, 'lines'),
                        mid = unit(0.1, 'lines'), long = unit(0.15, 'lines'),
                        size = 0.2) +
    scale_fill_gradientn('Density', colours = pal, breaks = brks, 
                         labels = c('min', 'max')) +
    scale_color_gradientn('Density', colours = pal, breaks = brks, 
                          labels = c('min', 'max')) +
    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
    boxed_theme() +
    theme(aspect.ratio = 1,
          legend.position = 'right',
          legend.justification = 'bottom',
          legend.key.size = unit(0.18, 'lines'))
  if (.x != 'Arabidopsis thaliana')
    p = p + theme(legend.position = 'none')
  p
}) %>% wrap_plots(nrow = 1)
p2
ggsave("fig/abundance/PaxDb-vs-n-experiments.pdf", p2,
       width = 12.5, height = 4.85, units = "cm", useDingbats = FALSE)

# plot number of fractions vs. abundance
expts_fractions = cdf %>% 
  group_by(species, gene) %>% 
  summarise(n = sum(n_fractions)) %>% 
  ungroup() %>% 
  left_join(abundance, by = c('species', 'gene')) %>% 
  drop_na(abundance) %>% 
  filter(abundance > 0) %>% 
  group_by(species) %>% 
  mutate(density = get_density(n, log(abundance))) %>% 
  ungroup() %>% 
  arrange(density) 
pal = pals::ocean.deep(100)
p3 = map(c('Homo sapiens', 'Mus musculus', 'Arabidopsis thaliana'), ~ {
  expts_fractions0 = filter(expts_fractions, species == .x)
  quantile = quantile(expts_fractions0$n, probs = 0.99)
  expts_fractions0 %<>% mutate(n = winsorize(n, limits = c(0, quantile)))
  range = range(expts_fractions0$density)
  brks = c(range[1] + 0.1 * diff(range),
           range[2] - 0.1 * diff(range))
  p = expts_fractions0 %>% 
    ggplot(aes(x = n, y = abundance, fill = density, color = density)) +
    facet_wrap(~ species, scales = 'free') +
    geom_point(shape = 21, size = 0.25, stroke = 0) +
    scale_x_continuous('# of fractions', expand = c(0.02, 0)) +
    scale_y_log10('Protein abundance (ppm)', labels = fancy_scientific,
                  expand = c(0.02, 0)) +
    annotation_logticks(sides = 'l', short = unit(0.05, 'lines'),
                        mid = unit(0.1, 'lines'), long = unit(0.15, 'lines'),
                        size = 0.2) +
    scale_fill_gradientn('Density', colours = pal, breaks = brks, 
                         labels = c('min', 'max')) +
    scale_color_gradientn('Density', colours = pal, breaks = brks, 
                          labels = c('min', 'max')) +
    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
    boxed_theme() +
    theme(aspect.ratio = 1,
          legend.position = 'right',
          legend.justification = 'bottom',
          legend.key.size = unit(0.18, 'lines'))
  if (.x != 'Arabidopsis thaliana')
    p = p + theme(legend.position = 'none')
  p
}) %>% wrap_plots(nrow = 1)
p3
ggsave("fig/abundance/PaxDb-vs-n-fractions.pdf", p3,
       width = 12.5, height = 4.85, units = "cm", useDingbats = FALSE)
