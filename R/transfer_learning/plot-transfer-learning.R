setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(PrInCE)
source("R/theme.R")

## PR curves
# list files
files = list.files("data/transfer_learning/networks_50",
                   pattern = 'rds', full.names = TRUE)
nets = map(files, readRDS) %>% 
  map(~ mutate(.x, idx = row_number())) %>% 
  setNames(basename(files)) %>% 
  bind_rows(.id = 'filename') %>% 
  mutate(prop_augmentation = gsub("^.*augmentation=|\\.rds$", "", filename)) %>% 
  dplyr::select(-filename) %>% 
  type_convert()

# focus on default
pal = c('grey80', brewer.pal(5, 'Spectral'))
max = max(nets$idx)
p1 = nets %>% 
  filter(prop_augmentation %in% c(0, 0.33)) %>% 
  ggplot(aes(x = idx, y = precision, color = factor(prop_augmentation))) +
  geom_path(size = 0.35) + 
  scale_color_manual('', values = pal,
                     labels = c('Default', 'Transfer learning')) +
  scale_y_continuous('Precision') +
  scale_x_continuous(expression('# of PPIs'~(10^3)), labels = ~ . / 1e3) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.key.width = unit(0.4, 'lines'),
        legend.key.height = unit(0.4, 'lines'),
        legend.position = c(0.97, 1.04),
        legend.justification = c(1, 1)) +
  ggh4x::force_panelsizes(rows = unit(30, 'mm'),
                          cols = unit(30, 'mm'))
p1
ggsave("fig/transfer_learning/PR-curves.pdf", p1,
       width = 6, height = 6, units = "cm", useDingbats = FALSE)

# facet by degree of augmentation
pal = c('grey80', brewer.pal(5, 'Spectral'))
props = unique(nets$prop_augmentation) %>% 
  extract(between(., 0.01, 0.51)) %>% 
  setdiff(0.2)
xlim = filter(nets, prop_augmentation <= 0.5) %$% max(idx)
plots = map(props, ~ {
  prop = .
  prop_idx = which(props == prop)
  pal = c('grey75', brewer.pal(5, 'Spectral')) # c(rev(brewer.pal(5, 'Spectral')))[prop_idx])
  p = nets %>% 
    mutate(title = paste0(100 * prop, '% augmentation')) %>% 
    filter(prop_augmentation %in% c(0, prop)) %>% 
    ggplot(aes(x = idx, y = precision, color = factor(prop_augmentation))) +
    facet_grid(~ title) +
    geom_path(size = 0.35) + 
    scale_color_manual('', values = pal,
                       labels = c('Default', 'Transfer learning')) +
    scale_y_continuous('Precision') +
    scale_x_continuous(expression('# of PPIs'~(10^3)), labels = ~ . / 1e3,
                       limits = c(0, xlim)) +
    boxed_theme() +
    theme(aspect.ratio = 1,
          legend.key.width = unit(0.4, 'lines'),
          legend.key.height = unit(0.4, 'lines'),
          legend.position = c(0.96, 1.06),
          legend.justification = c(1, 1))
  if (prop != props[1])
    p = p + theme(axis.text.y = element_blank(),
                  axis.title.y = element_blank(),
                  axis.ticks.y = element_blank())
  if (prop != dplyr::last(props))
    p = p + theme(legend.position = 'none')
  p
})
p2 = wrap_plots(plots, nrow = 1)
p2
ggsave("fig/transfer_learning/PR-curves-augmentation.pdf", p2,
       width = 18, height = 4.15, units = "cm", useDingbats = FALSE)

## coexpression
# set up function
plot_spectrum = function(dat, digits = 2, range = NULL) {
  spectrum = dat %>%
    # tag y-value of each
    group_by(network) %>%
    arrange(correlation) %>%
    mutate(y = row_number(),
           y = y / max(y),
           y = rescale(y, c(0, 1))) %>%
    ungroup()
  spectrum_neg = spectrum %>%
    group_by(network) %>%
    summarise(y = mean(correlation < 0, na.rm = TRUE), 
              median = median(correlation, na.rm = TRUE),
              mean = mean(correlation, na.rm = TRUE)) %>%
    arrange(y)
  
  max_dist = 3
  lvls = spectrum_neg$network
  
  if (is.null(range))
    range = range(spectrum$correlation) %>% round(2)
  
  p1a = spectrum %>%
    mutate(network = factor(network, levels = lvls)) %>%
    mutate(correlation = winsorize(correlation, range)) %>%
    ggplot(aes(x = network, y = y)) +
    facet_grid(network ~ ., scales = 'free') +
    geom_raster(aes(fill = correlation)) +
    geom_errorbar(data = spectrum_neg %>%
                    mutate(network = factor(network, levels = lvls)),
                  aes(ymin = y, ymax = y), width = 1, size = 0.3,
                  color = 'grey20') +
    scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                           name = 'Correlation',
                           limits = range,
                           breaks = range,
                           labels = range
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous('% of interactions', 
                       expand = c(0, 0), 
                       limits = c(-0.001, 1),
                       labels = function(x) paste0(x * 100, '%'),
                       breaks = seq(0, 1, 0.25)
    ) +
    guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE,
                                 title.pos = 'top', title.hjust = 0.5)) +
    coord_flip() +
    boxed_theme() +
    theme(axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          strip.text = element_blank(),
          panel.spacing = unit(0.1, 'lines'),
          plot.margin = margin(rep(0, 4)),
          legend.key.width = unit(0.25, 'lines'),
          legend.key.height = unit(0.25, 'lines'),
          legend.position = 'top',
          aspect.ratio = 0.15)
  p1a
  
  # median as the third column
  p1b = spectrum_neg %>%
    mutate(network = factor(network, levels = lvls)) %>%
    arrange(network) %>%
    ggplot(aes(y = network, x = 0)) +
    facet_grid(network ~ ., scales = 'free') +
    geom_text(aes(label = format(median, digits = digits)), 
              size = 2, hjust = 0, color = 'grey30') +
    # geom_tile(color = 'white') +
    scale_y_discrete(labels = function(x) gsub("^.*\\|", "", x),
                     expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 5)) +
    # scale_fill_manual(values = pal1) +
    # coord_fixed() +
    boxed_theme() +
    theme(strip.text = element_blank(),
          panel.spacing = unit(0.1, 'lines'),
          panel.border = element_blank(),
          plot.background = element_blank(),
          plot.margin = margin(rep(0, 4)),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          aspect.ratio = 0.75,
          legend.position = 'none')
  p1b
  
  p = p1a + p1b + plot_layout(widths = c(1, 0.2))
  return(p)
}

# read networks
v1 = readRDS("data/transfer_learning/networks_50/network-prop_augmentation=0.rds") %>% 
  mutate(version = 'Default')
v2 = readRDS("data/transfer_learning/networks_50/network-prop_augmentation=0.33.rds") %>% 
  mutate(version = 'Transfer learning')

# read protein groups
mat = readRDS("data/bee/protein-groups-egg.rds") %>% 
  replace(., !is.finite(.), 0)

# correlation
cor = cor(t(mat), use = 'p', method = 'pearson')

# extract networks and label coexpressions
df = bind_rows(
  v1 %>% mutate(network = 'Default'),
  v2 %>% mutate(network = 'Transfer learning')
)

# tag correlations for each PPI
cor_df = filter(df, protein_A %in% rownames(cor) &
                  protein_B %in% rownames(cor))
cor_idxs = cor_df %>% 
  dplyr::select(protein_A, protein_B) %>% 
  as.matrix()
cor_df$correlation = cor[cor_idxs]
cor_df %<>% filter(is.finite(correlation))

## spectrum
medians = cor_df %>% 
  group_by(network) %>%
  summarise(median = median(correlation, na.rm = TRUE)) %>%
  pull(median)
digits = ifelse(any(medians < 0.095), 1, 2)
p3 = plot_spectrum(cor_df, digits = digits)
p3
ggsave("fig/transfer_learning/coexpression.pdf", p3, width = 4.75, height = 5, 
       units = "cm", useDingbats = FALSE)

## cofractionation in Drosophila
# read networks
v1_bee = readRDS("data/transfer_learning/networks_50/network-prop_augmentation=0.rds") %>% 
  mutate(version = 'Default')
v2_bee = readRDS("data/transfer_learning/networks_50/network-prop_augmentation=0.33.rds") %>% 
  mutate(version = 'Transfer learning')

# map both to bee
ortho = read.delim("data/bee/ortholog/InParanoid-bee-fly.txt")
map = with(ortho, setNames(target_gene, source_gene))
v1 = mutate(v1_bee, protein_A = map[protein_A], protein_B = map[protein_B]) %>% 
  drop_na(protein_A, protein_B)
v2 = mutate(v2_bee, protein_A = map[protein_A], protein_B = map[protein_B]) %>% 
  drop_na(protein_A, protein_B)

# read chromatograms
expts = read.csv("~/git/CFdb-searches/data/experiments.csv")
expts0 = filter(expts, grepl('Drosophila', Species))
chrom_files = with(expts0, file.path("~/git/CFdb-searches/data/chromatograms_gene", 
                                     Accession, Replicate, 'default', 'iBAQ.rds'))
chroms = map(chrom_files, readRDS) %>% 
  setNames(expts0$Replicate) %>% 
  map(~ extract(.x, rowSums(is.finite(.) & . > 0) >= 4, ) %>% 
        log() %>% 
        replace(., !is.finite(.), NA))

## also create a merged dataset
proteins = map(chroms, rownames) %>% Reduce(union, .)
fractions = map(chroms, colnames) %>% unlist()
merged = matrix(NA, nrow = length(proteins), ncol = length(fractions),
                dimnames = list(proteins, fractions))
for (chrom in chroms) 
  merged[rownames(chrom), colnames(chrom)] = chrom
chroms[['merged']] = merged

# calculate coexpression
cors = map(chroms, ~ cor(t(.), use = 'p', method = 'pearson'))

# plot all together
cor_df = map_dfr(names(cors)[1:4], ~ {
  df = bind_rows(
    v1 %>% mutate(network = 'Default'),
    v2 %>% mutate(network = 'Transfer learning')
  )
  cor = cors[[.x]]
  cor_df = filter(df, protein_A %in% rownames(cor) &
                    protein_B %in% rownames(cor))
  cor_idxs = cor_df %>% 
    dplyr::select(protein_A, protein_B) %>% 
    as.matrix()
  cor_df$correlation = cor[cor_idxs]
  cor_df %<>% filter(is.finite(correlation))
  return(cor_df)
})
## spectrum
medians = cor_df %>% 
  group_by(network) %>%
  summarise(median = median(correlation, na.rm = TRUE)) %>%
  pull(median)
digits = ifelse(any(medians < 0.095), 1, 2)
p4 = plot_spectrum(cor_df, digits = digits)
p4
ggsave("fig/transfer_learning/cofractionation.pdf",
       p4, width = 4.75, height = 5, units = "cm", useDingbats = FALSE)

## GO
# set up function
plot_spectrum = function(dat, range = NULL) {
  spectrum = dat %>%
    # tag y-value of each
    group_by(network) %>%
    arrange(auroc) %>%
    mutate(y = row_number(),
           y = y / max(y),
           y = rescale(y, c(0, 1))) %>%
    ungroup()
  spectrum_neg = spectrum %>%
    group_by(network) %>%
    summarise(y = mean(auroc < 0.5), median = median(auroc),
              mean = mean(auroc)) %>%
    arrange(y)
  
  max_dist = 3
  digits = 3
  lvls = spectrum_neg$network
  
  if (is.null(range))
    range = range(spectrum$auroc) %>% round(2)
  
  p1a = spectrum %>%
    mutate(network = factor(network, levels = lvls)) %>%
    mutate(auroc = winsorize(auroc, range)) %>%
    ggplot(aes(x = network, y = y)) +
    facet_grid(network ~ ., scales = 'free') +
    geom_raster(aes(fill = auroc)) +
    geom_errorbar(data = spectrum_neg %>%
                    mutate(network = factor(network, levels = lvls)),
                  aes(ymin = y, ymax = y), width = 1, size = 0.3,
                  color = 'grey20') +
    scale_fill_gradientn(colors = jdb_palette("brewer_spectra"),
                         name = 'AUC',
                         limits = range,
                         breaks = range,
                         labels = range) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous('% of GO terms', 
                       expand = c(0, 0), 
                       limits = c(-0.001, 1),
                       labels = function(x) paste0(x * 100, '%'),
                       breaks = seq(0, 1, 0.25)
    ) +
    guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE,
                                 title.pos = 'top', title.hjust = 0.5)) +
    coord_flip() +
    boxed_theme() +
    theme(axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          strip.text = element_blank(),
          panel.spacing = unit(0.1, 'lines'),
          plot.margin = margin(rep(0, 4)),
          legend.key.width = unit(0.25, 'lines'),
          legend.key.height = unit(0.25, 'lines'),
          legend.position = 'top',
          aspect.ratio = 0.15) +
    ggh4x::force_panelsizes(cols = unit(25, 'mm'),
                            rows = unit(0.15 * 25, 'mm'))
  p1a
  
  # median as the third column
  p1b = spectrum_neg %>%
    mutate(network = factor(network, levels = lvls)) %>%
    arrange(network) %>%
    ggplot(aes(y = network, x = 0)) +
    facet_grid(network ~ ., scales = 'free') +
    geom_text(aes(label = format(median, digits = 2)), 
              size = 2, hjust = 0, color = 'grey30') +
    scale_y_discrete(labels = function(x) gsub("^.*\\|", "", x),
                     expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 5)) +
    boxed_theme() +
    theme(strip.text = element_blank(),
          panel.spacing = unit(0.1, 'lines'),
          panel.border = element_blank(),
          plot.background = element_blank(),
          plot.margin = margin(rep(0, 4)),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          aspect.ratio = 0.75,
          legend.position = 'none') +
    ggh4x::force_panelsizes(rows = unit(0.15 * 25, 'mm'))
  p1b
  
  p = p1a + p1b + plot_layout(widths = c(1, 0.2))
  return(p)
}

# read EGAD results
def = readRDS("data/transfer_learning/EGAD/network-prop_augmentation=0.rds") %>% 
  mutate(network = 'Default')
tl = readRDS("data/transfer_learning/EGAD/network-prop_augmentation=0.33.rds") %>% 
  mutate(network = 'Transfer learning')
dat = bind_rows(def, tl) %>% 
  # filter to between 10-100 proteins
  filter(between(n_proteins, 10, 100))
# plot
p = plot_spectrum(dat)
p
ggsave("fig/transfer_learning/EGAD-GO.pdf", p, 
       width = 5, height = 4, units = "cm", useDingbats = FALSE)
