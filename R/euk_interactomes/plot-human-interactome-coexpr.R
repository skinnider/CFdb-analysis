setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(flavin)
library(PrInCE)
source("R/theme.R")

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
v1 = read.delim("data/euk_interactomes/CF-MS-interactome.tsv.gz") %>% 
  mutate(version = 'Version 1')
v2 = readRDS("data/euk_interactomes/networks/network-species=Homo_sapiens-n_replicates=166-version=V2-n_features=1.rds") %>% 
  mutate(version = 'Version 2')

# plot RF against human high-throughput interactomes
hts_files = list.files("data/high_throughput", full.names = TRUE, 
                       pattern = '\\.gz$')
hts_ppis = map(hts_files, read.delim) %>%
  setNames(gsub("\\..*$", "", basename(hts_files)))

# set up coexpression files
coexpr_names = c('Kustatscher2019', 'Lapek2017')
coexpr_files = file.path('data/resources/coexpression', 
                         paste0(coexpr_names, '.txt.gz'))

# set up colocalization files
coloc_names = c('Geladaki2019', 'Orre2019')
coloc_files = file.path('data/resources/colocalization',
                        paste0(coloc_names, '.txt.gz'))

# read all matrices
mats1 = map(coexpr_files, ~ read.delim(.) %>% as.matrix()) %>%
  setNames(paste0('coexpression|', coexpr_names))
mats2 = map(coloc_files, ~ read.delim(.) %>% as.matrix()) %>%
  setNames(paste0('colocalization|', coloc_names))
mats = c(mats1, mats2)

# calculate coexpression
cors = map(mats, ~ cor(., use = 'p', method = 'pearson'))

# define plotting grid
grid = tidyr::crossing(matrix = names(cors), plot = c('v1_v2', 'hts'))
for (grid_idx in seq_len(nrow(grid))) {
  plot = grid$plot[grid_idx]
  matrix = grid$matrix[grid_idx]
  matrix_clean = chartr('|', '-', matrix) %>% gsub("protein/", "", .)
  
  # extract networks and label coexpressions
  if (plot == 'v1_v2') {
    v1_only = anti_join(v1, v2, by = c('protein_A', 'protein_B'))
    v2_only = anti_join(v2, v1, by = c('protein_A', 'protein_B'))
    ovr = inner_join(v1, v2, by = c('protein_A', 'protein_B'))
    df = bind_rows(
      v1_only %>% mutate(network = 'Version 1 only'),
      v2_only %>% mutate(network = 'Version 2 only'),
      ovr %>% mutate(network = 'Both')
    )
  } else if (plot == 'hts') {
    df = hts_ppis %>% 
      bind_rows(.id = 'network') %>% 
      # manually recode networks
      mutate(network = fct_recode(network,
                                  'HI-II-14' = 'Rolland2014',
                                  'HuRI' = 'Luck2019',
                                  'BioPlex' = 'Huttlin2015',
                                  'BioPlex 2' = 'Huttlin2017',
                                  'QUBIC' = 'Hein2015') %>% as.character()) %>% 
      dplyr::rename(protein_A = gene_A, protein_B = gene_B) %>% 
      bind_rows(mutate(v2, network = 'CF-MS'))
  } else {
    stop('not sure how to plot: `', plot, '` ...')
  }
  
  # tag correlations for each PPI
  cor = cors[[matrix]]
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
  p = plot_spectrum(cor_df, digits = digits)
  p
  ggsave(paste0("fig/euk_interactome/coexpression/human-", matrix_clean, 
                "-", plot, "-spectrum.pdf"),
         p, width = 4.25, height = 5, units = "cm", useDingbats = FALSE)
}
