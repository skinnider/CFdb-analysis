setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

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
    # scale_fill_paletteer_c("pals::kovesi.diverging_rainbow_bgymr_45_85_c67",
    #                        name = 'AUC',
    #                        limits = range,
    #                        breaks = range,
    #                        labels = range,
    #                        direction = -1
    # ) +
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
          aspect.ratio = 0.15)
  p1a
  
  # median as the third column
  p1b = spectrum_neg %>%
    mutate(network = factor(network, levels = lvls)) %>%
    arrange(network) %>%
    ggplot(aes(y = network, x = 0)) +
    facet_grid(network ~ ., scales = 'free') +
    geom_text(aes(label = format(median, digits = 2)), 
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

# read EGAD results
egad_v1 = readRDS("data/euk_interactomes/EGAD/network-species=Homo_sapiens-n_replicates=46-version=V1-n_features=1.rds")
egad_v2 = readRDS("data/euk_interactomes/EGAD/network-species=Homo_sapiens-n_replicates=166-version=V2-n_features=1.rds")
sizes_v1 = attr(egad_v1, 'gene_set_sizes')
sizes_v2 = attr(egad_v2, 'gene_set_sizes')

# read HTS results
hts_names = c('Hein2015', 'Huttlin2015', 'Huttlin2017', 'Huttlin2021-293T',
              'Huttlin2021-HCT116', 'Rolland2014', 'Luck2019')
hts_files = paste0("data/euk_interactomes/EGAD/", hts_names, ".rds")
hts = map(hts_files, readRDS) %>% setNames(hts_names)

# define plotting grid
grid = tidyr::crossing(terms = c('GO', 'disease_genes'),
                       plot = c('v1_v2_overlapping', 'hts'))
for (grid_idx in seq_len(nrow(grid))) {
  terms = grid$terms[grid_idx]
  plot = grid$plot[grid_idx]
  
  # extract data and filter breadth
  if (plot == 'v1_v2') {
    df = bind_rows(
      egad_v1 %>% filter(between(n_proteins, 10, 100)) %>% 
        mutate(network = 'Version 1'),
      egad_v2 %>% filter(between(n_proteins, 10, 100)) %>% 
        mutate(network = 'Version 2')
    )
    # overlapping GO terms only
    df %<>%
      group_by(term) %>% 
      filter(n_distinct(network) == 2) %>% 
      ungroup()
  } else if (plot == 'hts') {
    v2 = egad_v2 %>% filter(between(n_proteins, 10, 100)) %>% 
      mutate(network = 'CF-MS')
    df = map(hts, ~ filter(.x, between(n_proteins, 10, 100)) %>% 
               dplyr::select(-network)) %>% 
      bind_rows(.id = 'network') %>% 
      bind_rows(v2)
  } else {
    stop('not sure how to plot: `', plot, '` ...')
  }
  
  # filter to relevant category
  df %<>% mutate(category = gsub("\\|.*$", "", term))
  if (terms == 'GO') {
    df %<>% filter(category == 'GO')
  } else if (terms == 'disease_genes') {
    df %<>% filter(category == 'disease_gene')
  } else {
    stop('not sure how to handle category: `', terms, '` ...')
  }
  
  # post-processing
  df %<>% 
    # clean up network names
    mutate(network = fct_recode(network,
                                'HI-II-14' = 'Rolland2014',
                                'HuRI' = 'Luck2019',
                                'BioPlex' = 'Huttlin2015',
                                'BioPlex 2' = 'Huttlin2017',
                                'BioPlex 3-HCT116' = 'Huttlin2021-HCT116',
                                'BioPlex 3-293T' = 'Huttlin2021-293T',
                                'QUBIC' = 'Hein2015') %>% as.character())
  
  # plot
  p = plot_spectrum(df)
  # p
  # save
  output_file = paste0("fig/euk_interactomes/EGAD/human-plot=", plot, 
                       "-terms=", terms, ".pdf")
  ggsave(output_file, p, width = 4.25, height = 5, 
         units = "cm", useDingbats = FALSE)
} 
