setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(seqinr)
source("R/theme.R")

# read gene CDF
cdf = readRDS("data/QC/genes-CDF.rds")

# iterate through species
species_list = c('Homo sapiens', 'Mus musculus', 'Arabidopsis thaliana', 
                 'Saccharomyces cerevisiae')
for (species in species_list) {
  # filter to relevant species
  cdf0 = filter(cdf, species == !!species)
  # get clean name for plots
  species_clean = fct_recode(species, 
                             'human' = 'Homo sapiens',
                             'mouse' = 'Mus musculus',
                             'Arabidopsis' = 'Arabidopsis thaliana',
                             'yeast' = 'Saccharomyces cerevisiae') %>% 
    as.character()
  
  # read fasta
  fasta_file = file.path("~/git/CFdb-searches/data/fasta/filtered",
                         fct_recode(species,
                                    'UP000005640-H.sapiens' = 'Homo sapiens',
                                    'UP000000589-M.musculus' = 'Mus musculus',
                                    'UP000006548-A.thaliana' = 'Arabidopsis thaliana',
                                    'UP000002311-S.cerevisiae' = 'Saccharomyces cerevisiae') %>% 
                           as.character() %>% 
                           paste0('.fasta.gz'))
  fasta = read.fasta(fasta_file, as.string = TRUE, seqtype = "AA")
  genes = map_chr(fasta, ~ attr(., 'Annot')) %>% 
    extract(grepl('GN=', .)) %>% 
    map_chr(~ gsub("^.*GN=", "", .) %>% gsub(" .*$", "", .))
  genes_uniq = unique(genes)
  
  # count number of fractions per gene
  genes = cdf0 %>% 
    group_by(gene) %>% 
    summarise(n_fractions = sum(n_fractions)) %>% 
    ungroup() %>% 
    filter(n_fractions > 0)
  
  # add zeroes from complete proteome
  all_genes = data.frame(gene = genes_uniq) %>% 
    left_join(genes, by = 'gene') %>% 
    replace_na(list(n_fractions = 0))
  
  # plot 1: CDF of number of fractions in which each gene was detected
  p1 = genes %>% 
    arrange(desc(n_fractions)) %>% 
    mutate(gene = factor(gene, levels = unique(.$gene)),
           gene_idx = as.integer(gene)) %>% 
    ggplot(aes(x = gene_idx, y = n_fractions, group = '1')) +
    geom_path(size = 0.3) + 
    scale_x_continuous('Protein #') +
    scale_y_log10('# of fractions', labels = fancy_scientific,
                  expand = c(0.02, 0)) +
    annotation_logticks(sides = 'l', short = unit(0.06, 'lines'),
                        mid = unit(0.12, 'lines'), long = unit(0.18, 'lines'),
                        size = 0.2) +
    boxed_theme() +
    theme(aspect.ratio = 1)
  p1
  # inset pie chart
  title = round(100 * mean(all_genes$n_fractions > 0), digits = 1) %>% 
    paste0('% of ', species_clean, '\nproteins quantified')
  p1_inset = all_genes %>% 
    mutate(detected = n_fractions > 0,
           title = title) %>% 
    ggplot(aes(x = '1', fill = detected)) +
    facet_grid(~ title) +
    geom_bar(size = 0.15, color = 'grey20') + 
    coord_polar(theta = 'y') + 
    scale_fill_manual(values = c('TRUE' = 'grey60', 'FALSE' = 'white')) +
    scale_color_manual(values = c('TRUE' = 'grey50', 'FALSE' = 'grey50')) +
    boxed_theme() + 
    theme(plot.margin = margin(rep(0, 4)),
          plot.background = element_blank(),
          legend.position = 'none',
          legend.key.size = unit(0.35, 'lines'),
          panel.border = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())
  p1_inset
  # superimpose 
  p1_comb = p1 + inset_element(p1_inset, left = 0.55, bottom = 0.55, right = 1, top = 1,
                               align_to = 'panel')
  p1_comb
  ggsave(paste0("fig/QC/genes-cdf-", species_clean, ".pdf"), p1_comb, 
         width = 5.25, height = 5.4, units = "cm", useDingbats = FALSE)
  
  # plot 2: CDF of number of fractions in which each gene was detected, v1 vs. v2
  genes_v1 = cdf0 %>% 
    filter(version == 'V1') %>% 
    group_by(gene) %>% 
    summarise(n_fractions = sum(n_fractions)) %>% 
    ungroup() %>% 
    filter(n_fractions > 0) %>% 
    arrange(desc(n_fractions)) %>% 
    mutate(gene = factor(gene, levels = unique(.$gene)),
           gene_idx = as.integer(gene))
  genes_v2 = genes %>% 
    arrange(desc(n_fractions)) %>% 
    mutate(gene = factor(gene, levels = unique(.$gene)),
           gene_idx = as.integer(gene))
  genes_version = bind_rows(
    mutate(genes_v1, version = 'Version 1'),
    mutate(genes_v2, version = 'Version 2')
  )
  version_pal = c('grey88', pals::stepped()[3])
  p2 = genes_version %>% 
    ggplot(aes(x = gene_idx, y = n_fractions, color = version)) +
    geom_path(size = 0.3) + 
    scale_x_continuous('Protein #') +
    scale_y_log10('# of fractions', labels = fancy_scientific,
                  expand = c(0.02, 0)) +
    annotation_logticks(sides = 'l', short = unit(0.06, 'lines'),
                        mid = unit(0.12, 'lines'), long = unit(0.18, 'lines'),
                        size = 0.2) +
    scale_color_manual('', values = version_pal,
                       breaks = c('Version 1', 'Version 2')) +
    boxed_theme() +
    theme(aspect.ratio = 1,
          legend.key.width = unit(0.4, 'lines'),
          legend.key.height = unit(0.2, 'lines'))
  p2
  # inset pie chart
  n1 = nrow(genes_v1)
  n2 = nrow(genes_v2) - n1
  n0 = length(genes_uniq) - n2 - n1
  title = format(n2, big.mark = ',') %>% 
    paste('additional\nproteins detected')
  p2_inset = data.frame(n = c(n1, n2, n0),
                        fill = c('Version 1', 'Version 2', 'Undetected')) %>% 
    mutate(title = title) %>% 
    ggplot(aes(x = '1', y = n, fill = fill)) +
    facet_grid(~ title) +
    geom_col(size = 0.15, color = 'grey20') + 
    coord_polar(theta = 'y') + 
    scale_fill_manual('', values = c('Version 1' = version_pal[1],
                                     'Version 2' = version_pal[2],
                                     'Undetected' = 'white')) +
    # guides(fill = guide_legend())
    boxed_theme() + 
    theme(plot.margin = margin(rep(0, 4)),
          plot.background = element_blank(),
          legend.position = 'none',
          legend.key.size = unit(0.35, 'lines'),
          panel.border = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())
  p2_inset
  # superimpose 
  p2_comb = p2 + inset_element(p2_inset, left = 0.55, bottom = 0.55, right = 1, top = 1,
                               align_to = 'panel')
  p2_comb
  ggsave(paste0("fig/QC/genes-cdf-", species_clean, "-version.pdf"), p2_comb, 
         width = 5.25, height = 5.4, units = "cm", useDingbats = FALSE)
  
  # plot 3: fold-change in any given # of fractions
  fold_change = map_dfr(c(1, seq(10, max(genes_version$n_fractions), 10)), ~ {
    n1 = filter(genes_version, version == 'Version 1', n_fractions >= .x) %$% 
      n_distinct(gene)
    n2 = filter(genes_version, version == 'Version 2', n_fractions >= .x) %$% 
      n_distinct(gene)
    fc = n2 / n1
    data.frame(n_fractions = .x, n1 = n1, n2 = n2, fc = fc)
  }) %>% 
    mutate(pct_increase = (n2 - n1) / n1)
  lab = filter(fold_change, n_fractions == 500) %>% 
    mutate(label = round(100 * pct_increase, digits = 1) %>% 
             paste0('% more proteins\nin ', n_fractions, '+ fractions'))
  p3 = fold_change %>% 
    filter(n_fractions <= 1000) %>% 
    ggplot(aes(x = n_fractions, y = pct_increase)) +
    geom_path(size = 0.3) + 
    geom_point(data = lab, size = 0.5) +
    geom_label_repel(data = lab, aes(label = label), size = 1.75,
                     fill = NA, min.segment.length = 0.01, segment.size = 0.15,
                     label.size = NA, lineheight = 0.9, nudge_y = 0.8,
                     nudge_x = -250) +
    scale_x_continuous('Min. # of fractions') +
    scale_y_continuous('% increase', labels = ~ . * 100, limits = c(0, NA),
                       expand = c(0.02, 0)) +
    # scale_y_log10(limits = c(1, NA)) +
    boxed_theme() +
    theme(aspect.ratio = 0.8)
  p3
  ggsave(paste0("fig/QC/genes-cdf-", species_clean, "-FC.pdf"), p3, 
         width = 4.75, height = 5.4, units = "cm", useDingbats = FALSE)
}
