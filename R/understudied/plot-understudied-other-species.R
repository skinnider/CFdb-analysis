# Plot the number of understudied proteins quantified by CF-MS in three 
# other species.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(data.table)
source("R/theme.R")

# read taxid map
taxid = read.csv("data/GO/species-map.csv")
filter(taxid, grepl('musculus|cerevisiae', species))

# read gene2pubmed
g2p = fread("data/resources/uncharacterized/gene2pubmed.gz") %>% 
  filter(`#tax_id` %in% c(10090, 559292, 3702)) %>% 
  mutate(species = `#tax_id` %>% 
           as.character() %>% 
           fct_recode('Mus musculus' = '10090',
                      'Saccharomyces cerevisiae' = '559292') %>% 
           as.character()) %>% 
  dplyr::count(species, GeneID)
# map to gene names
map_patt = paste0(c(10090, 559292, 3702), collapse = '|')
map_files = list.files("data/identifier", pattern = map_patt, full.names = TRUE)
maps = map(map_files, read_tsv, col_names = c('uniprot', 'db', 'id'))
map = bind_rows(maps)
gn = filter(map, db == 'Gene_Name') %>% dplyr::select(-db)
id = filter(map, db == 'GeneID') %>% dplyr::select(-db)
gn2id = left_join(gn, id, by = 'uniprot') %>%
  dplyr::select(starts_with('id')) %>%
  set_colnames(c('gene', 'id'))
# execute mapping
g2p %<>%
  dplyr::rename(id = GeneID, n_pubs = n) %>% 
  mutate(id = as.character(id)) %>% 
  left_join(gn2id, by = 'id') %>% 
  drop_na(gene) %>% 
  distinct(species, gene, n_pubs)

# read CDF
cdf = readRDS("data/QC/genes-CDF.rds") %>% 
  filter(grepl('musculus|cerevisiae', species))

# iterate through species
for (species in unique(cdf$species)) {
  # extract genes ever detected by CF-MS
  cdf0 = cdf %>% 
    filter(species == !!species)
  cf_genes = unique(cdf0$gene)
  
  # get clean species name
  sp_clean = fct_recode(species,
                        'mouse' = 'Mus musculus',
                        'yeast' = 'Saccharomyces cerevisiae') %>% 
    as.character()
  
  # plot histogram: number of publications per gene
  g2p0 = g2p %>% 
    filter(gene %in% cf_genes)
  lab = format(length(cf_genes), big.mark = ',') %>% 
    paste0(' ', sp_clean, ' genes\ndetected by CF-MS') %>% 
    data.frame(x = Inf, y = Inf, label = .)
  pal = pals::kovesi.linear_gow_65_90_c35(100)
  p1 = g2p0 %>% 
    ggplot(aes(x = n_pubs)) +
    geom_rect(data = head(g2p0, 1),
              aes(xmin = 0, xmax = 10, ymin = -Inf, ymax = Inf),
              fill = pal[10], color = 'white', size = 0.3, alpha = 0.35) +
    geom_histogram(bins = 30, size = 0.3, color = 'grey76', fill = 'grey92') + 
    geom_label(data = lab, aes(x = x, y = y, label = label), hjust = 1, vjust = 1,
               size = 1.75, label.size = 0, label.padding = unit(0.6, 'lines'),
               lineheight = 1, fill = NA) +
    annotation_logticks(sides = 'b', short = unit(0.05, 'lines'),
                        mid = unit(0.1, 'lines'), long = unit(0.15, 'lines'),
                        size = 0.2) +
    scale_color_manual(values = col) +
    scale_fill_manual(values = col) +
    scale_x_log10("# of publications", labels = fancy_scientific) + 
    scale_y_continuous("Proteins", expand = expansion(c(0, 0.08))) +
    boxed_theme() +
    theme(legend.position = 'none',
          aspect.ratio = 0.7)
  p1
  
  # bar chart
  pal = pals::kovesi.linear_gow_65_90_c35(100)[c(10, 100)]
  ylim = switch(sp_clean,
                'mouse' = 11e3,
                'yeast' = 6e3)
  ymax = switch(sp_clean,
                'mouse' = 10e3,
                'yeast' = 5e3)
  ybrks = switch(sp_clean,
                 'mouse' = seq(0, 10, 2) * 1e3,
                 'yeast' = seq(0, 5, 1) * 1e3)
  p2 = g2p0 %>%
    # color
    replace_na(list(n_pubs = 9999)) %>% 
    mutate(color = ifelse(n_pubs <= 10, '1-10', '>10') %>% 
             fct_relevel('1-10', '>10')) %>% 
    ggplot(aes(x = '1', fill = color)) +
    geom_bar(color = 'grey20', size = 0.15) + 
    geom_segment(aes(x = 0, xend = 0, y = 0, yend = ymax),
                 color = 'grey50', size = 0.4) +
    scale_fill_manual('# of publications', values = pal) +
    scale_y_continuous(expression('Proteins'~(10^3)), breaks = ybrks,
                       limits = c(0, ylim), labels = ~ . / 1e3) +
    guides(fill = guide_legend(title.position = 'top', title.hjust = 0.5)) +
    coord_flip() +
    boxed_theme() +
    theme(panel.border = element_blank(),
          aspect.ratio = 0.2,
          legend.key.size = unit(0.3, 'lines'),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())
  p2
  
  # combine with p1
  p2_v2 = p2 + 
    coord_cartesian() +
    boxed_theme() +
    theme(panel.border = element_blank(),
          legend.position = 'right',
          legend.key.size = unit(0.3, 'lines'),
          aspect.ratio = 5,
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  p2_v2
  p1_comb = p1 | p2_v2
  p1_comb
  ggsave(paste0("fig/understudied/publications-per-gene-combined-", sp_clean, ".pdf"),
         p1_comb, width = 9, height = 4, units = "cm", useDingbats = FALSE)
}
