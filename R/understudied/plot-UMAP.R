setwd("~/git/CFdb-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read gene2pubmed
g2p = read.delim("data/resources/uncharacterized/gene2pubmed.gz") %>% 
  filter(X.tax_id == 9606) %>% 
  dplyr::count(GeneID)
# map to gene names
map = read_tsv("data/identifier/HUMAN_9606_idmapping.dat.gz",
               col_names = c("uniprot", "db", "id"))
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
  distinct(gene, n_pubs)

# read human CDF
cdf = readRDS("data/QC/genes-CDF.rds")
cdf_hs = filter(cdf, species == 'Homo sapiens')
# extract genes ever detected by CF-MS
cf_genes = unique(cdf_hs$gene)

# plot histogram: number of publications per gene
g2p0 = g2p %>% 
  filter(gene %in% cf_genes)
lab = format(length(cf_genes), big.mark = ',') %>% 
  paste0(' human genes\ndetected by CF-MS') %>% 
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
  scale_y_continuous("Proteins", expand = c(0, 0), limits = c(0, 1550)) +
  boxed_theme() +
  theme(legend.position = 'none',
        aspect.ratio = 0.7)
p1
ggsave("fig/understudied/publications-per-gene.pdf", p1,
       width = 5, height = 5, units = "cm", useDingbats = FALSE)

# bar chart
pal = pals::kovesi.linear_gow_65_90_c35(100)[c(10, 100)]
p2 = g2p0 %>%
  # color
  replace_na(list(n_pubs = 9999)) %>% 
  mutate(color = ifelse(n_pubs <= 10, '1-10', '>10') %>% 
           fct_relevel('1-10', '>10')) %>% 
  ggplot(aes(x = '1', fill = color)) +
  geom_bar(color = 'grey20', size = 0.15) + 
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 12e3),
               color = 'grey50', size = 0.4) +
  scale_fill_manual('# of publications', values = pal) +
  scale_y_continuous(expression('Proteins'~(10^3)), expand = c(0, 200),
                     limits = c(0, 12e3), labels = ~ . / 1e3) +
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
ggsave("fig/understudied/publications-per-gene-bar-chart.pdf", p2,
       width = 4, height = 4, units = "cm", useDingbats = FALSE)

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
ggsave("fig/understudied/publications-per-gene-combined.pdf", p1_comb,
       width = 9, height = 4, units = "cm", useDingbats = FALSE)

# UMAP
umap = readRDS("data/UMAP/coords-human.rds")
# convert to data frame
umap %<>% as.data.frame() %>% 
  set_colnames(c('x', 'y')) %>% 
  rownames_to_column('gene')
# merge in number of proteins
umap %<>% left_join(g2p, by = 'gene') %>% 
  # color
  replace_na(list(n_pubs = 9999)) %>% 
  mutate(color = ifelse(n_pubs <= 10, '1-10', '>10') %>% 
           fct_relevel('1-10', '>10'))
# plot
# pal = pals::kovesi.linear_gow_65_90_c35(100)[c(10, 30, 50, 100)]
pal = pals::kovesi.linear_gow_65_90_c35(100)[c(10, 100)]
color = c('black', pal[2])
p3 = umap %>%
  arrange(desc(color)) %>% 
  ggplot(aes(x = x, y = y, color = color, fill = color)) + 
  geom_point(shape = 21, size = 0.6, stroke = 0.1) + 
  scale_color_manual('# of publications', values = color) +
  scale_fill_manual('# of publications', values = pal) +
  # scale_color_gradientn(colours = pal) +
  # guides(color = guide_colorbar(frame.colour = 'black', ticks = FALSE)) + 
  labs(x = "UMAP1", y = "UMAP2") +
  guides(fill = guide_legend(override.aes = list(size = 1.1, stroke = 0.2))) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(0, 'lines'),
        axis.ticks.length.y = unit(0, 'lines'),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(hjust = 0, vjust = 1),
        axis.title.y = element_text(hjust = 0, vjust = 0),
        legend.key.size = unit(0.35, 'lines'),
        legend.position = 'right')
p3
ggsave("fig/understudied/umap-uncharacterized.pdf", p3, 
       width = 6, height = 5.4, units = "cm", useDingbats = FALSE)

# last plot: which of these were v1 vs. v2?
uncharacterized = g2p0 %>% filter(n_pubs <= 10) %>% pull(gene)
uncharacterized_v = cdf_hs %>% 
  filter(gene %in% uncharacterized) %>% 
  group_by(gene) %>% 
  arrange(version) %>% 
  dplyr::slice(1) %>% 
  ungroup()
table(uncharacterized_v$version)
## plot on UMAP
version_pal = c('grey88', pals::stepped()[3], 
                pals::kovesi.linear_gow_65_90_c35(100)[100]) %>% 
  setNames(c('V1', 'V2', '>10'))
version_color = c('black', 'black', pal[2]) %>% 
  setNames(c('V1', 'V2', '>10'))
p4 = umap %>%
  left_join(uncharacterized_v, by = 'gene') %>% 
  mutate(color = ifelse(gene %in% uncharacterized, version, '>10') %>% 
           fct_relevel('>10', 'V1', 'V2')) %>% 
  arrange(color) %>% 
  ggplot(aes(x = x, y = y, color = color, fill = color)) + 
  geom_point(shape = 21, size = 0.6, stroke = 0.1) + 
  scale_color_manual('Uncharacterized\nproteins', values = version_color,
                     breaks = c('V1', 'V2'),
                     labels = c('Version 1', 'Version 2')) +
  scale_fill_manual('Uncharacterized\nproteins', values = version_pal,
                    breaks = c('V1', 'V2'),
                    labels = c('Version 1', 'Version 2')) +
  labs(x = "UMAP1", y = "UMAP2") +
  guides(fill = guide_legend(override.aes = list(size = 1.1, stroke = 0.2))) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(0, 'lines'),
        axis.ticks.length.y = unit(0, 'lines'),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(hjust = 0, vjust = 1),
        axis.title.y = element_text(hjust = 0, vjust = 0),
        legend.key.size = unit(0.35, 'lines'),
        legend.position = 'right')
p4
## inset pie chart
umap0 = umap %>%
  left_join(uncharacterized_v, by = 'gene') %>% 
  drop_na(version)
title = round(100 * mean(umap0$version == 'V2'), digits = 1) %>% 
  paste0('% of uncharacterized\nproteins in version 2')
p4_inset = umap0 %>% 
  mutate(title = title) %>% 
  ggplot(aes(x = '1', fill = version)) +
  facet_grid(~ title) +
  geom_bar(size = 0.15, color = 'grey20') + 
  coord_polar(theta = 'y') + 
  scale_fill_manual(values = version_pal) +
  scale_color_manual(values = c('V2' = 'grey50', 'V1' = 'grey50')) +
  boxed_theme() + 
  theme(plot.margin = margin(rep(0, 4)),
        panel.background = element_blank(),
        # plot.background = element_rect(fill='pink'),
        plot.background = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        legend.key.size = unit(0.35, 'lines'),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
p4_inset
# superimpose 
p4_comb = p4 + inset_element(p4_inset, left = 0, bottom = 0.55, right = 0.35, top = 1,
                             align_to = 'panel')
p4_comb
ggsave("fig/understudied/umap-uncharacterized-by-version.pdf", p4_comb, 
       width = 6, height = 5.4, units = "cm", useDingbats = FALSE)

# compute p-value
versions = cdf_hs %>% 
  group_by(gene) %>% 
  arrange(version) %>% 
  dplyr::slice(1) %>% 
  ungroup() %>% 
  distinct(gene, version)
g2p1 = g2p0 %>% left_join(versions, by = 'gene') 
vec1 = g2p1$n_pubs <= 10
vec2 = g2p1$version == 'V2'
chisq.test(vec1, vec2)
