setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read human interactomes
hts_files = list.files("data/high_throughput", full.names = TRUE, 
                       pattern = '\\.gz$')
hts_ppis = map(hts_files, read.delim) %>%
  setNames(gsub("\\..*$", "", basename(hts_files)))
names(hts_ppis) %<>% fct_recode('BioPlex' = 'Huttlin2015',
                                'BioPlex 2' = 'Huttlin2017',
                                'BioPlex 3-HCT116' = 'Huttlin2021-HCT116',
                                'BioPlex 3-293T' = 'Huttlin2021-293T',
                                'BioPlex 2' = 'Huttlin2017',
                                'QUBIC' = 'Hein2015',
                                'HI-II-14' = 'Rolland2014',
                                'HuRI' = 'Luck2019')

# alphabetize all
hts_ppis %<>% map(~ {
  sorted = t(apply(.x[, 1:2], 1, sort))
  .x$gene_A = sorted[, 1]
  .x$gene_B = sorted[, 2]  
  return(.x)
})
map(hts_ppis, ~ table(.$gene_A < .$gene_B))

# calculate Jaccard indices
jaccard = tidyr::crossing(network1 = names(hts_ppis),
                          network2 = names(hts_ppis)) %>% 
  pmap_dfr(function(...) {
    current = tibble(...)
    net1 = hts_ppis[[current$network1]]
    net2 = hts_ppis[[current$network2]]
    intersect = inner_join(net1, net2) %>% nrow()
    union = full_join(net1, net2) %>% nrow()
    mutate(current, intersect = intersect, union = union)
  }) %>% 
  mutate(jaccard = intersect / union)

# cluster
mat = jaccard %>% 
  dplyr::select(-intersect, -union) %>% 
  spread(network2, jaccard) %>% 
  column_to_rownames('network1') %>% 
  as.matrix()
clust = hclust(dist(mat))
lvls = with(clust, labels[order])

# plot
jaccard %<>%
  mutate(label = round(jaccard, digits = 2),
         network1 = factor(network1, levels = lvls),
         network2 = factor(network2, levels = lvls))
text = filter(jaccard, network1 != network2)
p1 = jaccard %>% 
  ggplot(aes(x = network1, y = network2, fill = jaccard)) +
  geom_tile(color = 'white') +
  geom_text(data = text, aes(label = label), size = 1.75,
            # color = scico(n = 100)[100]
            ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_scico(name = 'Jaccard index    ', limits = c(0, 1), 
                   # breaks = c(0.05, 0.95), 
                   breaks = c(0, 1),
                   labels = c('0.0', '1.0')) +
  guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
  coord_fixed() + 
  boxed_theme() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'top',
        # legend.justification = 'bottom',
        legend.key.size = unit(0.18, 'lines'))
p1

# marginal plot: total number of interactions
n_ppis = data.frame(network = names(hts_ppis),
                    n_ppis = map_int(hts_ppis, nrow)) %>% 
  mutate(network = factor(network, levels = lvls))
p2 = n_ppis %>% 
  ggplot(aes(x = network, y = n_ppis)) +
  geom_col(width = 0.7, fill = 'grey88') +
  scale_y_continuous(expression(Interactions~(10^3)), 
                     expand = expansion(c(0, 0.15)),
                     labels = ~ . / 1e3) +
  coord_flip() + 
  boxed_theme(size_sm = 6) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = 1.75)
p2  

p = p1 + p2
ggsave("fig/human_hts/jaccard.pdf", p, width = 8.5, height = 8, units = "cm", 
       useDingbats = FALSE)

# separately, plot density plus ticks plus median
jaccard0 = filter(jaccard, as.character(network1) < as.character(network2))
median = median(jaccard0$jaccard)
mean = mean(jaccard0$jaccard)
p3 = jaccard0 %>% 
  ggplot(aes(x = jaccard)) +
  geom_density(color = 'black', fill = 'grey92', size = 0.2) +
  geom_errorbar(aes(ymin = 0, ymax = 1), width = 0, size = 0.25) +
  geom_vline(aes(xintercept = mean), linetype = 'dotted', size = 0.2) +
  geom_label(data = head(jaccard0, 1),
             aes(x = mean, y = Inf, label = format(mean, digits = 2) %>% 
                   paste0('mean = ', .)), 
             hjust = 0, vjust = 1, size = 1.75, fill = NA,
             label.size = NA, label.padding = unit(0.6, 'lines')) +
  scale_x_continuous('Jaccard index', expand = c(0, 0)) +
  scale_y_continuous('Density', expand = expansion(c(0, 0.15))) +
  boxed_theme() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = 0.9)
p3

p = p1 + p2 + p3 + plot_layout(nrow = 1)
p
ggsave("fig/human_hts/jaccard-plus-density.pdf", p, width = 12, height = 8,
       units = "cm", useDingbats = FALSE)
