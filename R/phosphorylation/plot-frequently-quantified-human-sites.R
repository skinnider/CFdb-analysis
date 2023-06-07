setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read data frame
dat = readRDS("data/phosphorylation/frequently-quantified-sites.rds") %>% 
  mutate(residue = strsplit(gene, '-') %>% 
           map_chr(~ extract(., length(.) - 1)),
         yval = paste0(gene_name, " p", residue, gsub("^.*-", "", gene)),
         xval = paste0(Accession, '-', Replicate, '-', fraction))

# subset to top-50 sites, not top-100, actually
top_sites = readRDS("data/QC/phosphosites-localization-scores.rds") %>% 
  filter(n_fractions > 0) %>% 
  filter(species == 'Homo sapiens') %>% 
  group_by(gene) %>% 
  summarise(fractions = sum(n_fractions)) %>% 
  ungroup() %>% 
  arrange(desc(fractions)) %>% 
  head(50) %>% 
  pull(gene)
dat %<>% filter(gene %in% top_sites)

# convert to matrix and cluster
mat = dat %>% 
  mutate(intensity = log(intensity)) %>% 
  dplyr::select(xval, yval, intensity) %>% 
  spread(yval, intensity) %>% 
  column_to_rownames('xval') %>% 
  as.matrix() %>% 
  # remove fractions never quantified
  extract(, colSums(is.finite(.) & . > 0) >= 1) %>% 
  # replace NAs with zeroes
  replace(!is.finite(.), 0)
dist_x = t(mat) %>% dist
dist_y = mat %>% dist
clust_x = hclust(dist_x)
clust_y = hclust(dist_y)
lvls_x = with(clust_x, labels[order])
lvls_y = with(clust_y, labels[order])
df = mat %>% 
  reshape2::melt(varnames = c('yval', 'xval'), value.name = 'intensity', 
                 as.is = TRUE) %>% 
  mutate(xval = factor(xval, levels = lvls_x), 
         yval = factor(yval, levels = lvls_y))
# limits = c(min, max)
limits = range(mat[mat > 0])
brks = c(limits[1] + 0.1 * diff(limits),
         limits[2] - 0.1 * diff(limits))
labels = format(limits, format = 'f', digits = 1)

# plot
p = df %>% 
  mutate(intensity = ifelse(intensity == 0, NA, intensity)) %>% 
  ggplot(aes(y = xval, x = yval, fill = intensity)) +
  geom_tile() +
  scale_x_discrete('Fraction', expand = c(0, 0)) +
  scale_y_discrete('Phosphosite', expand = c(0, 0)) +
  scale_fill_paletteer_c("pals::coolwarm",
                         name = 'Intensity, log10',
                         limits = limits, breaks = brks,
                         labels = labels,
                         na.value = 'grey96') +
  guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
  boxed_theme() +
  theme(legend.key.width = unit(0.18, 'lines'),
        legend.key.height = unit(0.18, 'lines'),
        aspect.ratio = 1.2,
        legend.position = 'right',
        legend.justification = 'bottom',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank())
p
ggsave("fig/phosphorylation/frequently-quantified-human-sites.pdf", p,
       width = 12, height = 10, units = "cm", useDingbats = FALSE)
