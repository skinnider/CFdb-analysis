setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(seqinr)
source("R/theme.R")

# read gene CDF
cdf = readRDS("data/QC/genes-CDF.rds") %>% 
  filter(species %in% c('Homo sapiens', 'Mus musculus', 'Arabidopsis thaliana'))

# bootstrap saturation curve
n_iter = 10
species = unique(cdf$species)
grid = tidyr::crossing(species = species, iteration = seq_len(n_iter))
saturation = pmap_dfr(grid, function(...) {
  current = tibble(...)
  print(current)
  
  # order experiments
  set.seed(current$iteration)
  expts0 = cdf %>% 
    filter(species == current$species) %>% 
    distinct(accession, experiment) %>% 
    sample_frac(1)
  # count genes
  map_dfr(seq_len(nrow(expts0)), ~ {
    genes = expts0 %>% 
      head(.x) %>% 
      left_join(cdf, by = c("accession", "experiment")) %>% 
      pull(gene) %>% 
      unique()
    data.frame(n_expts = .x, n_genes = n_distinct(genes))
  }) %>% 
    cbind(current, .)
})

# plot saturation curve
curve = saturation %>% 
  group_by(species, n_expts) %>% 
  summarise(mean = mean(n_genes), median = median(n_genes),
            sd = sd(n_genes), q1 = quantile(n_genes, probs = 0.25),
            q3 = quantile(n_genes, probs = 0.75)) %>% 
  ungroup()
blue = grafify::graf_palettes$bright[2] %>% unname
p = curve %>% 
  mutate(species = fct_relevel(species, 'Homo sapiens', 'Mus musculus',
                               'Arabidopsis thaliana')) %>% 
  ggplot(aes(x = n_expts, y = mean)) +
  facet_grid(~ species, scales = 'free_x', space = 'free_x') +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2,
              fill = 'blue') +
  geom_path(size = 0.35) + 
  scale_y_continuous(expression('# of proteins'~(10^3)), labels = ~ . / 1e3,
                     breaks = seq(0, 16, 2) * 1e3, limits = c(0, NA)) +
  scale_x_continuous('# of experiments', breaks = seq(0, 200, 25)) +
  boxed_theme()
p
ggsave("fig/QC/saturation-curve.pdf", p, width = 6.5, height = 4.5, 
       units = "cm", useDingbats = FALSE)

# fit saturation curve to human data
hs = filter(curve, species == 'Homo sapiens')
fit = nls( mean ~ alpha + beta * log(n_expts),
           start = list(alpha = 100, beta = 100),
           data = hs)
curve = data.frame(x = seq_len(max(hs$n_expts) * 2)) %>% 
  mutate(y = predict(fit, data.frame(n_expts = x)))
p2 = curve %>% 
  ggplot(aes(x = x, y = y)) +
  geom_ribbon(data = hs, aes(x = n_expts, y = mean, 
                             ymin = mean - sd, ymax = mean + sd), 
              alpha = 0.33, fill = 'blue') +
  geom_path(linetype = 'dotted', size = 0.3) + 
  scale_y_continuous(expression('# of proteins'~(10^3)), labels = ~ . / 1e3,
                     breaks = seq(0, 16, 2) * 1e3, limits = c(NA, 14e3)) +
  scale_x_continuous('# of experiments', breaks = seq(0, 400, 100)) +
  boxed_theme() +
  theme(aspect.ratio = 1)
p2
ggsave("fig/QC/saturation-curve-human-projected.pdf", p2, 
       width = 6.5, height = 4, units = "cm", useDingbats = FALSE)
