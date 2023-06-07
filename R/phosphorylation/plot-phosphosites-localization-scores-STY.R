setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read localization scores
dat = readRDS("data/QC/phosphosites-localization-scores.rds") %>% 
  filter(n_fractions > 0) %>% 
  # flag STY
  mutate(sty = strsplit(gene, '-') %>% map_chr(~ extract(., length(.) - 1))) %>% 
  # cut fractions at a max of 10
  # filter(n_fractions <= 10)
  mutate(n_fractions = winsorize(n_fractions, c(NA, 10)))

# plot by residue
pal = pals::stepped3()[c(1, 5) + 2] %>% c(pals::kelly()[1])
p = dat %>% 
  ggplot(aes(x = sty, y = localization_prob, fill = sty)) +
  geom_blank() +
  geom_rect(data = head(dat, 1), 
            aes(xmin = -Inf, xmax = Inf, ymin = 0.75, ymax = Inf),
            fill = 'grey97') +
  geom_hline(aes(yintercept = 0.75), linetype = 'dotted', size = 0.25) +
  # geom_violin(color = NA, fill = 'grey75', width = 1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.7, size = 0.3,
               color = 'grey30') + 
  scale_x_discrete('# of fractions', labels = ~ replace(., . == '10', '10+')) +
  scale_y_continuous('Localization probability') +
  scale_fill_manual('', values = pal) +
  boxed_theme() +
  theme(aspect.ratio = 2,
        legend.key.width = unit(0.35, 'lines'))
p
ggsave("fig/phosphorylation/phosphosite-localization-scores-STY.pdf", p,
       width = 5, height = 5, units = "cm", useDingbats = FALSE)
