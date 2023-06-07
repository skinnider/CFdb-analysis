setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read localization scores
dat = readRDS("data/QC/phosphosites-localization-scores.rds") %>% 
  filter(n_fractions > 0) %>% 
  # cut fractions at a max of 5
  mutate(n_fractions = winsorize(n_fractions, c(NA, 5)))

## number of fractions vs. localization ####
p1 = dat %>% 
  ggplot(aes(x = factor(n_fractions), y = localization_prob)) +
  geom_blank() +
  geom_rect(data = head(dat, 1), 
            aes(xmin = -Inf, xmax = Inf, ymin = 0.75, ymax = Inf),
            fill = 'grey97') +
  geom_hline(aes(yintercept = 0.75), linetype = 'dotted', size = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.6, size = 0.3,
               fill = pals::kelly()[1], color = 'grey30') + 
  scale_x_discrete('# of fractions', labels = ~ replace(., . == '5', '5+')) +
  scale_y_continuous('Localization probability') +
  boxed_theme() +
  theme(aspect.ratio = 2)
p1

## number of fractions vs. delta score ####
p2 = dat %>% 
  ggplot(aes(x = factor(n_fractions), y = delta_score)) +
  geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.6, size = 0.3,
               fill = pals::kelly()[1], color = 'grey30') + 
  scale_x_discrete('# of fractions', labels = ~ replace(., . == '5', '5+')) +
  scale_y_continuous('Delta score') +
  coord_cartesian(ylim = c(0, 250)) +
  boxed_theme() +
  theme(aspect.ratio = 2)
p2

## number of fractions vs. intensity ####
p3 = dat %>% 
  ggplot(aes(x = factor(n_fractions), y = intensity)) +
  geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.6, size = 0.3,
               fill = pals::kelly()[1], color = 'grey30') + 
  scale_x_discrete('# of fractions', labels = ~ replace(., . == '5', '5+')) +
  scale_y_log10('Intensity', labels = fancy_scientific) +
  coord_cartesian(ylim = c(1e4, 1e11)) +
  annotation_logticks(sides = 'l', short = unit(0.06, 'lines'),
                      mid = unit(0.12, 'lines'), long = unit(0.18, 'lines'),
                      size = 0.2) +
  boxed_theme() +
  theme(aspect.ratio = 2)
p3

## number of fractions vs. stoichiometry ####
p4 = dat %>% 
  drop_na(ratio_mod_base) %>% 
  ggplot(aes(x = factor(n_fractions), y = ratio_mod_base)) +
  geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.6, size = 0.3,
               fill = pals::kelly()[1], color = 'grey30') + 
  scale_x_discrete('# of fractions', labels = ~ replace(., . == '5', '5+')) +
  scale_y_log10('Stoichiometry', labels = fancy_scientific) +
  annotation_logticks(sides = 'l', short = unit(0.06, 'lines'),
                      mid = unit(0.12, 'lines'), long = unit(0.18, 'lines'),
                      size = 0.2) + 
  coord_cartesian(ylim = c(1e-5, 1e3)) +
  boxed_theme() +
  theme(aspect.ratio = 2)
p4

# combine and save
p = p1 | p2 | p3 | p4
p
ggsave("fig/phosphorylation/phosphosite-localization-scores.pdf", p,
       width = 18, height = 4, units = "cm", useDingbats = FALSE)

################################################################################
### repeat all plots by experiment 
################################################################################

# read localization scores
dat = readRDS("data/QC/phosphosites-localization-scores.rds") %>% 
  filter(n_fractions > 0) %>% 
  group_by(gene) %>% 
  mutate(n_expts = n()) %>% 
  ungroup() %>% 
  mutate(n_expts = winsorize(n_expts, c(NA, 5)))

## number of experiments vs. localization ####
p5 = dat %>% 
  ggplot(aes(x = factor(n_expts), y = localization_prob)) +
  geom_blank() +
  geom_rect(data = head(dat, 1), 
            aes(xmin = -Inf, xmax = Inf, ymin = 0.75, ymax = Inf),
            fill = 'grey97') +
  geom_hline(aes(yintercept = 0.75), linetype = 'dotted', size = 0.25) +
  # geom_violin(color = NA, fill = 'grey75', width = 1) +
  geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.6, size = 0.3,
               fill = pals::kelly()[1], color = 'grey30') + 
  scale_x_discrete('# of experiments', labels = ~ replace(., . == '5', '5+')) +
  scale_y_continuous('Localization probability') +
  boxed_theme() +
  theme(aspect.ratio = 2)
p5

## number of experiments vs. delta score ####
p6 = dat %>% 
  ggplot(aes(x = factor(n_expts), y = delta_score)) +
  geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.6, size = 0.3,
               fill = pals::kelly()[1], color = 'grey30') + 
  scale_x_discrete('# of experiments', labels = ~ replace(., . == '5', '5+')) +
  scale_y_continuous('Delta score') +
  coord_cartesian(ylim = c(0, 250)) +
  boxed_theme() +
  theme(aspect.ratio = 2)
p6

## number of experiments vs. intensity ####
p7 = dat %>% 
  ggplot(aes(x = factor(n_expts), y = intensity)) +
  geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.6, size = 0.3,
               fill = pals::kelly()[1], color = 'grey30') + 
  scale_x_discrete('# of experiments', labels = ~ replace(., . == '5', '5+')) +
  scale_y_log10('Intensity', labels = fancy_scientific) +
  coord_cartesian(ylim = c(1e4, 1e11)) +
  annotation_logticks(sides = 'l', short = unit(0.06, 'lines'),
                      mid = unit(0.12, 'lines'), long = unit(0.18, 'lines'),
                      size = 0.2) +
  boxed_theme() +
  theme(aspect.ratio = 2)
p7

## number of experiments vs. stoichiometry ####
p8 = dat %>% 
  drop_na(ratio_mod_base) %>% 
  ggplot(aes(x = factor(n_expts), y = ratio_mod_base)) +
  # geom_hline(aes(yintercept = 1), size = 0.44, color = 'grey80') +
  geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.6, size = 0.3,
               fill = pals::kelly()[1], color = 'grey30') + 
  scale_x_discrete('# of experiments', labels = ~ replace(., . == '5', '5+')) +
  scale_y_log10('Stoichiometry', labels = fancy_scientific) +
  annotation_logticks(sides = 'l', short = unit(0.06, 'lines'),
                      mid = unit(0.12, 'lines'), long = unit(0.18, 'lines'),
                      size = 0.2) + 
  coord_cartesian(ylim = c(1e-5, 1e3)) +
  boxed_theme() +
  theme(aspect.ratio = 2)
p8

# combine and save
p = p5 | p6 | p7 | p8
p
ggsave("fig/phosphorylation/phosphosite-localization-scores-experiments.pdf", p,
       width = 18, height = 4, units = "cm", useDingbats = FALSE)
