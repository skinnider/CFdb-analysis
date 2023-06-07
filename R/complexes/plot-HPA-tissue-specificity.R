# Plot the number of tissues in which detected in v1 vs. v2 are expressed. 
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read HPA
hpa = read_tsv("data/resources/HPA/normal_tissue.tsv.gz")

# calculate tissue-specificity
ts = hpa %>% 
  filter(Reliability != 'Uncertain', Level %in% c('High', 'Low', 'Medium')) %>% 
  distinct(`Gene name`, Tissue) %>% 
  dplyr::count(`Gene name`) %>% 
  dplyr::rename(gene = `Gene name`)

# read gene CDF
cdf = readRDS("data/QC/genes-CDF.rds")

# extract 2021 vs CFdb genes
genes_version = cdf %>% 
  filter(species == 'Homo sapiens') %>% 
  group_by(gene) %>% 
  arrange(version) %>% 
  dplyr::slice(1) %>% 
  ungroup()
# merge in abundance
ts %<>%
  left_join(genes_version, 'gene') %>% 
  replace_na(list(version = 'Never detected'))

# plot
version_pal = c('grey75', pals::stepped()[3]) %>%
  c('grey94')
p1 = ts %>% 
  mutate(version = fct_recode(version, 
                              'Version 1' = 'V1', 
                              'Version 2 only' = 'V2') %>% 
           fct_relevel('Version 1', 'Version 2 only', 'Undetected'),
         species = 'Homo sapiens') %>% 
  ggplot(aes(x = version, y = n, color = version, fill = version)) +
  facet_grid(~ species) +
  geom_violin(alpha = 0.3, color = NA, width = 0.65) +
  geom_boxplot(outlier.shape = NA, size = 0.3, width = 0.65, alpha = 0.4) + 
  scale_y_continuous('# of tissues expressed') +
  scale_color_manual('Detected in:', values = version_pal %>% 
                       replace(. == 'grey94', 'grey65')) +
  scale_fill_manual('Detected in:', values = version_pal) +
  boxed_theme() +
  theme(aspect.ratio = 1.4,
        legend.position = 'none',
        legend.key.width = unit(0.35, 'lines'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())
p1
ggsave("fig/abundance/HPA-tissue-specificity.pdf", p1,
       width = 8, height = 4.85, units = "cm", useDingbats = FALSE)
