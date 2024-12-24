library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations

source(here('config/color_palette.R'))
source(here('config/kit_order.R'))

# 3P ----
metadata <- read.csv(here('config/3p/metadata.csv')) %>%
  mutate(Kit = factor(Kit, levels = kit_order_3p))
sat_curves <- read.csv(here('data/3p/saturation_curves/manually_downsampled.csv')) |>
  select(-Kit) |>
  mutate(across(all_of(c('median_genes', 'median_umi', 'nreads')), as.numeric))
zeroes <- data.frame(Sample=metadata$Sample, nreads=0, median_genes=0, median_umi=0)
sat_curves <- rbind(sat_curves, zeroes)
sat_curves <- merge(sat_curves, metadata, by='Sample', all=TRUE)

samples <- c('Scale_F5A', 'Scale_F1A',
             'NextGEM3P_F1A', 'NextGEM3P_F1B',
             'GEMX3P_F1A', 'GEMX3P_F1B',
             'Parse_v3_F5A', 'Parse_v3_F1B',
             'Flex_F5A', 'Flex_F1B',
             'Fluent_V_F1A', 'Fluent_V_F5B',
             'Fluent_v4_F5A', 'Fluent_v4_F5B')
polygon_data <- sat_curves %>%
  filter(Sample %in% samples) %>%
  mutate(Sample = factor(Sample, levels=samples)) %>%
  mutate(sort_key = case_when(
    Sample %in% c('Scale_F5A', 'NextGEM3P_F1A', 'GEMX3P_F1A', 'Parse_v3_F5A', 'Flex_F5A', 'Fluent_V_F1A', 'Fluent_v4_F5A') ~ nreads,
    Sample %in% c('Scale_F1A', 'NextGEM3P_F1B', 'GEMX3P_F1B', 'Parse_v3_F1B', 'Flex_F1B', 'Fluent_V_F5B', 'Fluent_v4_F5B') ~ -nreads, 
    TRUE ~ NA
  )) %>%
  arrange(Kit, Sample, sort_key) |>
  filter(!is.na(median_genes))

ggplot(filter(sat_curves, !is.na(median_genes)),
       aes(x=nreads, y=median_genes, linetype = paste0(Individual, Replicate), color = Kit)) +
  geom_line() +
  scale_color_manual(values = unlist(color_palette$kits), labels = label_function) +
  geom_polygon(data = polygon_data, aes(x=nreads, y=median_genes, group=Kit, fill=Kit), alpha=0.1, inherit.aes = FALSE, show.legend=FALSE)+
  theme_bw() +
  labs(x='Average reads per cell', y='Median genes per cell', linetype='Sample', color = 'Kit')
ggsave('saturation_curves_genes.png', path = here('figures/3p'), width = unit(7, 'in'), height = unit(5, 'in'))

samples <- c('Scale_F5A', 'Scale_F1A',
             'NextGEM3P_F1A', 'NextGEM3P_F5B',
             'GEMX3P_F1A', 'GEMX3P_F1B',
             'Parse_v3_F5A', 'Parse_v3_F1B',
             'Flex_F5A', 'Flex_F1B',
             'Fluent_V_F1A', 'Fluent_V_F5B',
             'Fluent_v4_F5A', 'Fluent_v4_F5B')
polygon_data <- sat_curves %>%
  filter(Sample %in% samples) %>%
  mutate(Sample = factor(Sample, levels=samples)) %>%
  mutate(sort_key = case_when(
    Sample %in% c('Scale_F5A', 'NextGEM3P_F1A', 'GEMX3P_F1A', 'Parse_v3_F5A', 'Flex_F5A', 'Fluent_V_F1A', 'Fluent_v4_F5A') ~ nreads,
    Sample %in% c('Scale_F1A', 'NextGEM3P_F5B', 'GEMX3P_F1B', 'Parse_v3_F1B', 'Flex_F1B', 'Fluent_V_F5B', 'Fluent_v4_F5B') ~ -nreads, 
    TRUE ~ NA
  )) %>%
  arrange(Kit, Sample, sort_key) |>
  filter(!is.na(median_genes))

ggplot(filter(sat_curves, !is.na(median_umi)),
       aes(x=nreads, y=median_umi, linetype = paste0(Individual, Replicate), color = Kit)) +
  geom_line() +
  scale_color_manual(values = unlist(color_palette$kits), labels = label_function, na.translate = FALSE) +
  geom_polygon(data = polygon_data, aes(x=nreads, y=median_umi, group=Kit, fill=Kit), alpha=0.1, inherit.aes = FALSE, show.legend=FALSE)+
  theme_bw() +
  labs(x='Average reads per cell', y='Median transcripts per cell', linetype='Sample', color = 'Kit', caption = 'Flex data not generated')
ggsave('saturation_curves_umi.png', path = here('figures/3p'), width = unit(7, 'in'), height = unit(5, 'in'))
