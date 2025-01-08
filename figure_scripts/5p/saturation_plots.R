library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations

source(here('config/color_palette.R'))
source(here('config/kit_order.R'))
metadata <- read.csv(here('config/5p/metadata.csv')) %>%
  mutate(Kit = factor(Kit, levels = kit_order_5p))

sat_curves <- read.csv(here('data/5p/sequencing_saturation.csv')) |>
  select(-Kit) |>
  mutate(across(all_of(c('median_genes', 'median_umi', 'nreads')), as.numeric))
sat_curves <- merge(sat_curves, metadata, by='Sample', all=TRUE)

ggplot(filter(sat_curves, !is.na(median_genes) & nreads <= 25000),
       aes(x=nreads, y=median_genes, linetype = Individual, color = Kit)) +
  geom_line() +
  scale_color_manual(values = unlist(color_palette$kits)) +
  # geom_polygon(data = polygon_data, aes(x=nreads, y=median_genes, group=Kit, fill=Kit), alpha=0.1, inherit.aes = FALSE, show.legend=FALSE)+
  theme_bw() +
  labs(x='Average reads per cell', y='Median genes per cell', linetype='Sample', color = 'Kit')
ggsave('saturation_curves_genes.png', path = here('figures/5p'), width = unit(7, 'in'), height = unit(5, 'in'))
