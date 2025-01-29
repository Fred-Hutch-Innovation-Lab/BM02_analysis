# Setup ----
source(here('figure_scripts/utils.R'))

sat_curves <- read.csv(here('data/5p/sequencing_saturation.csv')) |>
  select(-Kit) |>
  mutate(across(all_of(c('median_genes', 'median_umi', 'nreads')), as.numeric))
sat_curves <- merge(sat_curves, metadata_5p, by='Sample', all=TRUE)

ggplot(filter(sat_curves, !is.na(median_genes)),
       aes(x=nreads, y=median_genes, linetype = Individual, color = Kit)) +
  geom_line() +
  scale_color_manual(values = unlist(color_palette$kits), labels = label_function) +
  guides(color = guide_legend(order = 1), 
         linetype = guide_legend(order = 2)) +
  # geom_polygon(data = polygon_data, aes(x=nreads, y=median_genes, group=Kit, fill=Kit), alpha=0.1, inherit.aes = FALSE, show.legend=FALSE)+
  labs(x='Average reads per cell', y='Median genes per cell', linetype='Sample', color = 'Kit') ->
  figures[['5p_sat_curves_genes']]
my_plot_save(figures[['5p_sat_curves_genes']], 
             here('figures/5p/saturation_curves_genes.svg'), 
             width = 7, height = 5)
write_plot_data(sat_curves, here('figure_data/5p/saturation_curves_genes.txt'))

ggplot(filter(sat_curves, !is.na(median_umi)),
       aes(x=nreads, y=median_umi, linetype = Individual, color = Kit)) +
  geom_line() +
  scale_color_manual(values = unlist(color_palette$kits), labels = label_function) +
  guides(color = guide_legend(order = 1), 
         linetype = guide_legend(order = 2)) +
  # geom_polygon(data = polygon_data, aes(x=nreads, y=median_genes, group=Kit, fill=Kit), alpha=0.1, inherit.aes = FALSE, show.legend=FALSE)+
  labs(x='Average reads per cell', y='Median UMI per cell', linetype='Sample', color = 'Kit') ->
  figures[['5p_sat_curves_UMI']]
my_plot_save(figures[['5p_sat_curves_UMI']], 
             here('figures/5p/saturation_curves_UMI.svg'), 
             width = 7, height = 5)
write_plot_data(sat_curves, here('figure_data/5p/saturation_curves_UMI.txt'))
