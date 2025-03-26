# Setup ----
source(here('figure_scripts/utils.R'))

# Load data ----
objs <- readRDS(here('rds/3p/pseudobulk_objs_split_by_kit_allgenes_downsamplecells_coarse_celltypes.Rds'))
varPart <- readRDS(here('rds/3p/varpart_results.Rds'))

# Generate plotdata ----
plotdata <- lapply(varPart, function(x) {
  as.data.frame(x) |> 
    as.data.table() |> 
    melt()
  }) |>
  rbindlist(idcol = 'Kit') |>
  mutate(Kit = factor(Kit, kit_order_3p), 
         variable = factor(variable, c('celltype', 'Individual', 'Replicate', 'Residuals')))
colnames(plotdata) <- c('Kit', 'Category', 'Value')

# Plot ----
ggplot(plotdata, aes(x=Kit, y=Value, fill = Kit))  +
  geom_boxplot(outliers = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(labels = label_function()) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values = color_palette$kits, labels = label_function()) + 
  facet_wrap(~ Category, scales='free_y', nrow=1) +
  labs(x='Kit', y='Percent of variance in expression explained') ->
  figures[['gene_variability']]

# Save plot
my_plot_save(figures[['gene_variability']], 
             path = here('figures/3p/variance_partitioning_boxplot.svg'),
             width = 18,
             height = 6)

# Save plot data
write_plot_data(plotdata,
                file = here('figure_data/3p/variance_partitioning_boxplot.txt'))

