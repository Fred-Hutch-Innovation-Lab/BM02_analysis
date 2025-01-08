# Setup ----
source(here('figure_scripts/utils.R'))

# Load data ----
comp_table <- read.table(here('rds/3p/composition_table_coarse.txt'), sep = '\t', header = TRUE)

# Plotting ----
comp_table |> 
  reshape2::melt(id.vars = c('Kit', 'Sample', 'Individual', 'Replicate')) |> #measure.vars=unique(coarse_mapping)) |> 
  mutate(variable = factor(variable, levels = c(
    'T', 'Monocyte', 'B', 'NK', 'Dendritic', 'Megakaryocyte', 'pDC', 
    'Granulocyte', 'Erythrocyte', 'Unknown'
  ))) |>
  group_by(Sample) |> 
  mutate(percent = 100*value / sum(value),
         Kit = factor(Kit, levels=kit_order_3p)) |>
  mutate(percent = ifelse(percent == 0, NA, percent)) |>
  ggplot(aes(x = Kit, y = percent,
             group = Kit, color = Kit, shape = paste0(Individual, Replicate))) + 
  geom_point() +
  labs(x = element_blank(), y = '% of sample',
       color = 'Kit', shape='Sample'#,
       # caption = paste(
         # 'Portions are of the total capture for a sample.',
         # '* = significant below 0.05 FDR',
         # sep='\n')
) +
  scale_shape_manual(values = c(1:5))+ 
  scale_color_manual(values=color_palette$kits, labels = label_function) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_grid(~ factor(variable), scales='free_y') ->
  figures[['composition_boxplot']]

  my_plot_save(image = figures[['composition_boxplot']], 
               path = here('figures/3p/composition_coarse_labels.svg'), 
               width = 13, height = 4)
  