# Setup ----
source(here('figure_scripts/utils.R'))

# Load data ----
# data <- data |>
  # filter(kit != 'Parse_v2') |>
  # mutate(approx_total_reads = passing_filter_reads / (percent_passing_filter_cluster / 100))
data <- read.csv(here('data/5p/sequencing_efficiency.txt'), sep='\t') 
# Prepare plotdata ----
plotdata <- data |>
  as.data.table() |>
  melt() |>
  mutate(variable = factor(variable, levels = c('reads_in_fastqs', 'reads_mapped_transcriptome', 'reads_in_matrix'))) |>
  group_by(Sample) |>
  arrange(Sample, variable) |> 
  mutate(value=as.numeric(value),
         next_remaining = lead(value, default = 0),  # Next step's remaining count (or 0 for the last step)
         difference = value - next_remaining         # Difference from the next entry
  ) |> 
  mutate(prop = difference / sum(difference)) |>
  # mutate(prop = reads_in_matrix / reads_in_fastqs) |>
  # select(c(Sample, prop)) |>
  merge(metadata_5p, by='Sample') |>
  mutate(Kit = factor(Kit, levels=kit_order_5p)) 

plotdata  |>
  ggplot(aes(x=Individual, y=prop, group = variable, fill=variable)) +
  geom_bar(stat='identity', position='stack') +
  scale_fill_manual(values = c('#D7191C', '#e9e29c', '#39B185'), #'#D7191C',
                    labels = c('Not aligned or\nin called cells',
                               'Aligned but\nnot in called cells',
                               'Aligned and\nassigned to cell')) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  # scale_x_discrete(labels = label_function) +
  facet_wrap(~ Kit, scales='free_x', nrow=1, labeller = labeller(Kit = label_function(mode='clean'))) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  # theme(axis.text.x = element_blank(), 
  #       axis.ticks.x = element_blank(),
  #       panel.spacing.x=unit(0, "lines")) + 
  labs(x='Kit', y='% of extracted reads', fill='Read utilization') ->
  figures[['seq_eff_prop']]

my_plot_save(image = figures[['seq_eff_prop']], 
             path = here('figures/5p/sequencing_efficiency/seq_eff_prop.svg'), 
             width = 6.5, height = 4)
write_plot_data(plotdata, here('figure_data/5p/wt/seq_eff_count.txt'))
