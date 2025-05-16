# Setup ----
source(here('figure_scripts/utils.R'))
library(readxl)

# Load data ----
data <- read.csv(here('data/3p/sequencing_efficiency_3.txt'), sep='\t') 
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
  merge(metadata_3p, by='Sample') |>
  mutate(Kit = factor(Kit, levels=kit_order_3p)) 

# Plot ----
plotdata |>
  ggplot(aes(x=paste0(Individual, Replicate), y=prop, group = variable, fill=variable)) +
  geom_bar(stat='identity', position='stack') +
  scale_fill_manual(values = c('#D7191C', '#e9e29c', '#39B185'), 
                    labels = c('Not aligned or\nin called cells',
                               'Aligned but\nnot in called cells',
                               'Aligned and\nassigned to cell')) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  facet_wrap(~ Kit, scales='free_x', nrow=1, labeller = labeller(Kit = label_function())) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  labs(x='Sample', y='% of extracted reads', fill='Read utilization') ->
  figures[['seq_eff_prop']]
figures[['seq_eff_prop']]

my_plot_save(image = figures[['seq_eff_prop']], 
             path = here('figures/3p/sequencing_efficiency/seq_eff_prop.svg'), 
             width = 10.5, height = 4)
write_plot_data(plotdata, here('figure_data/3p/seq_eff_count.txt'))
