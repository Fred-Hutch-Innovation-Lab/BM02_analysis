# Setup ----
source(here('figure_scripts/utils.R'))
library(readxl)

# Load data ----
data <- read.csv(here('data/3p/sequencing_efficiency.csv'), skip = 1) 
data <- data |>
  filter(kit != 'Parse_v2') |>
  mutate(approx_total_reads = passing_filter_reads / (percent_passing_filter_cluster / 100))
data2 <- read.csv(here('data/3p/sequencing_efficiency_3.txt'), sep='\t') 
# Prepare plotdata ----

expected_data <- data |>
  filter(library == 'all') |>
  select(c(kit, expected_reads)) |>
  filter(kit != 'Parse_v2') |>
  mutate(kit = factor(kit, levels = kit_order_3p))
  
plotdata <- data |>
  filter(library == 'all') |>
  select(-c(library, percent_passing_filter_cluster)) |>
  pivot_longer(cols=c(approx_total_reads,
                      # expected_reads,
                      passing_filter_reads,
                      reads_in_fastqs,
                      reads_in_cells_and_mapped),
               names_to = 'stage', values_to='reads') |>
  mutate(stage = factor(stage, 
                        levels = c('approx_total_reads',
                                   # 'expected_reads',
                                   'passing_filter_reads',
                                   'reads_in_fastqs',
                                   'reads_in_cells_and_mapped')),
         kit = factor(kit, levels=kit_order_3p)) |>
    filter(stage != 'expected_reads') |>
    group_by(kit) |> 
    arrange(kit, stage) |>
    mutate(difference = reads - lead(reads, default = 0)) #|>
  # na.omit() |>
  # mutate(prop = difference / sum(difference))

# Plot ----
plotdata |>
  ggplot(aes(x=kit, y=difference, group=stage, fill=stage)) +
  geom_bar(stat='identity', position='stack') +
  facet_wrap(~ kit, scales='free_x', nrow=1) +
  geom_hline(data=expected_data, aes(yintercept = expected_reads, group = kit), lty='dashed') +
  scale_fill_manual(values = c('#D7191C', 'orange', '#e9e29c', '#39B185'), #'#D7191C',
                     labels = c('Reads not passing filter',
                       'Not extracted to fastq',
                                'Not barcoded and aligned\nin cells', 
                                'Successfully aligned and\nassigned to cell')) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.spacing.x=unit(0, "lines")) + 
  scale_x_discrete(labels = label_function) +
  labs(x='Sample', y='Reads', fill='Read fate', caption = 'Line indicates expected read recovery based on sequencer loading') ->
  figures[['seq_eff_count']]

my_plot_save(image = figures[['seq_eff_count']], 
             path = here('figures/3p/sequencing_efficiency/seq_eff_count.svg'), 
             width = 10, height = 5)
write_plot_data(plotdata, here('figure_data/3p/seq_eff_count.txt'))

plotdata <- data2 |>
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
  merge(metadata_3p, by='Sample') |>
  mutate(Kit = factor(Kit, levels=kit_order_3p)) 

plotdata  |>
  ggplot(aes(x=paste0(Individual, Replicate), y=prop, group = variable, fill=variable)) +
  geom_bar(stat='identity', position='stack') +
  scale_fill_manual(values = c('#D7191C', '#e9e29c', '#39B185'), #'#D7191C',
                    labels = c('Not aligned or\nin called cells',
                               'Aligned but\nnot in called cells',
                               'Aligned and\nassigned to cell')) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  # scale_x_discrete(labels = label_function) +
  facet_wrap(~ Kit, scales='free_x', nrow=1, labeller = labeller(Kit = label_function)) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  # theme(axis.text.x = element_blank(), 
  #       axis.ticks.x = element_blank(),
  #       panel.spacing.x=unit(0, "lines")) + 
  labs(x='Sample', y='% of extracted reads', fill='Read utilization') ->
  figures[['seq_eff_prop']]

my_plot_save(image = figures[['seq_eff_prop']], 
             path = here('figures/3p/sequencing_efficiency/seq_eff_prop.svg'), 
             width = 10.5, height = 4)
