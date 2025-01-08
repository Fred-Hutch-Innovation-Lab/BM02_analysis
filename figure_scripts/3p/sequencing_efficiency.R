# Setup ----
source(here('figure_scripts/utils.R'))
library(readxl)

# Load data ----
data <- read.csv(here('data/3p/sequencing_efficiency.csv'), skip = 1) 

# Prepare plotdata ----

expected_data <- data |>
  filter(library == 'all') |>
  select(c(kit, expected_reads)) |>
  filter(kit != 'Parse_v2') |>
  mutate(kit = factor(kit, levels = kit_order_3p))
  
plotdata <- data |>
  filter(library == 'all') |>
  select(-c(library, percent_passing_filter_cluster)) |>
  pivot_longer(cols=c(expected_reads,
                      passing_filter_reads,
                      reads_in_fastqs,
                      reads_in_cells_and_mapped),
               names_to = 'stage', values_to='reads') |>
  mutate(stage = factor(stage, 
                        levels = c('expected_reads', 'passing_filter_reads', 
                                   'reads_in_fastqs', 'reads_in_cells_and_mapped')),
         kit = factor(kit, levels=kit_order_3p)) |>
    filter(stage != 'expected_reads') |>
    group_by(kit) |> 
    arrange(kit, stage) |>
    mutate(difference = reads - lead(reads, default = 0)) |>
  na.omit()

# Plot ----
plotdata |>
  ggplot(aes(x=kit, y=difference, group=stage, fill=stage)) +
  geom_bar(stat='identity', position='stack') +
  facet_wrap(~ kit, scales='free_x', nrow=1) +
  geom_hline(data=expected_data, aes(yintercept = expected_reads, group = kit), lty='dashed') +
  scale_fill_manual(values = c('#D7191C', '#e9e29c', '#39B185'),
                     labels = c('Not extracted to fastq',
                                'Not barcoded and aligned\nin cells', 
                                'Successfully aligned and\nassigned to cell')) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.spacing.x=unit(0, "lines")) + 
  labs(x='Kit', y='Reads', fill='Read fate', caption = 'Line indicates expected read recovery based on sequencer loading') ->
  figures[['seq_eff_count']]

my_plot_save(image = figures[['seq_eff_count']], 
             path = here('figures/3p/sequencing_efficiency/seq_eff_count.svg'), 
             width = 10, height = 5)

plotdata |>
  mutate(prop = difference / sum(difference)) |>
  ggplot(aes(x=kit, y=prop, group=stage, fill=stage)) +
  geom_bar(stat='identity', position='stack') +
  scale_fill_manual(values = c('#D7191C', '#e9e29c', '#39B185'),
                    labels = c('Not extracted to fastq',
                               'Not barcoded and\naligned in cells', 
                               'Successfully aligned and\nassigned to cell')) +
  scale_y_continuous(labels = scales::percent) +
  labs(x='Kit', y='% of reads', fill='Read fate') ->
  figures[['seq_eff_prop']]

my_plot_save(image = figures[['seq_eff_count']], 
             path = here('figures/3p/sequencing_efficiency/seq_eff_prop.svg'), 
             width = 10, height = 5)
