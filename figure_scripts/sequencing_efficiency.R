library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)    
library(readxl)

kit_order <- read.table(here('config/kit_order.txt'))$V1

data <- read.csv(here('data/sequencing_efficiency.csv'), skip = 1) 

plotdata1 <- data |>
  filter(library == 'all') |>
  select(c(kit, expected_reads)) |>
  filter(kit != 'Parse_v2') |>
  mutate(kit = factor(kit, levels = kit_order))
  
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
         kit = factor(kit, levels=kit_order)) |>
    filter(stage != 'expected_reads') |>
    group_by(kit) |> 
    arrange(kit, stage) |>
    mutate(difference = reads - lead(reads, default = 0)) |>
  na.omit()

plotdata |>
  ggplot(aes(x=kit, y=difference, group=stage, fill=stage)) +
  geom_bar(stat='identity', position='stack') +
  facet_wrap(~ kit, scales='free_x', nrow=1) +
  geom_hline(data=plotdata1, aes(yintercept = expected_reads, group = kit), lty='dashed') +
  scale_fill_manual(values = c('#D7191C', '#e9e29c', '#39B185'),
                     labels = c('Not extracted to fastq',
                                'Not barcoded and aligned\nin cells', 
                                'Successfully aligned and\nassigned to cell')) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.margin.x=unit(0, "lines")) + 
  labs(x='Kit', y='Reads', fill='Read fate', caption = 'Line indicates expected read recovery based on sequencer loading') 
ggsave(here('figures/sequencing_efficiency/seq_eff_count.png'), device = 'png', 
       width=unit(6.5, 'in'), height = unit(5, 'in'))

plotdata |>
  mutate(prop = difference / sum(difference)) |>
  ggplot(aes(x=kit, y=prop, group=stage, fill=stage)) +
  geom_bar(stat='identity', position='stack') +
  scale_fill_manual(values = c('#D7191C', '#e9e29c', '#39B185'),
                    labels = c('Not extracted to fastq',
                               'Not barcoded and aligned in cells', 
                               'Successfully aligned and assigned to cell')) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent) +
  labs(x='Kit', y='% of reads', fill='Read fate') 
ggsave(here('figures/sequencing_efficiency/seq_eff_prop.png'), device = 'png', 
       width=unit(6.5, 'in'), height = unit(5, 'in'))
