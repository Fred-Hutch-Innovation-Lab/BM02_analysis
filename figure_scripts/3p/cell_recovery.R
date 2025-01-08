# Setup ----
source(here('figure_scripts/utils.R'))
library(ggh4x)
library(readxl)
library(data.table)
library(ggbreak)

# Load data ----
fig_objs <- readRDS(here('rds/3p/02-objs_post_cell_filtering.rds'))
cell_loading_data <- read.csv(here('data/3p/loaded_cells.csv')) |>
  as.data.table() |>
  melt(id='Sample')

# Prepare plotdata ----
doublets <- lapply(fig_objs, function(x) {
  x@misc$filtering_receipts$after_doublet_filtering[[2]]
}) %>% stack()
colnames(doublets) <- c('value', 'Sample')
doublets$variable <- 'after_doublet_filtering'

bad_cells <- lapply(fig_objs, function(x) {
  x@misc$filtering_receipts$after_qc_filtering[[2]]
}) %>% stack() 
colnames(bad_cells) <- c('value', 'Sample')
bad_cells$variable <- 'after_qc_filtering'

target_cells <- filter(cell_loading_data, variable == 'target_cells') |>
  merge(metadata_3p, by = 'Sample')
target_cells_fraction <- filter(cell_loading_data, variable == 'target_cells_fraction') |>
  merge(metadata_3p, by = 'Sample')
plotdata <- do.call(rbind, list(cell_loading_data, doublets, bad_cells)) %>%
  filter(!variable %in% c('target_cells', 'target_cells_fraction')) |>
  filter(!is.na(value)) |>
  merge(metadata_3p, by='Sample') %>%
  mutate(Kit = factor(Kit, levels = kit_order_3p)) %>%
  arrange(Sample, variable) %>%
  mutate(variable = factor(variable, levels = c('cells_barcoded', 'cells_loaded', 'cells_recovered', 'after_qc_filtering', 'after_doublet_filtering'))) %>%
  group_by(Sample) %>%
  arrange(Sample, variable) %>%
  mutate(value=as.numeric(value),
         next_remaining = lead(value, default = 0),  # Next step's remaining count (or 0 for the last step)
         difference = value - next_remaining         # Difference from the next entry
  )
# Plot ----
ggplot(plotdata, aes(x=paste0(Individual, Replicate), y=difference, fill=variable)) +
  geom_col(position='stack') +
  geom_hline(data = target_cells, aes(yintercept = value), lty='dashed', color='black') +
  facet_wrap(~Kit, nrow=1, scales='free_x', labeller = labeller(Kit = label_function)) +
  theme(axis.text.x = element_text(angle=45, vjust=0.5)) + 
  ggbreak::scale_y_break(c(130000, 200000), ticklabels = c(200000, 230000))+ 
  scale_fill_manual(values = c('darkgrey', '#D7191C', '#FDAE61', '#e9e29c', '#39B185'),
                    labels = c('Remainder barcoded', 'Unrecovered cells', 'Low quality cells', 'Multiplets', 'High quality singlet')) +
  labs(x='Sample', y='Number of cells', fill='Classification', caption = 'Dashed lines indicate targeted cell recovery') ->
  figures[['cell_recovery_counts_break']] 
my_plot_save(image = figures[['cell_recovery_counts_break']] , 
             path = here('figures/3p/cell_recovery/cell_recovery_counts_break.svg'), 
             width = 7, height = 12)

ggplot(plotdata, aes(x=paste0(Individual, Replicate), y=difference, fill=variable)) +
  geom_col(position='stack') +
  geom_hline(data = target_cells, aes(yintercept = value), lty='dashed', color='black') +
  facet_wrap(~Kit, scales = 'free_x', nrow=1, labeller = labeller(Kit = label_function)) +
  theme(axis.text.x = element_text(angle=45, vjust=0.5)) + 
  scale_fill_manual(values = c('darkgrey', '#D7191C', '#FDAE61', '#e9e29c', '#39B185'),
                    labels = c('Remainder barcoded', 'Unrecovered cells', 'Low quality cells', 'Multiplets', 'High quality singlet')) +
  labs(x='Sample', y='Number of cells', fill='Classification') ->
  figures[['cell_recovery_counts']]
my_plot_save(image = figures[['cell_recovery_counts']] , 
             path = here('figures/3p/cell_recovery/cell_recovery_counts.svg'), 
             width = 7, height = 12)

ggplot(plotdata, aes(x=paste0(Individual, Replicate), y=difference, fill=variable)) +
  geom_col(position='stack') +
  geom_hline(data = target_cells, aes(yintercept = value), lty='dashed', color='black') +
  facet_wrap(~Kit, scales = 'free', nrow=1, labeller = labeller(Kit = label_function)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust=0.5)) + 
  scale_fill_manual(values = c('darkgrey', '#D7191C', '#FDAE61', '#e9e29c', '#39B185'),
                    labels = c('Remainder barcoded', 'Unrecovered cells', 'Low quality cells', 'Multiplets', 'High quality singlet')) +
  labs(x='Sample', y='Number of cells', fill='Classification') ->
  figures[['cell_recovery_counts_freey']] 
my_plot_save(image = figures[['cell_recovery_counts_freey']] , 
             path = here('figures/3p/cell_recovery/cell_recovery_counts_freey.svg'), 
             width = 7, height = 16)

plotdata %>% 
  group_by(Sample) %>%
  filter(variable != 'cells_barcoded') %>%
  mutate(difference = 100 * difference / sum (difference)) %>%
  ggplot(aes(x=paste0(Individual, Replicate), y=difference, fill=variable)) +
  geom_col(position='stack') +
  theme(axis.text.x = element_text(angle=45, vjust=0.5)) + 
  geom_hline(data = target_cells_fraction, aes(yintercept = 100*value), lty='dashed', color='black') +
  facet_wrap(~Kit, scales = 'free_x', nrow=1, labeller = labeller(Kit = label_function)) +
  scale_fill_manual(values = c('#D7191C', '#FDAE61', '#e9e29c', '#39B185'),
                    labels = c('Unrecovered cells', 'Low quality cells', 'Multiplets', 'High quality singlet')) +
  labs(x='Sample', y='% of capture', fill='Classification', caption = 'Dashed lines indicate expected cell recovery') ->
  figures[['cell_recovery_portions']] 
my_plot_save(image = figures[['cell_recovery_portions']] , 
             path = here('figures/3p/cell_recovery/cell_recovery_portions.svg'), 
             width = 7, height = 12)
