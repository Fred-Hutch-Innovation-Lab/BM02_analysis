library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)       ## parses config yaml
library(Seurat)
library(readxl)
library(data.table)
library(ggbreak)

fig_objs <- readRDS(here('rds/3p/02-objs_post_cell_filtering.rds'))
source(here('config/kit_order.R'))
metadata <- read.csv(here('config/3p/metadata.csv')) %>%
  mutate(Kit = factor(Kit, levels = kit_order_3p))

main <- function() {
  
  cell_recovery_figures <- list()
  
  target_cells <- list(
    'Flex' = 10000,
    'NextGEM3P' = 10000,
    'GEMX3P' = 20000,
    'Fluent_v4' = 20000,
    'Fluent_V' = 20000,
    'Parse_v3' = 25000,
    'Scale' = 125000 /4
  ) |> reshape2::melt() |>
    dplyr::rename(Kit = L1) |>
    mutate(Kit = factor(Kit, levels=kit_order_3p))
  
  target_cells_fraction <- list(
    'Flex' = 10000/16871,
    'NextGEM3P' = 10000/(1100*15),
    'GEMX3P' = 20000/(1300*22.3),
    'Fluent_v4' = .5,
    'Fluent_V' = .5,
    'Parse_v3' = 100000/(4*(520*14*12)),
    'Scale' = .25
  ) |> reshape2::melt() |>
    dplyr::rename(Kit = L1)  |>
    mutate(Kit = factor(Kit, levels=kit_order_3p))
  
  pipeline_data <- read_xlsx(here('data/3p/pipeline_summary_statistics/downsampled_data.xlsx'), skip = 1, .name_repair = 'minimal', col_names = TRUE) %>%
    filter(METRICS %in% c('# cells loaded', '# Cells recovered')) %>%
    as.data.table() |>
    melt(id.vars='METRICS') |>
    mutate(value = as.numeric(value)) %>%
    mutate(METRICS = case_when(
      METRICS == '# cells loaded' & grepl('Scale', variable) ~ '# cells barcoded',
      .default = METRICS
    ))
  colnames(pipeline_data) <- c('class', 'Sample', 'Freq')
  scale_used_cells <- 125000
  pipeline_data <- rbind(pipeline_data,
    data.frame(class = '# cells loaded', 
               Sample = metadata$Sample[metadata$Kit=='Scale'], 
             Freq = scale_used_cells
  ))
    
  doublets <- lapply(fig_objs, function(x) {
    x@misc$filtering_receipts$after_doublet_filtering[[2]]
  }) %>% stack()
  colnames(doublets) <- c('Freq', 'Sample')
  doublets$class <- 'after_doublet_filtering'
  
  bad_cells <- lapply(fig_objs, function(x) {
    x@misc$filtering_receipts$after_qc_filtering[[2]]
  }) %>% stack() 
  colnames(bad_cells) <- c('Freq', 'Sample')
  bad_cells$class <- 'after_qc_filtering'
  
  
  plotdata <- do.call(rbind, list(pipeline_data, doublets, bad_cells)) %>%
    merge(metadata, by='Sample') %>%
    mutate(Kit = factor(Kit, levels = kit_order_3p)) %>%
    arrange(Sample, class) %>%
    mutate(class = factor(class, levels = c('# cells barcoded', '# cells loaded', '# Cells recovered', 'after_qc_filtering', 'after_doublet_filtering'))) %>%
    group_by(Sample) %>%
    arrange(Sample, class) %>%
    mutate(Freq=as.numeric(Freq),
             next_remaining = lead(Freq, default = 0),  # Next step's remaining count (or 0 for the last step)
             difference = Freq - next_remaining         # Difference from the next entry
           )
  ggplot(plotdata, aes(x=paste0(Individual, Replicate), y=difference, fill=class)) +
    geom_col(position='stack') +
    geom_hline(data = target_cells, aes(yintercept = value), lty='dashed', color='black') +
    facet_wrap(~Kit, nrow=1, scales='free_x', labeller = labeller(Kit = label_function)) +
    theme_bw() +
    theme(text = element_text(size = 16),
          axis.text.x = element_text(angle=45, vjust=0.5)) + 
    ggbreak::scale_y_break(c(130000, 200000), ticklabels = c(200000, 230000))+ 
    scale_fill_manual(values = c('darkgrey', '#D7191C', '#FDAE61', '#e9e29c', '#39B185'),
                      labels = c('Remainder barcoded', 'Unrecovered cells', 'Low quality cells', 'Multiplets', 'High quality singlet')) +
    labs(x='Sample', y='Number of cells', fill='Classification', caption = 'Dashed lines indicate targeted cell recovery') ->
    cell_recovery_figures[['cell_recovery_counts_break']] 
  ggsave(filename = 'cell_recovery_counts_break.png', path = here('figures/3p/cell_recovery'),
         device = 'png', 
         height = unit(7, 'in'), width = unit(12, 'in'))
  
  ggplot(plotdata, aes(x=paste0(Individual, Replicate), y=difference, fill=class)) +
    geom_col(position='stack') +
    geom_hline(data = target_cells, aes(yintercept = value), lty='dashed', color='black') +
    facet_wrap(~Kit, scales = 'free_x', nrow=1, labeller = labeller(Kit = label_function)) +
    theme_bw() +
    theme(text = element_text(size = 16),
          axis.text.x = element_text(angle=45, vjust=0.5)) + 
    scale_fill_manual(values = c('darkgrey', '#D7191C', '#FDAE61', '#e9e29c', '#39B185'),
                      labels = c('Remainder barcoded', 'Unrecovered cells', 'Low quality cells', 'Multiplets', 'High quality singlet')) +
    labs(x='Sample', y='Number of cells', fill='Classification') ->
    cell_recovery_figures[['cell_recovery_counts']]
  ggsave(filename = 'cell_recovery_counts.png', path = here('figures/3p/cell_recovery'),
         device = 'png',
         height = unit(7, 'in'), width = unit(12, 'in'))
  
  ggplot(plotdata, aes(x=paste0(Individual, Replicate), y=difference, fill=class)) +
    geom_col(position='stack') +
    geom_hline(data = target_cells, aes(yintercept = value), lty='dashed', color='black') +
    facet_wrap(~Kit, scales = 'free', nrow=1, labeller = labeller(Kit = label_function)) +
    theme_bw() +
    theme(text = element_text(size = 16),
          axis.text.x = element_text(angle=45, vjust=0.5)) + 
    scale_fill_manual(values = c('darkgrey', '#D7191C', '#FDAE61', '#e9e29c', '#39B185'),
                      labels = c('Remainder barcoded', 'Unrecovered cells', 'Low quality cells', 'Multiplets', 'High quality singlet')) +
    labs(x='Sample', y='Number of cells', fill='Classification') ->
    cell_recovery_figures[['cell_recovery_counts']] 
  ggsave(filename = 'cell_recovery_counts_freey.png', path = here('figures/3p/cell_recovery'),
         device = 'png', 
         height = unit(7, 'in'), width = unit(16, 'in'))
  
  plotdata %>% 
    group_by(Sample) %>%
    filter(class != '# cells barcoded') %>%
    mutate(difference = 100 * difference / sum (difference)) %>%
    ggplot(aes(x=paste0(Individual, Replicate), y=difference, fill=class)) +
    geom_col(position='stack') +
    theme_bw() +
    theme(text = element_text(size = 16),
          axis.text.x = element_text(angle=45, vjust=0.5)) + 
    geom_hline(data = target_cells_fraction, aes(yintercept = 100*value), lty='dashed', color='black') +
    facet_wrap(~Kit, scales = 'free_x', nrow=1, labeller = labeller(Kit = label_function)) +
    scale_fill_manual(values = c('#D7191C', '#FDAE61', '#e9e29c', '#39B185'),
                      labels = c('Unrecovered cells', 'Low quality cells', 'Multiplets', 'High quality singlet')) +
    labs(x='Sample', y='% of capture', fill='Classification', caption = 'Dashed lines indicate expected cell recovery') ->
    cell_recovery_figures[['cell_recovery_portions_nocellbarcoded']] 
  ggsave(filename = 'cell_recovery_portions_noremainderbarcoded.png', path = here('figures/3p/cell_recovery'),
         device = 'png', 
         height = unit(7, 'in'), width = unit(12, 'in')) 
  
  return(cell_recovery_figures)
}

cell_recovery_figures <- main()

rm(main)