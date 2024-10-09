library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)       ## parses config yaml
library(Seurat)

main <- function() {
  fig_objs <- readRDS(here('rds/02-filtered_objs.rds'))
  metadata <- read.csv(here('config/metadata.csv'))
  cell_recovery_figures <- list()
  
  doublets <- lapply(fig_objs, function(x) {
    x$scDblFinder.class %>% 
      table() %>%
      as.data.frame()
  }) %>% 
    data.table::rbindlist(idcol = 'Sample') %>%
    dplyr::rename('class' = '.') %>%
    group_by(Sample)
  
  bad_cells <- lapply(fig_objs, function(x) {
    x@misc$filtering_receipts$original_capture[[2]] - x@misc$filtering_receipts$after_qc_filtering[[2]]
  }) %>% stack()
  colnames(bad_cells) <- c('Freq', 'Sample')
  bad_cells$class <- 'low_quality_cell'
  
  plotdata <- rbind(doublets, bad_cells) %>%
    merge(metadata, by='Sample') %>%
    arrange(Sample, class) %>%
    mutate(class = factor(class, levels = c('low_quality_cell', 'doublet', 'singlet')))
  
  
  ggplot(plotdata, aes(x=paste0(Individual, Replicate), y=Freq, fill=class)) +
    geom_col(position='stack') +
    facet_wrap(~Kit, scales = 'free_x', nrow=1) +
    scale_fill_manual(values = c('darkred', 'goldenrod', 'darkgreen'), labels = c('Low quality cell', 'Multiplet', 'Singlet')) +
    labs(x='Sample', y='Number of cells', fill='Classification') ->
    cell_recovery_figures[['cell_recovery_counts']] 
  ggsave(filename = 'cell_recovery_counts.png', path = here('figures'),
         device = 'png', 
         height = unit(7, 'in'), width = unit(9, 'in'))
  
  plotdata %>% 
    group_by(Sample) %>%
    mutate(Freq = Freq / sum (Freq)) %>%
    ggplot(aes(x=paste0(Individual, Replicate), y=Freq, fill=class)) +
    geom_col(position='stack') +
    facet_wrap(~Kit, scales = 'free_x', nrow=1) +
    scale_fill_manual(values = c('darkred', 'goldenrod', 'darkgreen'), labels = c('Low quality cell', 'Multiplet', 'Singlet')) +
    labs(x='Sample', y='Portion of capture', fill='Classification') ->
    cell_recovery_figures[['cell_recovery_portions']] 
  ggsave(filename = 'cell_recovery_portions.png', path = here('figures'),
         device = 'png', 
         height = unit(7, 'in'), width = unit(9, 'in')) 
  
  return(cell_recovery_figures)
}

cell_recovery_figures <- main()
rm(main)