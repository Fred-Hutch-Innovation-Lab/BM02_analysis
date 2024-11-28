library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)       ## parses config yaml
library(Seurat)
# library(patchwork)  ## arrange plots


# main <- function() {
  # fig_objs <- readRDS(here('rds/05_objs_post_clustering.rds'))
  comp_table <- read.table(here('rds/composition_table_coarse.txt'), sep = '\t', header = TRUE)
  kit_order <- read.table(here('config/kit_order.txt'))$V1
  metadata <- read.csv(here('config/metadata.csv')) %>%
    mutate(Kit = factor(Kit, levels = kit_order))
  
  
  # comp_table <-  lapply(objs, function(x) {
  #   y <- tabyl(x$cell_labels.coarse)
  #   colnames(y) <- c('Celltype', 'n', 'percent')
  #   y
  # }) |> rbindlist(idcol = 'Sample') |>
  #   mutate_all(~replace(., is.na(.), 0)) 
  # 
  # ## Adding missing counts
  # if (FALSE){
  #   all_combinations <- expand.grid(Sample = unique(comp_table$Sample),
  #                                   Celltype = unique(comp_table$Celltype))
  #   
  #   comp_table <- merge(all_combinations, comp_table,
  #                       by = c("Sample", "Celltype"), all.x = TRUE)
  #   comp_table$percent[is.na(comp_table$percent)] <- 0
  # }

  comp_table |> 
    melt(measure.vars=unique(coarse_mapping)) |> 
    group_by(Sample) |> 
    mutate(percent = 100*value / sum(value)) |>
    ggplot(aes(x = Kit, y = percent,
               group = Kit, color = Kit, shape = paste0(Individual, Replicate))) + 
    geom_point() +
    labs(x = element_blank(), y = 'Proportion of sample',
         color = 'Kit', shape='Sample',
         caption = paste(
           'Portions are of the total capture for a sample.',
           # '* = significant below 0.05 FDR',
           sep='\n')) +
    theme_bw() +
    scale_shape_manual(values = c(1:5))+ 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    facet_grid(~ factor(variable), scales='free_y')
  ggsave(path= here('figures'), filename='composition_coarse_labels.png', device = 'png', 
         width = unit(13, 'in'), height = unit(4, 'in'))
  