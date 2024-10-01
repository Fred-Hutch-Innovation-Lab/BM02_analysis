library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)       ## parses config yaml
library(Seurat)


objs <- readRDS(here('rds/01_raw_objs_10x.rds'))
metadata <- read.csv(here('config/metadata.csv'))


QC_metric_VlnPlot <- function(objs, metric, ylab = 'metric_value', kit_level_cutoffs = kit_level_cutoffs_df, global_cutoffs = global_cutoffs_df) {
  plotdata <- lapply(objs, function(x){unname(x@meta.data[[metric]])}) %>%
    melt(idcol='Sample') %>%
    merge(metadata, by.x='L1', by.y='Sample') 
  plotdata2 <- group_by(plotdata, Kit) %>% 
    summarize(med = median(value)) 
  
  ggplot(plotdata, aes(y=value, x=paste0(Individual, Replicate))) +
    geom_violin(draw_quantiles = 0.5) +
    facet_wrap(~ Kit, scales='free_y', ncol = 1) +
    labs(x='Sample', y = ylab, title = metric,
         # title='', 
         caption = paste(#'Kit median shown in blue',
           #'Sample median shown in black',
           'Global filtering thresholds in yellow',
           'Kit level filtering thresholds in orange', 
           sep = '\n')) +
    coord_flip()
}

QC_metric_VlnPlot(objs[metadata$Sample[metadata$BM_project=='BM02']], 'nFeature_RNA', 'Feature Count') +
  scale_y_log10()
