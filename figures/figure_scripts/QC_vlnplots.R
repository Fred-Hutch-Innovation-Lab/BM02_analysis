library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)       ## parses config yaml
library(Seurat)

QC_metric_VlnPlot <- function(objs, metric, ylab = 'metric_value') {
  plotdata <- lapply(objs, function(x){unname(x@meta.data[[metric]])}) %>%
    melt(idcol='Sample') %>%
    merge(metadata, by.x='L1', by.y='Sample') 
  plotdata2 <- group_by(plotdata, Kit) %>% 
    summarize(med = median(value)) 
  
  ggplot(plotdata, aes(y=value, x=paste0(Individual, Replicate))) +
    geom_violin(draw_quantiles = 0.5) +
    # theme_bw() +
    facet_wrap(~ Kit, scales='free_y', ncol = 1) +
    labs(x='Sample', y = ylab) +
    coord_flip()
}

main <- function(objs){
  metadata <- read.csv(here('config/metadata.csv'))
  
  qc_figures <- list()
  
  labels <- data.frame(
    metric = c('nFeature_RNA', 'nCount_RNA', 'rbRatio', 'mtRatio'),
    label  = c('Feature count', 'UMI count', 'Ribosomal read portion', 'Mitochondrial read portion'),
    logscale = c(TRUE, TRUE, FALSE, FALSE)
  )
  
  for (bm in unique(metadata$BM_project)) {
    samples <- metadata$Sample[metadata$BM_project==bm]
    for (i in 1:nrow(labels)) {
      metric <- labels[i, 'metric']
      label <- labels[i, 'label']
      plt <- QC_metric_VlnPlot(objs[samples], metric, label)
      if (labels$logscale[i]) {
        plt <- plt + scale_y_log10()
      }
      filename <- paste0(c(metric, 'vln'), collapse='_')
      qc_figures[[filename]] <- plt
      ggsave(filename = paste0(filename, '.pdf'), path = here(file.path('figures', bm, 'qc_vln_plots')),
             device = 'pdf', 
             height = unit(2+0.3*length(samples), 'in'), width = unit(7, 'in'))
    }
  }
  return(qc_figures)
}

# objs <- readRDS(here('rds/01_raw_objs_10x.rds'))
qc_figures <- main(objs)

