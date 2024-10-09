library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)       ## parses config yaml
library(ggrastr)    ## rasterized point layers for reduced image size
library(Seurat)
library(reshape2)  ## melt

## wrapped in function to avoid filling environment if sourced
main <- function(figures_list, filtering_stage = 'after', outdir='qc_vln_plots') {
  
  QC_metric_VlnPlot <- function(objs, metric, ylab = 'metric_value', caption=NA) {
    plotdata <- lapply(objs, function(x){unname(x@meta.data[[metric]])}) %>%
      melt(idcol='Sample') %>%
      merge(metadata, by.x='L1', by.y='Sample') 
    plotdata2 <- group_by(plotdata, Kit) %>% 
      summarize(med = median(value)) 
    
    ggplot(plotdata, aes(y=value, x=paste0(Individual, Replicate))) +
      rasterize(geom_jitter(alpha=0.1, size=.1), dpi=300) +
      geom_violin(draw_quantiles = 0.5) +
      # theme_bw() +
      facet_wrap(~ Kit, nrow = 1) +
      labs(x='Sample', y = ylab, caption=caption) #+
    # coord_flip()
  }
  
  if (filtering_stage == 'after') {
    fig_objs <- readRDS(here('rds/02-filtered_objs.rds'))
    caption <- 'QC metrics per cell on data after filtering by global thresholds'
  } else if (filtering_stage == 'before') {
    fig_objs <- readRDS(here('rds/01_raw_objs_10x.rds'))
    caption <- 'QC metrics per cell on data output from pipelines without additional filtering'
  }
  metadata <- read.csv(here('config/metadata.csv'))
  
  labels <- data.frame(
    metric = c('nFeature_RNA', 'nCount_RNA', 'rbRatio', 'mtRatio'),
    label  = c('Feature count', 'UMI count', 'Ribosomal read portion', 'Mitochondrial read portion'),
    logscale = c(TRUE, TRUE, FALSE, FALSE)
  )
  
  for (i in 1:nrow(labels)) {
    metric <- labels[i, 'metric']
    label <- labels[i, 'label']
    plt <- QC_metric_VlnPlot(fig_objs, metric, label, caption)
    if (labels$logscale[i]) {
      plt <- plt + scale_y_log10()
    }
    filename <- paste0(c(metric, 'vln'), collapse='_')
    figures_list[[filename]] <- plt
    ggsave(plt, filename = paste0(filename, '.png'), path = here(file.path('figures', outdir)),
           device = 'png', 
           height = unit(7, 'in'), width = unit(15, 'in'))
  }
  return(figures_list)
}

if (!exists(quote(figures))) {
  figures <- list()
} else {
  message('Figures object exists in environment, adding figures to existing obj. Be careful of overwrites.')
}

figures <- main(figures, filtering_stage = 'after', outdir='qc_vln_plots_after_filtering')

rm(main)
