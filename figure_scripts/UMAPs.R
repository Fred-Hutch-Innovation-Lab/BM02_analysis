library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)       ## parses config yaml
library(Seurat)
library(patchwork)  ## arrange plots


main <- function() {
  fig_objs <- readRDS('rds/05_merged_objs_post_clustering.rds')
  fig_objs <- lapply(fig_objs, function(x){
    x$individual <- gsub('.+_([^_]+)$', '\\1', x$orig.ident)
    x
  })
  kit_order <- read.table(here('config/kit_order.txt'))$V1
  metadata <- read.csv(here('config/metadata.csv')) %>%
    mutate(Kit = factor(Kit, levels = kit_order))
  
  ## returns color_palette obj
  source(here('config/color_palette.R'))
  
  umap_figures <- list()
  
  my_dimplot <- function(obj,
                         reduction='UMAP_allgenes',
                         key='UMAPallgenes_', 
                         dims=c(1,2),
                         group.by, 
                         colors = color_palette$cell_colors,
                         color_label='Cell type', 
                         alpha=1,
                         shuffle=TRUE){
    dimdata <- obj@reductions[[reduction]]@cell.embeddings[,dims]
    dimdata[,1] <- scale(dimdata[,1])
    dimdata[,2] <- scale(dimdata[,2])
    plotdata <- merge(dimdata, obj@meta.data, by='row.names')
    plotdata[[group.by]] <- factor(plotdata[[group.by]], levels=names(colors))
    if (shuffle) {
      plotdata <- plotdata[sample(x = 1:nrow(x = plotdata)), ]
    }
    ggplot(plotdata, aes(x=.data[[paste0(key, dims[1])]],
                         y=.data[[paste0(key, dims[2])]],
                         color=.data[[group.by]])) +
      geom_point(size=0.1, alpha=alpha, show.legend = TRUE) +
      scale_color_manual(values = colors, breaks = names(colors), drop=TRUE) +
      labs(x=paste0(reduction, '_', dims[1]), y=paste0(reduction, '_', dims[2]), color=color_label) +
      theme_classic() +
      theme(axis.text = element_blank(), 
            axis.title = element_blank(),
            axis.ticks = element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=3))) #+
      # ggtitle(gsub('.+_(.+)$', '\\1', obj@\))
      
  }
  
  
  for (kit in unique(metadata$Kit)) {
    umap_figures[[paste0(kit, '_fine')]] <- my_dimplot(fig_objs[[kit]], group.by = 'cell_labels.fine') +
      ggtitle(kit)
    ggsave(plot = umap_figures[[paste0(kit, '_fine')]],
           path= here('figures/UMAPs'), filename=paste0(kit, '_fine_labels.png'), device = 'png', 
           width = unit(6, 'in'), height = unit(4, 'in'), )
  }
  
  for (kit in unique(metadata$Kit)) {
    umap_figures[[paste0(kit, '_fine')]] <- 
      my_dimplot(fig_objs[[kit]], group.by = 'cell_labels.fine') +
      my_dimplot(fig_objs[[kit]], group.by = 'individual', 
                 colors = color_palette$samples, color_label='Sample', alpha=0.5)
    ggsave(plot = umap_figures[[paste0(kit, '_fine')]],
           path= here('figures/UMAPs'), filename=paste0(kit, '_sample_&_label.png'), device = 'png', 
           width = unit(12, 'in'), height = unit(4, 'in'), )
  }
  return(umap_figures)
}

umap_figures <- main()
rm(main)