library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)       ## parses config yaml
library(Seurat)
library(patchwork)  ## arrange plots


main <- function() {
  fig_objs <- readRDS(here('rds/05_objs_post_clustering.rds'))
  kit_order <- read.table(here('config/kit_order.txt'))$V1
  metadata <- read.csv(here('config/metadata.csv')) %>%
    mutate(Kit = factor(Kit, levels = kit_order))
    
  umap_figures <- list()
  
  cell_colors <- list(coarse=c(     
    "T" = "blue",
    "B" = "red",
    "Monocyte" = "purple",
    "NK" = "#4DAF4A",
    "Megakaryocyte" = "#A65628",
    "Dendritic" = "orange",
    "pDC" = "gold",
    'Unknown' = 'grey'
    ), fine = c(
    'T' = "blue",
    'CD8+ T' = "#377EB8",
    'CD4+ T' = "#80B1D3",
    'B naive' = "#FB8072",
    'B memory' = "#E41A1C",
    'Monocyte' = "purple",
    'Classical monocyte' = "#BC80BD",
    'Non-classical monocyte' = "#984EA3",
    'NK' = "#4DAF4A",
    'Megakaryocyte' = "#A65628",
    'Dendritic' =  "orange", 
    'pDC' = "gold", 
    'Unknown' = 'grey'
    )
  ) 
  my_dimplot <- function(obj, reduction='umap', dims=c(1,2), group.by, colors = cell_colors$fine){
    dimdata <- obj@reductions[[reduction]]@cell.embeddings[,dims]
    dimdata[,1] <- scale(dimdata[,1])
    dimdata[,2] <- scale(dimdata[,2])
    plotdata <- merge(dimdata, obj@meta.data, by='row.names')
    plotdata[[group.by]] <- factor(plotdata[[group.by]], levels=names(colors))
    ggplot(plotdata, aes(x=.data[[paste0(reduction, '_', dims[1])]],
                         y=.data[[paste0(reduction, '_', dims[2])]],
                         color=.data[[group.by]])) +
      geom_point(size=0.1, show.legend = TRUE) +
      scale_color_manual(values = colors, breaks = names(colors), drop=FALSE) +
      labs(x=paste0(reduction, '_', dims[1]), y=paste0(reduction, '_', dims[2]), color='Cell type') +
      theme_classic() +
      theme(axis.text = element_blank(), 
            axis.title = element_blank(),
            axis.ticks = element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=3))) +
      ggtitle(gsub('.+_(.+)$', '\\1', obj@project.name))
      
  }
  
  extras <- c("#CCEBC5",  "#FB8072",  "#FCCDE5", "#4DAF4A")
  for (kit in unique(metadata$Kit)) {
    samples <- metadata$Sample[metadata$Kit==kit]
    ps <- lapply(samples, function(sample){
      my_dimplot(fig_objs[[sample]], group.by = 'cell_labels.fine')
    })
    umap_figures[[paste0(kit, '_fine')]] <- ps[[1]] + ps[[2]] + ps[[3]] + ps[[4]] + 
      plot_layout(guides='collect')+ 
      plot_annotation(title = kit,
                      theme = theme(plot.title = element_text(size = 24, hjust = 0.5, face='bold')))
    ggsave(plot = umap_figures[[paste0(kit, '_fine')]],
           path= here('figures/UMAPs'), filename=paste0(kit, '_fine_labels.png'), device = 'png', 
           width = unit(10, 'in'), height = unit(8, 'in'), )
  }
  for (kit in unique(metadata$Kit)) {
    samples <- metadata$Sample[metadata$Kit==kit]
    ps <- lapply(samples, function(sample){
      my_dimplot(fig_objs[[sample]], group.by = 'cell_labels.coarse', colors = cell_colors$coarse)
    })
    umap_figures[[paste0(kit, '_coarse')]] <- ps[[1]] + ps[[2]] + ps[[3]] + ps[[4]] + 
      plot_layout(ncol = 2, guides='collect') + 
      plot_annotation(title = kit, 
                      theme = theme(plot.title = element_text(size = 24, hjust = 0.5, face='bold')))
    ggsave(plot = umap_figures[[paste0(kit, '_coarse')]],
           path= here('figures/UMAPs'), filename=paste0(kit, '_coarse_labels.png'), device = 'png', 
           width = unit(10, 'in'), height = unit(8, 'in'), )
  }
  return(umap_figures)
}

umap_figures <- main()
rm(main)