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
  
  cell_colors <- list(
    'T' = "#E41A1C",
    'B' = "#377EB8",
    'Dendritic' = "#4DAF4A",
    'Monocyte' = "#984EA3",
    'NK' = "#FF7F00",
    'Platelet' = "#A65628",
    'Unknown' = "darkgrey"
  )
  for (kit in unique(metadata$Kit)) {
    # for (sample in metadata$Sample[metadata$Kit==kit]) {
    samples <- metadata$Sample[metadata$Kit==kit]
    ps <- lapply(samples, function(sample){
      ind <- metadata$Individual[metadata$Sample==sample]
      rep <- metadata$Replicate[metadata$Sample==sample]
      DimPlot(fig_objs[[sample]], group.by = 'cell_labels.coarse', label=FALSE, order=rev(names(cell_colors))) + 
        scale_color_manual(values = cell_colors, drop=FALSE, show.legend=TRUE) +
        ggtitle(paste0(ind, rep)) + 
        theme(axis.text = element_blank(), 
              axis.title = element_blank(),
              axis.ticks = element_blank())
    })
    umap_figures[[kit]] <- ps[[1]] + ps[[2]] + ps[[3]] + ps[[4]] + 
      plot_layout(ncol = 2, guides='collect') + 
      plot_annotation(title = kit, 
                      theme = theme(plot.title = element_text(size = 24, hjust = 0.5, face='bold')))
    ggsave(plot = umap_figures[[kit]],
           path= here('figures/UMAPs'), filename=paste0(kit, '_coarse_labels.png'), device = 'png', 
           width = unit(9, 'in'), height = unit(8, 'in'), )
  }
  return(umap_figures)
}

umap_figures <- main()
rm(main)