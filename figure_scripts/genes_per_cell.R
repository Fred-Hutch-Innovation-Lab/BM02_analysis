library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)       ## parses config yaml
library(Seurat)
library(UpSetR)

main <- function(figures_list) {
  fig_objs <- readRDS(here('rds/02-objs_post_cell_filtering.rds'))
  metadata <- read.csv(here('config/metadata.csv'))
  
  # Sample levels -----
  detected_genes_sample_level <- lapply(fig_objs, function(obj) {
    genes <- obj@assays$RNA@features@.Data 
    genes <- genes %>% 
      as.data.frame() %>%
      filter(counts == TRUE) %>% 
      rownames()
    detected_genes <- rowSums(obj@assays$RNA@layers$counts > 0) > 10 
    detected_genes <- genes[detected_genes]
    detected_genes
  })
  
  figures_list <- list()
  figures_list[['usable_genes_sample']] <-
    lapply(detected_genes_sample_level, length) %>%
    stack() %>%
    merge(metadata, by.x='ind', by.y='Sample') %>%
    ggplot(aes(x=Kit, y=values)) +
    geom_boxplot() +
    geom_point(aes(shape=paste0(Individual, Replicate))) +
    # facet_wrap(~ Kit, scales='free_x') + 
    labs(x = 'Kit', y='Genes detected', shape='Sample', caption = 'Genes detected in at least 10 cells')
  ggsave(plot=figures_list[['usable_genes_sample']], 
         here('figures/gene_recovery/usable_genes_sample.png'), 
         width = unit(7, 'in'), height = unit(5, 'in'))
  
   # Kit level ----
  detected_genes_kit_level <- lapply(unique(metadata$Kit), function(kit) {
    samples <- metadata$Sample[metadata$Kit == kit]
    sets <- list()
    for (sample in samples) {
      sets <- list(sets, detected_genes_sample_level[[sample]])
    }
    result <- table(unlist(sets))
    names(result[result >= 3])
  })
  names(detected_genes_kit_level) <- unique(metadata$Kit)
  
  figures_list[['usable_genes_kit']] <-
    lapply(detected_genes_kit_level, length) %>%
    stack() %>%
    ggplot(aes(x=ind, y=values)) +
    geom_col() +
    # geom_point(aes(shape=paste0(Individual, Replicate))) +
    # facet_wrap(~ Kit, scales='free_x') + 
    labs(x = 'Kit', y='Genes detected', caption = 'Genes detected in at least 10 cells in 3+ samples')
  ggsave(plot=figures_list[['usable_genes_kit']],
         here('figures/gene_recovery/usable_genes_kit.png'),
         width = unit(7, 'in'), height = unit(6, 'in'))
  
  # Euler plot ---- 
  figures_list[['gene_overlap']] <-
    plot(eulerr::euler(detected_genes_kit_level, shape='ellipse'))
  ggplot2::ggsave(plot=figures_list[['gene_overlap']], 
                  filename = here('figures/gene_recovery/usable_genes_set_overlap.png'),
                  width = unit(4, 'in'), height = unit(4, 'in'))
  
  # Upset plot ----
  kits <- levels(metadata$Kit)[levels(metadata$Kit) != 'Flex']
  intersections <- lapply(1:length(kits), function(x) kits[-x])
  intersections <- c(
    levels(metadata$Kit),
    list(list('Fluent_v4', 'Fluent_V')),
    list(list('GEMX3P', 'NextGEM3P')),
    intersections,
    list(kits),
    list(levels(metadata$Kit))
  )
  upset(fromList(detected_genes_kit_level), 
        sets=levels(metadata$Kit),
        # nsets=8,
        # nintersects = 140, 
        keep.order = TRUE,
        intersections = intersections,
        order.by='degree',
        # cutoff=100,
        decreasing = FALSE,
        mainbar.y.label='Gene count',
        number.angles = 0,
        empty.intersections=FALSE) ->
    figures_list[['gene_overlap_upset']]
  figures_list[['gene_overlap_upset']]
  ggplot2::ggsave(plot=figures_list[['gene_overlap_upset']], 
                  filename = here('figures/gene_recovery/usable_genes_set_overlap.png'),
                  width = unit(12, 'in'), height = unit(5, 'in'))
  
  return(figures_list)
}
  
if (!exists(quote(figures))) {
  figures <- list()
} else {
  message('Figures object exists in environment, adding figures to existing obj. Be careful of overwrites.')
}

figures <- main(figures)
rm(main)
