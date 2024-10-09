library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)       ## parses config yaml
library(Seurat)

main <- function(figures_list) {
  fig_objs <- readRDS(here('rds/02-filtered_objs.rds'))
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
  
  
  figures_list[['usable_genes_sample']] <-
    lapply(detected_genes_sample_level, length) %>%
    stack() %>%
    merge(metadata, by.x='ind', by.y='Sample') %>%
    ggplot(aes(x=Kit, y=values)) +
    geom_boxplot() +
    geom_point(aes(shape=paste0(Individual, Replicate))) +
    # facet_wrap(~ Kit, scales='free_x') + 
    labs(x = 'Kit', y='Genes detected', shape='Sample', caption = 'Genes detected in at least 10 cells')
  ggsave(plot=figures_list[['usable_genes_sample']], here('figures/usable_genes_sample.png'), width = unit(7, 'in'), height = unit(5, 'in'))
  
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
  ggsave(plot=figures_list[['usable_genes_kit']], here('figures/usable_genes_kit.png'), width = unit(7, 'in'), height = unit(6, 'in'))
  
  # Euler plot ---- 
  figures_list[['gene_overlap']] <-
    plot(eulerr::euler(detected_genes_kit_level, shape='ellipse'))
  # png(filename = here('figures/usable_genes_set_overlap.png'), width = unit(4, 'in'), height = unit(4, 'in'))
  ggplot2::ggsave(plot=figures_list[['gene_overlap']], filename = here('figures/usable_genes_set_overlap.png'), width = unit(4, 'in'), height = unit(4, 'in'))
  # dev.off()
  # ggsave()
  
  return(figures_list)
}
  
if (!exists(quote(figures))) {
  figures <- list()
} else {
  message('Figures object exists in environment, adding figures to existing obj. Be careful of overwrites.')
}

figures <- main(figures)
rm(main)
