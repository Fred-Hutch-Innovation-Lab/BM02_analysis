library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)       ## parses config yaml
library(Seurat)
library(UpSetR)

main <- function(figures_list) {
  fig_objs <- readRDS(here('rds/3p/02-objs_post_cell_filtering.rds'))
  source(here('config/kit_order.R'))
  metadata <- read.csv(here('config/3p/metadata.csv')) %>%
    mutate(Kit = factor(Kit, levels = kit_order_3p))
  
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
    scale_x_discrete(labels = label_function) +
    # facet_wrap(~ Kit, scales='free_x') + 
    labs(x = 'Kit', y='Genes detected', shape='Sample', caption = 'Genes detected in at least 10 cells')
  ggsave(plot=figures_list[['usable_genes_sample']], 
         here('figures/3p/gene_recovery/usable_genes_sample.png'), 
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
    scale_x_discrete(labels = label_function) +
    labs(x = 'Kit', y='Genes detected', caption = 'Genes detected in at least 10 cells in 3+ samples')
  ggsave(plot=figures_list[['usable_genes_kit']],
         here('figures/3p/gene_recovery/usable_genes_kit.png'),
         width = unit(7, 'in'), height = unit(6, 'in'))
  
  # Euler plot ---- 
  plotdata <- detected_genes_kit_level
  names(plotdata) <- label_function(names(detected_genes_kit_level))
  figures_list[['gene_overlap']] <-
    plot(eulerr::euler(plotdata, shape='ellipse'))
  ggplot2::ggsave(plot=figures_list[['gene_overlap']], 
                  filename = here('figures/3p/gene_recovery/usable_genes_set_overlap.png'),
                  width = unit(4, 'in'), height = unit(4, 'in'))
  
  # Upset plot ----
  label_function2 <- function(kits) {
    lapply(kits, function(kit) {
      case_when(
        kit == 'NextGEM3P' ~ "NextGEM_3P",
        kit == 'GEMX3P' ~ "GEMX_3P",
        kit == 'Fluent_v4' ~ "Fluent_v4",
        kit == 'Fluent_V' ~ "Fluent_V",
        kit == 'Parse_v3' ~ "Parse_v3",
        .default = kit
      )
    }) |>
      unlist() |>
      unname()
  }
  kits <- levels(metadata$Kit)[levels(metadata$Kit) != 'Flex'] |>
    label_function2()
  intersections <- lapply(1:length(kits), function(x) kits[-x])
  intersections <- c(
    c('Flex', kits),
    list(list('Fluent_v4', 'Fluent_V')),
    list(list("GEMX_3P", "NextGEM_3P")),
    intersections,
    list(kits),
    list(c(kits, 'Flex'))
  )
  plotdata <- detected_genes_kit_level
  names(plotdata) <- label_function2(names(detected_genes_kit_level))
  
  upset(fromList(plotdata), 
        sets=kits,
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
    figures[['gene_overlap_upset']]
  # figures[['gene_overlap_upset']] + scale_y_discrete(labels = label_function)
  # figures[['gene_overlap_upset']]$Set_names <- label_function(figures[['gene_overlap_upset']]$Set_names)
  figures[['gene_overlap_upset']]
  png(here('figures/3p/gene_recovery/usable_genes_upset.png'), width = 10, height = 6, units = 'in', res=300)
  figures[['gene_overlap_upset']]
  dev.off()
  
  return(figures_list)
}
  
if (!exists(quote(figures))) {
  figures <- list()
} else {
  message('Figures object exists in environment, adding figures to existing obj. Be careful of overwrites.')
}

figures <- main(figures)
rm(main)
