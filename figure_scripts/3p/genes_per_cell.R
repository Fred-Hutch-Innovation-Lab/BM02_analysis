# Setup ----
source(here('figure_scripts/utils.R'))
library(UpSetR)

# Load data ----
fig_objs <- readRDS(here('rds/3p/02-objs_post_cell_filtering.rds'))
 
# Prepare plotdata -----
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

detected_genes_kit_level <- lapply(unique(metadata_3p$Kit), function(kit) {
  samples <- metadata_3p$Sample[metadata_3p$Kit == kit]
  sets <- list()
  for (sample in samples) {
    sets <- list(sets, detected_genes_sample_level[[sample]])
  }
  result <- table(unlist(sets))
  names(result[result >= 3])
})
names(detected_genes_kit_level) <- unique(metadata_3p$Kit)

# Plot ----
lapply(detected_genes_sample_level, length) %>%
  stack() %>%
  merge(metadata_3p, by.x='ind', by.y='Sample') %>%
  ggplot(aes(x=Kit, y=values)) +
  geom_boxplot() +
  geom_point(aes(shape=paste0(Individual, Replicate))) +
  scale_x_discrete(labels = label_function()) +
  labs(x = 'Kit', y='Genes detected', shape='Sample', caption = 'Genes detected in at least 10 cells') ->
  figures[['usable_genes_sample']]
my_plot_save(image = figures[['usable_genes_sample']], 
             path = here('figures/3p/gene_recovery/usable_genes_sample.svg'), 
             width = 9.5, height = 5)


lapply(detected_genes_kit_level, length) %>%
  stack() %>%
  ggplot(aes(x=ind, y=values)) +
  geom_col() +
  # geom_point(aes(shape=paste0(Individual, Replicate))) +
  # facet_wrap(~ Kit, scales='free_x') + 
  scale_x_discrete(labels = label_function()) +
  labs(x = 'Kit', y='Genes detected', caption = 'Genes detected in at least 10 cells in 3+ samples') ->
  figures[['usable_genes_kit']]
my_plot_save(image = figures[['usable_genes_kit']], 
             path = here('figures/3p/gene_recovery/usable_genes_kit.svg'), 
             width = 9, height = 6)
  
## Euler plot ---- 
plotdata <- detected_genes_kit_level
names(plotdata) <- label_function(names(detected_genes_kit_level))
figures[['gene_overlap']] <-
  plot(eulerr::euler(plotdata, shape='ellipse'))
my_plot_save(image = figures[['gene_overlap']], 
             path = here('figures/3p/gene_recovery/usable_genes_set_overlap.svg'), 
             width = 4, height = 4)
  
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
kits <- levels(metadata_3p$Kit)[levels(metadata_3p$Kit) != 'Flex'] |>
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
svglite::svglite(filename = here('figures/3p/gene_recovery/usable_genes_upset.svg'),
                 width = unit(10, 'in'), height = unit(6, 'in'))
figures[['gene_overlap_upset']]
dev.off()
