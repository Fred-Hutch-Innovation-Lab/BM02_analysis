library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)       ## parses config yaml
library(Seurat)
library(reshape2)
library(data.table)
library(patchwork)

# Config ----

label_order <- list(
  # coarse = c(
  #   "B", "T", "NK", "Monocyte",
  #   "Dendritic", "pDC",
  #   "Megakaryocyte", "Unknown"),
  fine = c(
    "B","B naive", "B memory",
    "T", "CD4+ T", "CD8+ T", 
    "NK",
    'Granulocyte', 'Erythrocyte',
    "Classical monocyte", "Non-classical monocyte", "Monocyte", 
    "Dendritic", "pDC",
    "Megakaryocyte", 'Unknown'
  )
)

# Load data ----

fig_objs <- readRDS(here('rds/07_post_module_scoring.rds'))
kit_order <- read.table(here('config/kit_order.txt'))$V1
metadata <- read.csv(here('config/metadata.csv')) %>%
  mutate(Kit = factor(Kit, levels = kit_order))

module_figures <- list()

# Module dotplots ----

# ## Coarse ----
# 
signature.names <- colnames(fig_objs$Flex@meta.data)[grepl('_filtered_UCell', colnames(fig_objs$Flex@meta.data))]
# plotdata <- lapply(fig_objs, function(x){
#   x@meta.data |>
#     select(cell_labels.coarse, all_of(signature.names)) |>
#     group_by(cell_labels.coarse) |>
#     summarize(across(everything(),\(x) mean(x, na.rm = TRUE)))
# }) |> 
#   rbindlist(idcol = 'Sample') |>
#   melt() |>
#   mutate(module=gsub('(.+)_filtered_UCell', '\\1', variable)) |>
#   select(-variable) |>
#   mutate(module=factor(module, levels=c(label_order$fine)),
#          cell_labels.coarse = factor(cell_labels.coarse, levels = label_order$coarse))
# 
# module_figures[['module_scores_kit_coarse']] <-
#   plotdata |> filter(
#     as.character(cell_labels.coarse) == as.character(module)
#   ) |> merge(metadata, by='Sample') |>
#   ggplot(aes(x=Kit, y=module, fill=value, label=round(value, 2))) +
#   geom_tile() +
#   geom_text() +
#   scale_fill_gradient(low = "white", high = "red", limits=c(0,1)) +
#   theme_minimal() +
#   labs(x='Kit', y='Marker module', fill='Average\nscore') +
#   facet_wrap(~ Individual + Replicate) +
#   theme(axis.text.x = element_text(angle=45, hjust=1))
# ggsave(plot = module_figures[['module_scores_kit_coarse']],
#        path= here('figures/module_expression'), filename='modules_coarse.png', device = 'png', 
#        width = unit(9, 'in'), height = unit(8, 'in'))
# 

## Fine ----

plotdata <- lapply(fig_objs, function(x){
  x@meta.data |>
    select(cell_labels.fine, all_of(signature.names)) |>
    group_by(cell_labels.fine) |>
    summarize(across(everything(),\(x) mean(x, na.rm = TRUE)))
}) |> 
  rbindlist(idcol = 'Kit') |>
  melt(id.vars=c('Kit','cell_labels.fine')) |>
  mutate(module=gsub('(.+)_filtered_UCell', '\\1', variable)) |>
  select(-variable) |>
  mutate(module=factor(module, levels=c(label_order$fine)),
         cell_labels.fine = factor(cell_labels.fine, levels = label_order$fine))

module_figures[['module_scores_kit_fine']] <-
  plotdata |> filter(
    as.character(cell_labels.fine) == as.character(module)
  ) |> 
  filter(!module %in% c('B', 'T', 'Monocyte')) |> 
  # merge(metadata, by='Sample') |>
  ggplot(aes(x=Kit, y=module, fill=value, label=round(value, 2))) +
  geom_tile() +
  geom_text() +
  scale_fill_gradient(low = "white", high = "red", limits=c(0,1)) +
  theme_minimal() +
  labs(x='Kit', y='Marker module', fill='Average\nscore') +
  # facet_wrap(~ Individual + Replicate) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave(plot = module_figures[['module_scores_kit_fine']],
       path= here('figures/module_expression'), filename='modules_fine.png', device = 'png', 
       width = unit(9, 'in'), height = unit(8, 'in'))


# Confusion matrices ----

confusion_plot <- function(plotdata, kit){
  plotdata |>
    filter(Kit == kit) |>
    ggplot(aes(x=cell_labels.fine, y=module, fill=value, label=round(value, 2))) +
    geom_tile() +
    geom_text() +
    scale_fill_gradient(low = "white", high = "red", limits=c(0,1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(x='Cell label', y='Marker module', fill='Average\nscore')
}

# ## Coarse ----
# 
# 
# plotdata <- lapply(fig_objs, function(x){
#   x@meta.data |>
#     select(cell_labels.coarse, all_of(signature.names)) |>
#     group_by(cell_labels.coarse) |>
#     summarize(across(everything(),\(x) mean(x, na.rm = TRUE)))
# }) |> 
#   rbindlist(idcol = 'Sample') |>
#   melt(id.vars = c('Sample', 'cell_labels.coarse')) |>
#   mutate(module=gsub('(.+)_filtered_UCell', '\\1', variable)) |>
#   select(-variable) |>
#   mutate(module=factor(module, levels=c(label_order$fine)),
#          cell_labels.coarse = factor(cell_labels.coarse, levels = label_order$coarse))
# 
# 
# for (kit in unique(metadata$Kit)) {
#   plots <- list()
#   for (sample in metadata$Sample[metadata$Kit==kit]) {
#     ind <- metadata$Individual[metadata$Sample==sample]
#     rep <- metadata$Replicate[metadata$Sample==sample]
#     plots[[sample]] <- confusion_plot(plotdata, sample) +
#       ggtitle(paste0(ind,rep))
#   }
#   module_figures[['confusion_mat_coarse']][[kit]] <- 
#     plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] +
#     plot_layout(guides = 'collect', axes='keep')
#   ggsave(plot = module_figures[['confusion_mat_coarse']][[kit]],
#          path= here('figures/module_expression'), filename=paste0(kit, '_confusion_mat_coarse.png'), device = 'png', 
#          width = unit(9, 'in'), height = unit(8, 'in'))
# }
# 
# 
# 
# 

## Fine ----


plotdata <- lapply(fig_objs, function(x){
  x@meta.data |>
    select(cell_labels.fine, all_of(signature.names)) |>
    group_by(cell_labels.fine) |>
    summarize(across(everything(),\(x) mean(x, na.rm = TRUE)))
}) |> 
  rbindlist(idcol = 'Kit') |>
  melt(id.vars = c('Kit', 'cell_labels.fine')) |>
  mutate(module=gsub('(.+)_filtered_UCell', '\\1', variable)) |>
  select(-variable) |>
  mutate(module=factor(module, levels=c(label_order$fine)),
         cell_labels.fine = factor(cell_labels.fine, levels = label_order$fine))

for (kit in unique(metadata$Kit)) {
  plots <- list()
  # for (sample in metadata$Kit[metadata$Kit==kit]) {
  #   ind <- metadata$Individual[metadata$Kit==sample]
  #   rep <- metadata$Replicate[metadata$Kit==sample]
    # plots[[sample]] <- confusion_plot(plotdata, kit) +
    #   ggtitle(paste0(kit))
  # }
  module_figures[['confusion_mat_fine']][[kit]] <- 
    confusion_plot(plotdata, kit) +
    ggtitle(paste0(kit))
    # plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] +
    # plot_layout(guides = 'collect', axes='keep')
  ggsave(plot = module_figures[['confusion_mat_fine']][[kit]],
         path= here('figures/module_expression'), filename=paste0(kit, '_confusion_mat_fine.png'), device = 'png', 
         width = unit(9, 'in'), height = unit(8, 'in'))
}

