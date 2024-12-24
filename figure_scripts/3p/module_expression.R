library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)       ## parses config yaml
library(Seurat)
library(reshape2)
library(data.table)
library(patchwork)

# Config ----

label_order <- list(
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

fig_objs <- readRDS(here('rds/3p/07_post_module_scoring.rds'))
source(here('config/kit_order.R'))
metadata <- read.csv(here('config/3p/metadata.csv')) %>%
  mutate(Kit = factor(Kit, levels = kit_order_3p))

module_figures <- list()

# Module dotplots ----
## Fine ----
signature.names <- colnames(fig_objs$Flex@meta.data)[grepl('_filtered_UCell$', colnames(fig_objs$Flex@meta.data))]
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
  mutate(Kit = factor(Kit, levels = kit_order_3p), 
         module=factor(module, levels=c(label_order$fine)),
         cell_labels.fine = factor(cell_labels.fine, levels = label_order$fine)) 

module_figures[['module_scores_kit_fine']] <-
  plotdata |> filter(
    as.character(cell_labels.fine) == as.character(module)
  ) |> 
  filter(!module %in% c('B', 'T', 'Monocyte')) |> 
  ggplot(aes(x=Kit, y=module, fill=value, label=round(value, 2))) +
  geom_tile() +
  geom_text() +
  scale_fill_gradient(low = "white", high = "red", limits=c(0,1)) +
  scale_x_discrete(labels = label_function) +
  theme_minimal() +
  labs(x='Kit', y='Marker module', fill='Average\nscore') +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave(plot = module_figures[['module_scores_kit_fine']],
       path= here('figures/3p/module_expression'), filename='modules_fine.png', device = 'png', 
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
  module_figures[['confusion_mat_fine']][[kit]] <- 
    confusion_plot(plotdata, kit) +
    ggtitle(label_function(kit))
  ggsave(plot = module_figures[['confusion_mat_fine']][[kit]],
         path= here('figures/3p/module_expression'), filename=paste0(kit, '_confusion_mat_fine.png'), device = 'png', 
         width = unit(9, 'in'), height = unit(8, 'in'))
}

