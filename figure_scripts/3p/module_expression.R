# Setup ----
source(here('figure_scripts/utils.R'))
# library(reshape2)
# library(data.table)

# Load data ----


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

# Plot functions ----
confusion_plot <- function(plotdata, kit){
  plotdata |>
    filter(Kit == kit) |>
    ggplot(aes(x=cell_labels.fine, y=module, fill=value, label=ifelse(is.na(value), 'NA', value))) +
    geom_tile() +
    geom_text() +
    scale_fill_gradient(low = "white", high = "red", limits=c(0,1)) +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(x='Cell label', y='Marker module', fill='Average\nmarker gene\nmodule score')
}

# Generate plotdata ----
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
         cell_labels.fine = factor(cell_labels.fine, levels = label_order$fine),
         value = round(value, 2)) |>
  complete(Kit, module, cell_labels.fine, fill = list(value = NA)) |>
  filter(cell_labels.fine != 'Unknown')

# Plot ----
plotdata |> filter(
    as.character(cell_labels.fine) == as.character(module)
  ) |> 
  filter(!module %in% c('B', 'T', 'Monocyte')) |> 
  ggplot(aes(x=Kit, y=module, fill=value, label=ifelse(is.na(value), 'NA', value))) +
  geom_tile() +
  geom_text() +
  scale_fill_gradient(low = "white", high = "red", limits=c(0,1)) +
  scale_x_discrete(labels = label_function) +
  labs(x='Kit', y='Celltype', fill='Average\nmarker gene\nmodule score') +
  theme(axis.text.x = element_text(angle=45, hjust=1)) ->
figures[['module_scores_kit_fine']] 

my_plot_save(image = figures[['module_scores_kit_fine']],
       path= here('figures/3p/module_expression/module_scores_kit.svg'), 
       width = 9, height=8)


# Confusion matrices ----
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

for (kit in unique(metadata_3p$Kit)) {
  plots <- list()
  figures[['confusion_mat_fine']][[kit]] <- 
    confusion_plot(plotdata, kit) +
    ggtitle(label_function(kit))
  my_plot_save(image = figures[['confusion_mat_fine']][[kit]],
               path= here('figures/3p/module_expression', paste0(kit, '_confusion_mat_fine.svg')),
               width = 9, height=8)
}
