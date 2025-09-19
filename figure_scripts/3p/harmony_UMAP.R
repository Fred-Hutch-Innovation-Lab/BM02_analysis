# Setup ----
source(here('figure_scripts/utils.R'))

# Load data ----
cell_types <- readRDS('rds/3p/05_merged_objs_post_clustering.rds')
cell_types <- lapply(cell_types, \(x){
  x$cell_labels.fine |>
    as.data.frame() |>
    rownames_to_column('cell_id') |>
    rename('cell_labels.fine'=`x$cell_labels.fine`)
}) |> rbindlist()
harmony_umap_loadings <- readRDS(here('rds/3p/harmony_umap_loadings.rds'))
plotdata <- 
  harmony_umap_loadings@cell.embeddings |>
    as.data.frame() |>
  rownames_to_column('cell_id') |>
  # separate_wider_delim(cols=cell_id, names=c('x', 'cell_id'), delim='_', too_many='merge') #|>
  mutate(Kit=case_when(
    grepl('Flex', cell_id) ~ 'Flex',
    grepl('NextGEM', cell_id) ~ 'NextGEM3P',
    grepl('GEMX', cell_id) ~ 'GEMX3P',
    grepl('Fluent_v4', cell_id) ~ 'Fluent_v4',
    grepl('Fluent_V', cell_id) ~ 'Fluent_V',
    grepl('Parse_v3', cell_id) ~ 'Parse_v3',
    grepl('Scale', cell_id) ~ 'Scale',
  ), sample=case_when(
    grepl('F1A', cell_id) ~ 'F1A',
    grepl('F1B', cell_id) ~ 'F1B',
    grepl('F5A', cell_id) ~ 'F5A',
    grepl('F5B', cell_id) ~ 'F5B'
  )) |>
  mutate(cell_id = case_when(
    grepl('Flex', cell_id) ~ gsub(paste0('^Flex_'), '', cell_id),
    grepl('NextGEM', cell_id) ~ gsub(paste0('^NextGEM3P_'), '', cell_id),
    grepl('GEMX', cell_id) ~ gsub(paste0('^GEMX3P_'), '', cell_id),
    grepl('Fluent_v4', cell_id) ~ gsub(paste0('^Fluent_v4_'), '', cell_id),
    grepl('Fluent_V', cell_id) ~ gsub(paste0('^Fluent_V_'), '', cell_id),
    grepl('Parse_v3', cell_id) ~ gsub(paste0('^Parse_v3_'), '', cell_id),
    grepl('Scale', cell_id) ~ gsub(paste0('^Scale_'), '', cell_id)
  )) |>
  mutate(Sample = paste0(Kit, '_', sample)) |>
  select(-c(Kit)) |>
  merge(select(metadata_3p, -c(Sensitivity, File)), by='Sample') 
plotdata |> head()  

my_dimplot <- function(dimdata,
                       group.by, 
                       colors = color_palette$cell_colors,
                       color_label='Cell type', 
                       alpha=1,
                       drop=TRUE,
                       shuffle=TRUE, 
                       big=FALSE){
  dimdata$umap_1 <- scale(dimdata$umap_1)
  dimdata$umap_2 <- scale(dimdata$umap_2)
  plotdata[[group.by]] <- factor(plotdata[[group.by]], levels=names(colors))
  if (shuffle) {
    plotdata <- plotdata[sample(x = 1:nrow(x = plotdata)), ]
  }
  ggplot(plotdata, aes(x=.data[['umap_1']],
                       y=.data[['umap_2']],
                       color=.data[[group.by]])) +
    scale_color_manual(values = colors, breaks = names(colors), drop=drop, labels = label_function()) +
    rasterize(geom_point(size=0.001, alpha=alpha, show.legend = TRUE), dpi=300) +
    labs(x=NULL, y=NULL, color=color_label) +
    theme_classic() +
    theme(text = element_text(size=16, hjust = 0),
          axis.text = element_blank(), 
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 30, face = "bold")) +
    guides(colour = guide_legend(override.aes = list(size=ifelse(big, 10, 3)))) +
    (if (big) {theme(legend.title = element_text(size = 30), 
                     legend.text  = element_text(size = 24),
                     legend.key.size = unit(0.5, "in"))} else NULL)
}

figures[['umap_harmony_kit']] <- my_dimplot(plotdata, group.by='Kit', colors = color_palette$kits, color_label = 'Kit')
figures[['umap_harmony_sample']] <- my_dimplot(plotdata, group.by='sample', colors = color_palette$samples, color_label = 'Sample')
plotdata <- merge(plotdata, cell_types, by= 'cell_id')
figures[['umap_harmony_cell']] <- my_dimplot(plotdata, group.by='cell_labels.fine', colors = color_palette$cell_colors, color_label = 'Cell type')

my_plot_save(image = figures[['umap_harmony_kit']],
             path = here('figures/3p/UMAPs', 'integrated_kit.svg'),
             width = 6.5, height = 5)

my_plot_save(image = figures[['umap_harmony_sample']],
             path = here('figures/3p/UMAPs', 'integrated_sample.svg'),
             width = 6.5, height = 5)

my_plot_save(image = figures[['umap_harmony_cell']],
             path = here('figures/3p/UMAPs', 'integrated_cell.svg'),
             width = 7.5, height = 5)

