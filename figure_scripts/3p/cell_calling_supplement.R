# Setup ----
source(here('figure_scripts/utils.R'))
library(janitor)
library(gt)
library(ggforestplot) ## Striped lines for dotplots

# Load data ----
fig_objs <- readRDS('rds/3p/05_merged_objs_post_clustering.rds')

# Plotting functions ----
my_dimplot <- function(obj,
                       reduction='UMAP_allgenes',
                       key='UMAPallgenes_', 
                       dims=c(1,2),
                       group.by, 
                       colors,
                       color_label='Cell type', 
                       alpha=1,
                       shuffle=TRUE,
                       title){
  dimdata <- obj@reductions[[reduction]]@cell.embeddings[,dims]
  dimdata[,1] <- scale(dimdata[,1])
  dimdata[,2] <- scale(dimdata[,2])
  plotdata <- merge(dimdata, obj@meta.data, by='row.names')
  # plotdata[[group.by]] <- factor(plotdata[[group.by]], levels=names(colors))
  if (shuffle) {
    plotdata <- plotdata[sample(x = 1:nrow(x = plotdata)), ]
  }
  ggplot(plotdata, aes(x=.data[[paste0(key, dims[1])]],
                       y=.data[[paste0(key, dims[2])]],
                       color=.data[[group.by]])) +
    rasterize(geom_point(size=0.1, alpha=alpha, show.legend = TRUE), dpi=300) +
    (if (!missing(colors)) {scale_color_manual(values = colors, breaks = names(colors), drop=TRUE)} else NULL) +
    labs(x=paste0(reduction, '_', dims[1]), y=paste0(reduction, '_', dims[2]), color=color_label, title = title) +
    theme(axis.text = element_blank(), 
          axis.title = element_blank(),
          axis.ticks = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=3))) 
}
my_dotplot <-  function(x) {
  z <- DotPlot(x, celltype_markers, group.by = 'clusters') + 
    theme_grey() +
    theme(axis.text.x = element_text(angle=45, hjust=1), text = element_text(size=16, family = 'Arial'))
  geom_stripes() 
  z$data$feature.groups <- factor(z$data$feature.groups,
                                  levels=gtools::mixedsort(as.character(unique(z$data$feature.groups))))
  z
}
my_table <- function(x) {
  x@meta.data$cell_labels.fine |> 
    tabyl() |> 
    arrange(desc(n)) |>
    adorn_rounding(digits=2) |>
    # knitr::kable(col.names=c("Cell label", "Count", "Proportion"))
    gt() |>
    cols_label(`x@meta.data$cell_labels.fine`="Cell label", n="Count", percent="Proportion") |>
    opt_table_font(size = 40) |> 
    wrap_table() +
    theme(plot.margin = unit(c(60,0,0,60), 'pt'))
}

# Prepare plotdata ----
## Dotplots ----
celltype_markers <- read.csv('/fh/fast/_IRC/FHIL/grp/FHIL_knowledgebase/biology/celltype_markers.csv')
parse_marker_table <- function(celltype_dataframe) {
  celltype_dataframe %>%
    filter(expression_level == 'Increased' & tissue == 'PBMC') %>%
    filter(confidence == 'high' | 
             (celltype %in% c('Dendritic','Monocyte', 
                              'Erythrocyte', 'Granulocyte',
                              'B naive', 'B memory') 
              & confidence == 'med')) %>%
    filter(celltype %in% c('T', 'CD4+ T', 'CD8+ T', 
                           'B', 'B naive', 'B memory',
                           'Monocyte', 'Non-classical monocyte', 'Classical monocyte',
                           'NK', 'Dendritic', 'pDC', 'Neutrophil',  
                           'Erythrocyte', 'Granulocyte',
                           'Megakaryocyte', 'Lymphocyte progenitor', 'HSPC'
    )) %>%
    select(gene_symbol, celltype) %>%
    unstack()
}

prune_marker_list <- function(input_list) {
  all_items <- unlist(input_list, use.names = TRUE)
  dup_items <- names(table(all_items))[table(all_items) > 1]
  input_list_cleaned <- lapply(input_list, function(x) setdiff(x, dup_items))
  
  new_groups <- sapply(dup_items, function(item) {
    found_in <- names(input_list)[sapply(input_list, function(x) item %in% x)]
    new_group_name <- paste(found_in, collapse = " &\n")
    return(new_group_name)
  }, USE.NAMES = TRUE)
  
  new_list <- input_list_cleaned  
  for (i in seq_along(dup_items)) {
    new_list[[new_groups[i]]] <- c(new_list[[new_groups[i]]], dup_items[i])
  }
  
  wrap_name <- function(name) {
    wrapped_name <- strwrap(name, width = 8, simplify = TRUE)
    wrapped_name <- paste(wrapped_name, collapse = "\n")
    return(wrapped_name)
  }
  
  names(new_list) <- sapply(names(new_list), wrap_name)# {
  new_list
}

celltype_markers <- parse_marker_table(celltype_markers)
celltype_markers <- prune_marker_list(celltype_markers)

# Plot ----
for (kit in unique(metadata_3p$Kit)) {
  figures[[paste0(kit, '_fine')]] <- 
    ((my_dimplot(fig_objs[[kit]], group.by = 'renamed.predicted.pbmcsca_seurat_annotations', title='pbmc_sca') +
     my_dimplot(fig_objs[[kit]], group.by = 'renamed.predicted.pbmc3k_seurat_annotations', title= 'pbmc_3k') +
     my_dimplot(fig_objs[[kit]], group.by = 'renamed.Mona.main.labels', title= 'Monaco') +
     my_dimplot(fig_objs[[kit]], group.by = 'renamed.HPCA.main.labels', title='HPCA') +
       plot_layout(ncol=4)) /
    (my_dimplot(fig_objs[[kit]], group.by = 'clusters', title='Leiden clustering') +
     my_dimplot(fig_objs[[kit]], group.by = 'cell_labels.fine', colors = color_palette$cell_colors,  title='Final annotations') +
     my_table(fig_objs[[kit]])
     ) /
    my_dotplot(fig_objs[[kit]])) +
    plot_annotation(label_function(kit))
  
  my_plot_save(image = figures[[paste0(kit, '_fine')]], 
               path = here('figures/3p/cell_calling_supplement', paste0(kit, '_cell_calling_supplement.svg')), 
               width = 25, height = 18)
}
