# Setup ----
source(here('figure_scripts/utils.R'))

# Load data ----
fig_objs <- readRDS('rds/3p/05_merged_objs_post_clustering.rds')
fig_objs <- lapply(fig_objs, function(x){
  x$individual <- gsub('.+_([^_]+)$', '\\1', x$orig.ident)
  x
})

# Plotting functions ----
my_dimplot <- function(obj,
                       reduction='UMAP_allgenes',
                       key='UMAPallgenes_', 
                       dims=c(1,2),
                       group.by, 
                       colors = color_palette$cell_colors,
                       color_label='Cell type', 
                       alpha=1,
                       drop=TRUE,
                       shuffle=TRUE){
  dimdata <- obj@reductions[[reduction]]@cell.embeddings[,dims]
  dimdata[,1] <- scale(dimdata[,1])
  dimdata[,2] <- scale(dimdata[,2])
  plotdata <- merge(dimdata, obj@meta.data, by='row.names')
  plotdata[[group.by]] <- factor(plotdata[[group.by]], levels=names(colors))
  if (shuffle) {
    plotdata <- plotdata[sample(x = 1:nrow(x = plotdata)), ]
  }
  ggplot(plotdata, aes(x=.data[[paste0(key, dims[1])]],
                       y=.data[[paste0(key, dims[2])]],
                       color=.data[[group.by]])) +
    scale_color_manual(values = colors, breaks = names(colors), drop=drop) +
    rasterize(geom_point(size=0.1, alpha=alpha, show.legend = TRUE), dpi=300) +
    labs(x=paste0(reduction, '_', dims[1]), y=paste0(reduction, '_', dims[2]), color=color_label) +
    theme_classic() +
    theme(text = element_text(size=16, hjust = 0),
          axis.text = element_blank(), 
          axis.title = element_blank(),
          axis.ticks = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=3)))
}

# Generate plots ----
for (kit in unique(metadata_3p$Kit)) {
  ## Cell labels only
  figures[[paste0(kit, '_fine')]] <- 
    my_dimplot(fig_objs[[kit]],
               group.by = 'cell_labels.fine') +
    ggtitle(label_function(kit))
  
  my_plot_save(image = figures[[paste0(kit, '_fine')]], 
               path = here('figures/3p/UMAPs', paste0(kit, '_cell_labels.svg')), 
               width = 6, height = 4)
  
  ## Cell labels and invidiual labels
  figures[[paste0(kit, '_fine_&_individual')]] <- 
    my_dimplot(fig_objs[[kit]], group.by = 'cell_labels.fine') +
    my_dimplot(fig_objs[[kit]], group.by = 'individual', 
               colors = color_palette$samples, color_label='Sample', alpha=0.5) +
    ggtitle(label_function(kit))
  my_plot_save(image = figures[[paste0(kit, '_fine_&_individual')]], 
               path = here('figures/3p/UMAPs', paste0(kit, '_cell_and_individual_labels.svg')), 
               width = 12, height = 4)
  figures[[paste0(kit, '_for_merging')]] <- 
    my_dimplot(fig_objs[[kit]], group.by = 'cell_labels.fine', drop=FALSE) +
    my_dimplot(fig_objs[[kit]], group.by = 'individual', 
               colors = color_palette$samples, color_label='Sample', alpha=0.5) 
  ggtitle(label_function(kit))
}
p1 <- 
  figures[['Flex_for_merging']] +
figures[['NextGEM3P_for_merging']] +
figures[['GEMX3P_for_merging']] +
figures[['Fluent_v4_for_merging']] +
figures[['Fluent_V_for_merging']]  +
figures[['Parse_v3_for_merging']] +
figures[['Scale_for_merging']] + 
  plot_layout(guides = 'collect', ncol=1, design = 'AB
              CC
              DD
              EE
              FF
              GG
              HH')
my_plot_save(image = p1, 
             path = here('figures/3p/UMAPs/all_kits_label_&_sample.svg'), 
             width = 20, height = 28)
  