# Setup ----
source(here('figure_scripts/utils.R'))

# Load data ----
fig_objs <- readRDS(here('rds/3p/02-objs_post_cell_filtering.rds'))
caption <- NULL #QC metrics per cell on data after filtering by global thresholds'

# Plotting functions ----
QC_metric_VlnPlot <- function(objs, metric, ylab = 'metric_value', caption=NA) {
  plotdata <- lapply(objs, function(x){unname(x@meta.data[[metric]])}) %>%
    melt(idcol='Sample') %>%
    merge(metadata_3p, by.x='L1', by.y='Sample') 
  plotdata2 <- group_by(plotdata, Kit) %>% 
    summarize(med = median(value)) 
  
  ggplot(plotdata, aes(y=value, x=paste0(Individual, Replicate), fill = Kit)) +
    rasterize(geom_jitter(alpha=0.1, size=.1), dpi=100) +
    geom_violin(draw_quantiles = 0.5) +
    scale_fill_manual(values = color_palette$kits, labels = label_function) +
    theme(text = element_text(size = 20),
          panel.spacing=unit(0, "lines"),
          axis.text.x = element_text(angle=45, vjust=0.5)) + 
    facet_wrap(~ Kit, nrow = 1, scales='free_x', labeller = labeller(Kit = label_function)) +
    labs(x='Sample', y = ylab, caption=caption) #+
  # coord_flip()
}
# Plotting ----
labels <- data.frame(
  metric = c('nFeature_RNA', 'nCount_RNA', 'rbRatio', 'mtRatio'),
  label  = c('Feature count', 'UMI count', 'Ribosomal read %', 'Mitochondrial read %'),
  logscale = c(FALSE, FALSE, FALSE, FALSE),
  pctscale = c(FALSE, FALSE, TRUE, TRUE)
)

for (i in 1:nrow(labels)) {
  metric <- labels[i, 'metric']
  label <- labels[i, 'label']
  plt <- QC_metric_VlnPlot(fig_objs, metric, label, caption)
  if (labels$logscale[i]) {
    plt <- plt + scale_y_log10()
  }
  if (labels$pctscale[i]) {
    plt <- plt + scale_y_continuous(labels = scales::percent)
  }
  filename <- paste0(c(metric, 'vln'), collapse='_')
  figures[[filename]] <- plt
  my_plot_save(image = plt, 
               path = here('figures/3p/qc_vln_plots_after_filtering/', paste0(filename, '.svg')), 
               width = 17, height = 7)
}
plt <- figures[['nFeature_RNA_vln']] + figures[['nCount_RNA_vln']] + figures[['rbRatio_vln']] + figures[['mtRatio_vln']] +
  plot_layout(nrow=2, ncol=2, guides = 'collect', axes = 'collect_x', axis_titles = 'collect_x')
my_plot_save(image = plt, 
             path = here('figures/3p/qc_vln_plots_after_filtering/QC_metrics_combined.svg'), 
             width = 24, height = 12)
