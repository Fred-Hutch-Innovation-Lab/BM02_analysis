library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)       ## parses config yaml
library(ggrastr)    ## rasterized point layers for reduced image size
library(Seurat)
library(reshape2)  ## melt
library(patchwork)

source(here('config/color_palette.R'))
source(here('config/kit_order.R'))

fig_objs <- readRDS(here('rds/5p/wt/02-objs_post_cell_filtering.rds'))
metadata <- read.csv(here('config/5p/metadata.csv')) %>%
  mutate(Kit = factor(Kit, levels = kit_order_5p))
figures_list <- list()
  
QC_metric_VlnPlot <- function(objs, metric, ylab = 'metric_value', caption=NA) {
  plotdata <- lapply(objs, function(x){unname(x@meta.data[[metric]])}) %>%
    melt(idcol='Sample') %>%
    merge(metadata, by.x='L1', by.y='Sample') 
  plotdata2 <- group_by(plotdata, Kit) %>% 
    summarize(med = median(value)) 
  
  ggplot(plotdata, aes(y=value, x=Individual, fill = Kit)) +
    rasterize(geom_jitter(alpha=0.1, size=.1), dpi=100) +
    geom_violin(draw_quantiles = 0.5) +
    scale_fill_manual(values = color_palette$kits, labels = label_function) +
    theme_bw() +
    theme(text = element_text(size = 20),
          panel.spacing=unit(0, "lines"),
          axis.text.x = element_text(angle=45, vjust=0.5)) + 
    facet_grid(~ Kit, scales='free_x', space='free_x', labeller = labeller(Kit = label_function)) +
    labs(x='Sample', y = ylab)
}

labels <- data.frame(
  metric = c('nFeature_RNA', 'nCount_RNA', 'rbRatio', 'mtRatio'),
  label  = c('Feature count', 'UMI count', 'Ribosomal read %', 'Mitochondrial read %'),
  logscale = c(FALSE, FALSE, FALSE, FALSE),
  pctscale = c(FALSE, FALSE, TRUE, TRUE)
)

for (i in 1:nrow(labels)) {
  metric <- labels[i, 'metric']
  label <- labels[i, 'label']
  plt <- QC_metric_VlnPlot(fig_objs, metric, label)
  if (labels$logscale[i]) {
    plt <- plt + scale_y_log10()
  }
  if (labels$pctscale[i]) {
    plt <- plt + scale_y_continuous(labels = scales::percent)
  }
  filename <- paste0(c(metric, 'vln'), collapse='_')
  figures_list[[filename]] <- plt
  ggsave(plt, filename = paste0(filename, '.png'), path = here(file.path('figures/5p/wt', 'qc_vln_plots_after_filtering')),
         device = 'png', 
         height = unit(7, 'in'), width = unit(12, 'in'))
}
plt <- figures_list[[1]] + figures_list[[2]] + figures_list[[3]] + figures_list[[4]] +
  plot_layout(nrow=2, ncol=2, guides = 'collect', axes = 'collect_x', axis_titles = 'collect_x')
ggsave(plt, filename = 'QC_metrics_combined.png', path = here(file.path('figures/5p/wt', 'qc_vln_plots_after_filtering')),
       device = 'png', 
       height = unit(12, 'in'), width = unit(18, 'in'))
