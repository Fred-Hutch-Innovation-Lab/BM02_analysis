library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)       ## parses config yaml
library(ggrastr)    ## rasterized point layers for reduced image size
library(Seurat)
library(reshape2)  ## melt
library(patchwork)

## wrapped in function to avoid filling environment if sourced
main <- function(figures_list, filtering_stage = 'after', outdir='qc_vln_plots') {
  ## returns color_palette obj
  source(here('config/color_palette.R'))
  
  QC_metric_VlnPlot <- function(objs, metric, ylab = 'metric_value', caption=NA) {
    plotdata <- lapply(objs, function(x){unname(x@meta.data[[metric]])}) %>%
      melt(idcol='Sample') %>%
      merge(metadata, by.x='L1', by.y='Sample') 
    plotdata2 <- group_by(plotdata, Kit) %>% 
      summarize(med = median(value)) 
    
    ggplot(plotdata, aes(y=value, x=paste0(Individual, Replicate), fill = Kit)) +
      rasterize(geom_jitter(alpha=0.1, size=.1), dpi=300) +
      geom_violin(draw_quantiles = 0.5) +
      scale_fill_manual(values = color_palette$kits) +
      theme_bw() +
      theme(text = element_text(size = 20),
            panel.spacing=unit(0, "lines"),
            axis.text.x = element_text(angle=45, vjust=0.5)) + 
      facet_wrap(~ Kit, nrow = 1, scales='free_x') +
      labs(x='Sample', y = ylab, caption=caption) #+
    # coord_flip()
  }
  
  if (filtering_stage == 'after') {
    fig_objs <- readRDS(here('rds/3p/02-objs_post_cell_filtering.rds'))
    caption <- 'QC metrics per cell on data after filtering by global thresholds'
    caption <- NULL
  } else if (filtering_stage == 'before') {
    fig_objs <- readRDS(here('rds/3p/01_raw_objs.rds'))
    caption <- 'QC metrics per cell on data output from pipelines without additional filtering'
  }
  
  kit_order <- read.table(here('config/3p/kit_order.txt'))$V1
  metadata <- read.csv(here('config/3p/metadata.csv')) %>%
    mutate(Kit = factor(Kit, levels = kit_order))
  
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
    figures_list[[filename]] <- plt
    ggsave(plt, filename = paste0(filename, '.png'), path = here(file.path('figures/3p', outdir)),
           device = 'png', 
           height = unit(7, 'in'), width = unit(17, 'in'))
  }
  # legend <- get_legend(figures_list[[1]] + )
  # plot_grid(figures_list[[1]] + theme(legend.position = 'none'),
  #           figures_list[[2]] + theme(legend.position = 'none'),
  #           figures_list[[3]] + theme(legend.position = 'none'),
  #           figures_list[[4]] + theme(legend.position = 'none'),
  #           nrow=2, ncol=2, labels=c('A', 'B', 'C', 'D'))
  plt <- figures_list[[1]] + figures_list[[2]] + figures_list[[3]] + figures_list[[4]] +
    plot_layout(nrow=2, ncol=2, guides = 'collect', axes = 'collect_x', axis_titles = 'collect_x')
  ggsave(plt, filename = 'QC_metrics_combined.png', path = here(file.path('figures/3p', outdir)),
         device = 'png', 
         height = unit(12, 'in'), width = unit(23, 'in'))
  return(figures_list)
}

if (!exists(quote(figures))) {
  figures <- list()
} else {
  message('Figures object exists in environment, adding figures to existing obj. Be careful of overwrites.')
}

figures <- main(figures, filtering_stage = 'after', outdir='qc_vln_plots_after_filtering')

rm(main)
