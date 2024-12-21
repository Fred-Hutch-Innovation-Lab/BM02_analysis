#!/usr/bin/env Rscript

# library(argparser)  ## Command line arguments
library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(ggh4x)
rbindlist <- data.table::rbindlist

main <- function() {
  # argv <- arg_parser('Soupx barchart') %>%
  #   add_argument('--outfile', help='output file', default='ambient_rna.pdf') %>%
  #   add_argument('--metadata', help='metadata file', default='config/metadata.csv') %>%
  #   add_argument('--input', help='input file', default='rds/01-soup_channels.rds') %>%
  #   add_argument('--save_plot', help='Overwrite saved file?', default=FALSE) %>%
  #   parse_args()
  # 
  kit_order <- read.table(here('config/kit_order.txt'))$V1
  metadata <- read.csv(here('config/3p/metadata.csv')) %>%
    mutate(Kit = factor(Kit, levels = kit_order))
  
  objs <- readRDS(here('rds/3p/01-soup_channels.rds'))
  
  extractSoupXContamEst <- function(sc){
    rho <- sc$fit$rhoEst
    rho_low <- sc$fit$rhoFWHM[1]
    rho_high <- sc$fit$rhoFWHM[2]
    return(list(rho_low = rho_low, rho = rho, rho_high = rho_high))
  }
  
  soupx_results <- data.table::rbindlist(lapply(objs, extractSoupXContamEst), idcol = 'Sample')
  soupx_plot <- soupx_results %>%
    merge(metadata, by = 'Sample', all.y=TRUE) %>%
    ggplot(aes(x=paste0(Individual, Replicate), y=rho)) +
    geom_col() +
    geom_errorbar(aes(ymin=rho_low, ymax=rho_high), width=0.4) +
    facet_manual(~ Kit, scales='free_x',
                 design="#AABBCC#
                 DDEEFFGG"
                 # design=rbind(c(1,2,3, NA), c(4,5,6,7))
                 ) +
    # theme(axis.text.x = element_text(angle=45, vjust=1)) +
    scale_y_continuous(labels = scales::percent) +
    labs(x='Sample', y='Contamination estimate', caption = 'Range bars indicate upper and lower bounds of contamination estimate')
  
  # if (argv$save_plot) {
    ggsave(plot=soupx_plot, here('figures/3p/ambient_rna.png'), width = unit(7, 'in'), height = unit(5, 'in'))
  # }
  return(soupx_plot)
}

soupx_plot <- main()
