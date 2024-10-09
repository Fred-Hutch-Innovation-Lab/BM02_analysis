#!/usr/bin/env Rscript

library(argparser)  ## Command line arguments
library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
rbindlist <- data.table::rbindlist

main <- function() {
  argv <- arg_parser('Soupx barchart') %>%
    add_argument('--outfile', help='output file', default='ambient_rna.png') %>%
    add_argument('--metadata', help='metadata file', default='config/metadata.csv') %>%
    add_argument('--input', help='input file', default='rds/01-soup_channels.rds') %>%
    add_argument('--save_plot', help='Overwrite saved file?', default=FALSE) %>%
    parse_args()
  
  metadata <- read.csv(here(file.path(argv$metadata)))
  objs <- readRDS(here(file.path(argv$input)))
  
  extractSoupXContamEst <- function(sc){
    rho <- sc$fit$rhoEst
    rho_low <- sc$fit$rhoFWHM[1]
    rho_high <- sc$fit$rhoFWHM[2]
    return(list(rho_low = rho_low, rho = rho, rho_high = rho_high))
  }
  
  soupx_results <- data.table::rbindlist(lapply(objs, extractSoupXContamEst), idcol = 'Sample')
  soupx_plot <- soupx_results %>%
    merge(metadata, by = 'Sample') %>%
    ggplot(aes(x=paste0(Individual, Replicate), y=rho)) +
    geom_col() +
    geom_errorbar(aes(ymin=rho_low, ymax=rho_high), width=0.4) +
    facet_wrap(~ Kit, scales='free_x') +
    labs(x='Sample', y='Contamination fraction')
  
  if (argv$save_plot) {
    ggsave(soupx_plot, here(file.path(argv$outfile)), width = unit(5, 'in'), height = unit(5, 'in'))
  }
  return(soupx_plot)
}

soupx_plot <- main()
