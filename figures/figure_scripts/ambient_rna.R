library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
rbindlist <- data.table::rbindlist

metadata <- read.csv(here('config/metadata.csv'))
objs <- readRDS(here('rds/01-soup_channels.rds'))

extractSoupXContamEst <- function(sc){
  rho <- sc$fit$rhoEst
  rho_low <- sc$fit$rhoFWHM[1]
  rho_high <- sc$fit$rhoFWHM[2]
  return(list(rho_low = rho_low, rho = rho, rho_high = rho_high))
}
soupx_results <- data.table::rbindlist(lapply(objs, extractSoupXContamEst), idcol = 'Sample')

soupx_results %>%
  merge(metadata, by = 'Sample') %>%
  ggplot(aes(x=paste0(Individual, Replicate), y=rho)) +
  geom_col() +
  geom_errorbar(aes(ymin=rho_low, ymax=rho_high), width=0.4) +
  facet_wrap(~ Kit, scales='free_x') +
  labs(x='Sample', y='Contamination fraction')

ggsave(here('figures/ambient_rna.png'), width = unit(5, 'in'), height = unit(5, 'in'))
