# Setup ----
source(here('figure_scripts/utils.R'))
# library(ggh4x)

# Load data ----
fig_objs <- readRDS(here('rds/3p/01-soup_channels.rds'))

# Prepare plotdata ----
extractSoupXContamEst <- function(sc){
  rho <- sc$fit$rhoEst
  rho_low <- sc$fit$rhoFWHM[1]
  rho_high <- sc$fit$rhoFWHM[2]
  return(list(rho_low = rho_low, rho = rho, rho_high = rho_high))
}
plotdata <- data.table::rbindlist(lapply(fig_objs, extractSoupXContamEst), idcol = 'Sample') %>%
  merge(metadata_3p, by = 'Sample', all.y=TRUE)

# Plot ----
ggplot(plotdata, aes(x=paste0(Individual, Replicate), y=rho)) +
  geom_col() +
  geom_errorbar(aes(ymin=rho_low, ymax=rho_high), width=0.4) +
  facet_grid(~ Kit, scales='free_x',
             labeller = labeller(Kit = label_function())
               ) +
  scale_y_continuous(labels = scales::percent) +
  labs(x='Sample', y='Ambient RNA') ->
figures[['soupx']]
my_plot_save(image = figures[['soupx']], 
             path = here('figures/3p/ambient_rna.svg'), 
             width = 12, height = 3.5)
write_plot_data(plotdata, here('figure_data/3p/ambient_rna.txt'))
