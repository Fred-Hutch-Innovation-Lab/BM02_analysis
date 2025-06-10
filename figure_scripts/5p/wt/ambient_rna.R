# Setup ----
source(here('figure_scripts/utils.R'))
# library(ggh4x)

# Load data ----
fig_objs <- readRDS(here('rds/5p/wt/01-soup_channels.rds'))
fig_objs$NextGEM5P_F5 <- NULL
fig_objs$NextGEM5P_F4 <- NULL
capitalize_first <- function(s) {
  unname(sapply(s, function(str) {
    if (nchar(str) == 0) return(str)  # Handle empty strings
    paste0(toupper(substr(str, 1, 1)), substr(str, 2, nchar(str)))
  }))
}
names(fig_objs)[grepl('parse', names(fig_objs), ignore.case=TRUE)] <-
  capitalize_first(names(fig_objs)[grepl('parse', names(fig_objs), ignore.case=TRUE)])
# Prepare plotdata ----
extractSoupXContamEst <- function(sc){
  rho <- sc$fit$rhoEst
  rho_low <- sc$fit$rhoFWHM[1]
  rho_high <- sc$fit$rhoFWHM[2]
  return(list(rho_low = rho_low, rho = rho, rho_high = rho_high))
}
plotdata <- data.table::rbindlist(lapply(fig_objs, extractSoupXContamEst), idcol = 'Sample') %>%
  merge(metadata_5p, by = 'Sample')

# Plot ----
ggplot(plotdata, aes(x=paste0(Individual, Replicate), y=rho)) +
  geom_col() +
  geom_errorbar(aes(ymin=rho_low, ymax=rho_high), width=0.4) +
  facet_grid(~ Kit, scales='free_x', space = 'free_x',
             labeller = labeller(Kit = label_function(mode='clean'))
  ) +
  scale_y_continuous(labels = scales::percent) +
  labs(x='Sample', y='Ambient RNA') ->
  figures[['soupx']]
my_plot_save(image = figures[['soupx']], 
             path = here('figures/5p/wt/ambient_rna.svg'), 
             width = 8, height = 3.5)
write_plot_data(plotdata, here('figure_data/5p/wt/ambient_rna.txt'))
