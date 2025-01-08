# Setup ----
source(here('figure_scripts/utils.R'))

# Load data ----
sat_curves <- read.csv(here('data/3p/saturation_curves/manually_downsampled.csv')) |>
  select(-Kit) |>
  mutate(across(all_of(c('median_genes', 'median_umi', 'nreads')), as.numeric))

# Prepare plotdata ----
zeroes <- data.frame(Sample=metadata_3p$Sample, nreads=0, median_genes=0, median_umi=0)
sat_curves <- rbind(sat_curves, zeroes)
sat_curves <- merge(sat_curves, metadata_3p, by='Sample', all=TRUE)

# Plot ----
ggplot(filter(sat_curves, !is.na(median_genes)),
       aes(x=nreads, y=median_genes, linetype = paste0(Individual, Replicate), color = Kit)) +
  geom_line() +
  scale_color_manual(values = unlist(color_palette$kits), labels = label_function) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  labs(x='Average reads per cell', y='Median genes per cell', linetype='Sample', color = 'Kit') ->
  figures[['sat_curves_genes']]
my_plot_save(image = figures[['sat_curves_genes']], 
             path = here('figures/3p/saturation_curves/saturation_curves_genes.svg'), 
             width = 7, height = 5)

ggplot(filter(sat_curves, !is.na(median_umi)),
       aes(x=nreads, y=median_umi, linetype = paste0(Individual, Replicate), color = Kit)) +
  geom_line() +
  scale_color_manual(values = unlist(color_palette$kits), labels = label_function, na.translate = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  labs(x='Average reads per cell', y='Median transcripts per cell', linetype='Sample', color = 'Kit') ->
  figures[['sat_curves_umi']]
my_plot_save(image = figures[['sat_curves_umi']], 
             path = here('figures/3p/saturation_curves/saturation_curves_umi.svg'), 
             width = 7, height = 5)
