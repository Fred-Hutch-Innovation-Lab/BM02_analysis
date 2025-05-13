# Setup ----
source(here('figure_scripts/utils.R'))

# Load data ----
sat_curves <- read.csv(here('data/5p/sequencing_saturation.csv')) |>
  select(-Kit) |>
  mutate(across(all_of(c('median_genes', 'median_umi', 'nreads')), as.numeric))
models <- readRDS(here('rds/5p/sat_models.rds'))
models_u <- models$models_u
models_g <- models$models_g
rm(models)

# Prepare plotdata ----
sat_curves <- merge(sat_curves, metadata_5p, by='Sample', all=TRUE)

# Modeling ----

model_function <- function(nreads, a, c) {
  # michaelis menten
  (a * nreads) / (c + nreads)
}

## predict data ----
readrange <- seq(0, 25000, 1000)
predicted_g <- lapply(models_g, function(model) {
  predict(model$model, list(nreads=readrange))
}) |> as.data.table() |>
  mutate(nreads=readrange) |>
  melt(id.vars = 'nreads', variable.name = 'kit', value.name = 'median_genes')
predicted_u <- lapply(models_u, function(model) {
  predict(model$model, list(nreads=readrange))
}) |> as.data.table() |>
  mutate(nreads=readrange) |>
  melt(id.vars = 'nreads', variable.name = 'kit', value.name = 'median_umi')

## Extract RD50 ----
model_coef_g <- lapply(models_g, function(x) summary(x$model)$coefficients)
model_coef_u <- lapply(models_u, function(x) summary(x$model)$coefficients)

rd50_point_u <- lapply(names(model_coef_u), function(kit) {
  x <- model_coef_u[[kit]]['c', 'Estimate']
  y <- predict(models_u[[kit]]$model, list(nreads=x))
  return(list(Kit=kit, nreads=x, umi=y))
}) |> rbindlist()

rd50_point_g <- lapply(names(model_coef_g), function(kit) {
  x <- model_coef_g[[kit]]['c', 'Estimate']
  y <- predict(models_g[[kit]]$model, list(nreads=x))
  return(list(Kit=kit, nreads=x, genes=y))
}) |> rbindlist()
# 
rd50_ribbon_g <- lapply(names(model_coef_g), function(kit) {
  x <- model_coef_g[[kit]]['c', 'Estimate']
  se <- model_coef_g[[kit]]['c', 'Std. Error']
  x <- seq(x-se,x+se, length.out=10)
  y <- predict(models_g[[kit]]$model, list(nreads=x))
  return(list(Kit=kit, x=x, y=y))
}) |> rbindlist()

rd50_ribbon_u <- lapply(names(model_coef_g), function(kit) {
  x <- model_coef_u[[kit]]['c', 'Estimate']
  se <- model_coef_u[[kit]]['c', 'Std. Error']
  x <- seq(x-se,x+se, length.out=10)
  y <- predict(models_u[[kit]]$model, list(nreads=x))
  return(list(Kit=kit, x=x, y=y))
}) |> rbindlist()

## Plot ----
# ggplot() +
#   geom_line(data = predicted_g, aes(x=nreads, y=median_genes, color=kit)) +
#   # geom_point(data = filter(sat_curves, !is.na(median_genes)),
#   #            aes(x=nreads, y=median_genes, color = Kit), size=.5) +
#   geom_point(data = rd50_point_g, aes(x=nreads, y=genes, color=Kit), shape=108, size=9) +
#   geom_ribbon(data = rd50_ribbon_g, aes(ymin=y-100, ymax=y+100, x=x, color=Kit, fill=Kit), linetype='dotted', alpha=0.25, show.legend = FALSE) +
#   scale_color_manual(values = unlist(color_palette$kits), labels = label_function(mode='WT'), na.translate = FALSE) +
#   scale_fill_manual(values = unlist(color_palette$kits)) +
#   guides(color = guide_legend(override.aes = list(size = 2))) +
#   labs(x='Average reads per cell', y='Median genes per cell', color = 'Kit') 

ggplot(filter(sat_curves, !is.na(median_genes)),
       aes(x=nreads, y=median_genes, linetype = paste0(Individual, Replicate), color = Kit)) +
  geom_line() +
  lims(x=c(0,35000), y=c(0,3000)) +
  geom_point(data = rd50_point_g, aes(x=nreads, y=genes, color=Kit),
             shape=43, size=14,
             inherit.aes = FALSE, show.legend = FALSE) +
  scale_color_manual(values = unlist(color_palette$kits), labels = label_function(mode='clean')) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  labs(x='Average reads per cell', y='Median genes per cell', linetype='Sample', color = 'Kit') ->
  figures[['sat_curves_genes']]
figures[['sat_curves_genes']]
my_plot_save(image = figures[['sat_curves_genes']], 
             path = here('figures/5p/saturation_curves/saturation_curves_genes.svg'), 
             width = 7, height = 5)


ggplot(filter(sat_curves, !is.na(median_umi)),
       aes(x=nreads, y=median_umi, linetype = paste0(Individual, Replicate), color = Kit)) +
  geom_line() +
  geom_point(data = rd50_point_u, aes(x=nreads, y=umi, color=Kit),
             shape=43, size=14,
             inherit.aes = FALSE, show.legend = FALSE) +
  lims(x=c(0,35000), y=c(0,9000)) +
  scale_color_manual(values = unlist(color_palette$kits), labels = label_function(mode='clean'), na.translate = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  labs(x='Average reads per cell', y='Median transcripts per cell', linetype='Sample', color = 'Kit') ->
  figures[['sat_curves_umi']]
figures[['sat_curves_umi']]
my_plot_save(image = figures[['sat_curves_umi']], 
             path = here('figures/5p/saturation_curves/saturation_curves_umi.svg'), 
             width = 7, height = 5)

ggplot() +
  geom_line(data = predicted_g, aes(x=nreads, y=median_genes, color=kit)) +
  geom_point(data = filter(sat_curves, !is.na(median_genes)),
             aes(x=nreads, y=median_genes, color = Kit), size=.5) +
  geom_point(data = rd50_point_g, aes(x=nreads, y=genes, color=Kit), shape=108, size=8) +
  geom_ribbon(data = rd50_ribbon_g, aes(ymin=y-100, ymax=y+100, x=x, color=Kit, fill=Kit), linetype='dotted', alpha=0.25, show.legend = FALSE) +
  scale_color_manual(values = unlist(color_palette$kits), labels = label_function(mode='clean'), na.translate = FALSE) +
  scale_fill_manual(values = unlist(color_palette$kits)) +
  lims(x=c(0,35000), y=c(0,3000)) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  labs(x='Average reads per cell', y='Median genes per cell', color = 'Kit', 
       caption = 'Points indicate raw data, lines show modeled fit.
       Vertical lines indicate RD50.
       Shaded region indicates 95% confidence estimates.') ->
  figures[['sat_curves_model_gene']]
figures[['sat_curves_model_gene']]
my_plot_save(image = figures[['sat_curves_model_gene']], 
             path = here('figures/5p/saturation_curves/saturation_curves_modeled_gene.svg'), 
             width = 7, height = 5)

ggplot() +
  geom_line(data = predicted_u, aes(x=nreads, y=median_umi, color=kit)) +
  geom_point(data = filter(sat_curves, !is.na(median_umi)),
             aes(x=nreads, y=median_umi, color = Kit), size=.5) +
  geom_point(data = rd50_point_u, aes(x=nreads, y=umi, color=Kit), shape=108, size=8) +
  geom_ribbon(data = rd50_ribbon_u, aes(ymin=y-350, ymax=y+350, x=x, color=Kit, fill=Kit), linetype='dotted', alpha=0.25, show.legend = FALSE) +
  scale_color_manual(values = unlist(color_palette$kits), labels = label_function(mode='WT'), na.translate = FALSE) +
  scale_fill_manual(values = unlist(color_palette$kits)) +
  lims(x=c(0,35000), y=c(0,9000)) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  labs(x='Average reads per cell', y='Median transcripts per cell', color = 'Kit', 
       caption = 'Points indicate raw data, lines show modeled fit.
       Vertical lines indicate RD50.
       Shaded region indicates 95% confidence estimates.') ->
  figures[['sat_curves_model_umi']]
figures[['sat_curves_model_umi']]
my_plot_save(image = figures[['sat_curves_model_umi']], 
             path = here('figures/5p/saturation_curves/saturation_curves_modeled_umi.svg'), 
             width = 7, height = 5)

write_plot_data(sat_curves, here('figure_data/5p/sequencing_saturation/sat_curves_5p.txt'))

# rd50_g <- lapply(names(model_coef_g), function(kit) {
#   x <- model_coef_g[[kit]]['c', 'Estimate']
#   se <- model_coef_g[[kit]]['c', 'Std. Error']
#   return(list(x=x, se=se))
# }) 
# names(rd50_g) <- names(model_coef_g)
# rd50_g <- rbindlist(rd50_g, idcol='kit')
# write_plot_data(rd50_g, here('figure_data/5p/sequencing_saturation/rd50_g.txt'))
# 
# rd50_u <- lapply(names(model_coef_u), function(kit) {
#   x <- model_coef_u[[kit]]['c', 'Estimate']
#   se <- model_coef_u[[kit]]['c', 'Std. Error']
#   return(list(x=x, se=se))
# }) 
# names(rd50_u) <- names(model_coef_u)
# rd50_u <- rbindlist(rd50_u, idcol='kit')
# write_plot_data(rd50_u, here('figure_data/5p/sequencing_saturation/rd50_u.txt'))
