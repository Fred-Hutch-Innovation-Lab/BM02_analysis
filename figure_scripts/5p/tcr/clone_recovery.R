# Setup ----
source(here('figure_scripts/utils.R'))
library(immunarch)  ## repertoire analyses
library(ggbeeswarm) ## swarm plots

# Load data ----
fig_objs <- readRDS(here('rds/5p/tcr/tcr_paired.Rds'))
metric_summaries <- read.csv(here('rds/5p/tcr/metric_summaries.csv')) |>
  mutate(Metric.Name = case_when(
    Metric.Name == 'TRA & TRB' ~ 'TRA+TRB',
    .default = Metric.Name
  )) |>
  mutate(Metric.Name = factor(Metric.Name, levels = c(
    'Total cells', 'TRA+TRB', 'TRA', 'TRB'
  )))

# Plot ----
## Reported chain recovery ----
plotdata <- merge(metric_summaries, metadata_5p, by='Sample')
ggplot(filter(plotdata, Metric.Name != 'Total cells'), aes(x=Metric.Name, y=pct)) +
  geom_boxplot() +
  geom_beeswarm(aes(shape=Individual)) +
  facet_wrap(~ Kit, labeller = labeller(Kit = label_function)) +
  lims(y=c(0,100)) +
  labs(x='Chain(s)', y='Recovery (% of T cells)') ->
  figures[['reported chain recovery']]
figures[['reported chain recovery']]
my_plot_save(image = figures[['reported chain recovery']], 
             path = here('figures/5p/tcr/reported_chain_recovery.svg'), 
             width = 7, height = 3.5)
write_plot_data(plotdata, file=here('figure_data/5p/tcr/reported_chain_recovery.txt'))

return_total_cells <- function(sample){
  metric_summaries$Metric.Value[
    metric_summaries$Sample %in% sample &
      metric_summaries$Metric.Name == 'Total cells'
  ] 
}

plotdata <- data.table::rbindlist(lapply(fig_objs, function(x) {
  x$data |>
    group_by(chain) |>
    summarize(n=sum(Clones))
}), idcol = 'Sample') |>
  # merge(metadata_5p, by='Sample') |>
  # mutate(chain = factor(chain, levels = c(
  #   'TRA;TRB',
  #   'TRA',
  #   'TRB',
  #   'TRA;TRA;TRB',
  #   'TRA;TRB;TRB',
  #   'TRA;TRA',
  #   'TRB;TRB',
  #   'TRA;TRA;TRB;TRB'
  # ))) |>
  mutate(chaingroup = case_when(
    grepl('TRA;TRB', chain) ~ 'TRA+TRB',
    grepl('TRA', chain) ~ 'TRA',
    grepl('TRB', chain) ~ 'TRB',
  )) |> 
  group_by(Sample, chaingroup) |>
    summarize(n=sum(n)) |>
  as.data.table() |> 
  dcast(Sample ~ chaingroup) |>
  mutate(TRA = TRA + `TRA+TRB`,
         TRB = TRB + `TRA+TRB`) |>
  melt() |>
  mutate(prop = value / return_total_cells(Sample),
         variable = factor(variable, c('TRA+TRB', 'TRA', 'TRB'))) |>
  merge(metadata_5p, by='Sample')
# plotdata$total_cells <- unlist(lapply(plotdata$Sample, return_total_cells))
ggplot(plotdata, aes(x=variable, y=prop)) +
  geom_boxplot() +
  geom_point(aes(shape = Individual)) +
  scale_y_continuous(labels=scales::percent, limits = c(0,1)) +
  facet_wrap(~ Kit, scales='free_x', labeller = labeller(Kit = label_function)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(x='Chain(s)', y='% of cells', caption='% of all detected T cells regardless of chain recovery') ->
  figures[['chain recovery % all T']]
figures[['chain recovery % all T']]
my_plot_save(image = figures[['chain recovery % all T']], 
             path = here('figures/5p/tcr/all_chain_recovery.svg'), 
             width = 8, height = 3.5)
write_plot_data(plotdata, file=here('figure_data/5p/tcr/all_chain_recovery.txt'))

# Simpson index ----

plotdata <- lapply(fig_objs, function(x) repDiversity(x$data, 'inv.simp')) |>
  as.data.table() |>
  melt() |>
  merge(metadata_5p, by.x='variable', by.y='Sample') 
ggplot(plotdata, aes(x=Individual, y=value)) +
  geom_col() +
  facet_grid(~ Kit, scales='free_x', space='free_x', labeller = labeller(Kit = label_function)) +
  labs(x='Individual', y='Inverse Simpson index') ->
  figures[['simpson index']]
figures[['simpson index']]
my_plot_save(image = figures[['simpson index']], 
             path = here('figures/5p/tcr/simpson_index.svg'), 
             width = 5, height = 3.5)
write_plot_data(plotdata, file=here('figure_data/5p/tcr/simpson_index.txt'))