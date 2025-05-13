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

return_total_cells <- function(sample){
  metric_summaries$Metric.Value[
    metric_summaries$Sample %in% sample &
      metric_summaries$Metric.Name == 'Total cells'
  ] 
}

# Plot ----
# ## Reported chain recovery ---
# plotdata <- merge(metric_summaries, metadata_5p, by='Sample')
# ggplot(filter(plotdata, Metric.Name != 'Total cells'), aes(x=Metric.Name, y=pct)) +
#   geom_boxplot() +
#   geom_beeswarm(aes(shape=Individual)) +
#   facet_wrap(~ Kit, labeller = labeller(Kit = label_function(mode='TCR'))) +
#   lims(y=c(0,100)) +
#   labs(x='Chain(s)', y='Recovery (% of T cells)') ->
#   figures[['reported chain recovery']]
# figures[['reported chain recovery']]
# my_plot_save(image = figures[['reported chain recovery']], 
#              path = here('figures/5p/tcr/reported_chain_recovery.svg'), 
#              width = 7, height = 3.5)
# write_plot_data(plotdata, file=here('figure_data/5p/tcr/reported_chain_recovery.txt'))
# 
# plotdata <- data.table::rbindlist(lapply(fig_objs, function(x) {
#   x$data |>
#     group_by(chain) |>
#     summarize(n=sum(Clones))
# }), idcol = 'Sample') |>
#   # merge(metadata_5p, by='Sample') |>
#   # mutate(chain = factor(chain, levels = c(
#   #   'TRA;TRB',
#   #   'TRA',
#   #   'TRB',
#   #   'TRA;TRA;TRB',
#   #   'TRA;TRB;TRB',
#   #   'TRA;TRA',
#   #   'TRB;TRB',
#   #   'TRA;TRA;TRB;TRB'
#   # ))) |>
#   mutate(chaingroup = case_when(
#     grepl('TRA;TRB', chain) ~ 'TRA+TRB',
#     grepl('TRA', chain) ~ 'TRA',
#     grepl('TRB', chain) ~ 'TRB',
#   )) |> 
#   group_by(Sample, chaingroup) |>
#     summarize(n=sum(n)) |>
#   as.data.table() |> 
#   dcast(Sample ~ chaingroup) |>
#   mutate(TRA = TRA + `TRA+TRB`,
#          TRB = TRB + `TRA+TRB`) |>
#   melt() |>
#   mutate(prop = value / return_total_cells(Sample),
#          variable = factor(variable, c('TRA+TRB', 'TRA', 'TRB'))) |>
#   merge(metadata_5p, by='Sample')
# # plotdata$total_cells <- unlist(lapply(plotdata$Sample, return_total_cells))
# ggplot(plotdata, aes(x=variable, y=prop)) +
#   geom_boxplot() +
#   geom_point(aes(shape = Individual)) +
#   scale_y_continuous(labels=scales::percent, limits = c(0,1)) +
#   facet_wrap(~ Kit, scales='free_x', labeller = labeller(Kit = label_function(mode='TCR'))) +
#   theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
#   labs(x='Chain(s)', y='% of T cells') ->
#   figures[['chain recovery % all T']]
# figures[['chain recovery % all T']]
# my_plot_save(image = figures[['chain recovery % all T']], 
#              path = here('figures/5p/tcr/all_chain_recovery.svg'), 
#              width = 8, height = 3.5)
# write_plot_data(plotdata, file=here('figure_data/5p/tcr/all_chain_recovery.txt'))

## Simpson index ----

plotdata <- lapply(fig_objs, function(x) repDiversity(x$data |> filter(chain == 'TRA;TRB'), 'inv.simp')) |>
  as.data.table() |>
  melt() |>
  merge(metadata_5p, by.x='variable', by.y='Sample') 
ggplot(plotdata, aes(x=Individual, y=value, fill=Kit)) +
  geom_col() +
  scale_fill_manual(values = color_palette$kits, labels = label_function(mode='clean'), breaks = kit_order_5p) +
  facet_grid(~ Kit, scales='free_x', space='free_x',
             labeller = labeller(Kit = label_function(mode='clean', wrap = TRUE, width = 10))) +
  labs(x='Individual', y='Inverse Simpson index', caption = 'TRA;TRB clones only') ->
  figures[['simpson index']]
figures[['simpson index']]
my_plot_save(image = figures[['simpson index']], 
             path = here('figures/5p/tcr/simpson_index.svg'), 
             width = 9, height = 4)
write_plot_data(plotdata, file=here('figure_data/5p/tcr/simpson_index.txt'))

## recovery proportion stacked ----
plotdata <- data.table::rbindlist(lapply(fig_objs, function(x) {
  x$data |>
    group_by(chain) |>
    summarize(n=n())
}), idcol = 'Sample') |>
  merge(metadata_5p, by='Sample') |>
  group_by(Sample)
plotdata$prop <- plotdata$n / unlist(lapply(plotdata$Sample, return_total_cells))

plotdata <- plotdata |> 
  # group_by(Sample, File, Kit, Individual, Replicate) |>
  # summarize(prop = 1-sum(prop)) |>
  # mutate(chain = 'None') |>
  # rbind(plotdata) |>
  mutate(chain = factor(chain, levels = c(
    'TRA;TRB',
    'TRA',
    'TRB',
    'TRA;TRA;TRB',
    'TRA;TRB;TRB',
    'TRA;TRA',
    'TRB;TRB',
    'TRA;TRA;TRB;TRB',
    'None'
  ))) |>
  arrange(Sample, chain)
ggplot(plotdata, aes(x=Individual, y=prop, group = Individual, fill = chain)) +
  geom_col() +
  # geom_point(aes(shape = Individual)) +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_brewer(type = 'qual') +
  facet_wrap(~ Kit, scales='free_x', labeller = label_function(mode='clean', wrap=TRUE, width=10)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(x='Chain(s)', y='% of T cells') ->
  figures[['chain recovery stacked']]
figures[['chain recovery stacked']]
my_plot_save(image = figures[['chain recovery stacked']], 
             path = here('figures/5p/tcr/expanded_clone_recovery_proportions.svg'), 
             width = 7, height = 4)
write_plot_data(plotdata, file=here('figure_data/5p/tcr/expanded_clone_recovery_proportions.txt'))

## Richness ----
plotdata <- data.table::rbindlist(lapply(fig_objs, function(x) {
  x$data |>
    group_by(chain) |>
    summarize(n=n())
}), idcol = 'Sample') |>
  merge(metadata_5p, by='Sample') |>
  group_by(Sample) |>
  mutate(chain = factor(chain, levels = c(
    'TRA;TRB',
    'TRA',
    'TRB',
    'TRA;TRA;TRB',
    'TRA;TRB;TRB',
    'TRA;TRA',
    'TRB;TRB',
    'TRA;TRA;TRB;TRB'
  ))) |>
  arrange(Sample, desc(chain))
plotdata$total_cells <- unlist(lapply(plotdata$Sample, return_total_cells))

ggplot(plotdata, aes(x=Individual)) +
  geom_col(aes(y=n, fill=chain), position = position_stack(reverse = TRUE)) + 
  geom_point(aes(y=total_cells), shape=95, size=15, color = 'darkblue') +
  scale_fill_brewer(type = 'qual') +
  facet_grid(~ Kit, scales='free_x',space = 'free_x' , labeller = label_function(mode='clean', wrap=TRUE, width=10)) + 
  labs(x= 'Sample', y='Unique clonotypes', caption='Blue line shows total T cells recovered') ->
  figures[['unique clone recovery count']]
figures[['unique clone recovery count']]

my_plot_save(image = figures[['unique clone recovery count']], 
             path = here('figures/5p/tcr/expanded_clone_recovery_counts.svg'), 
             width = 7, height = 4)
write_plot_data(plotdata, file=here('figure_data/5p/tcr/expanded_clone_recovery_counts.txt'))
