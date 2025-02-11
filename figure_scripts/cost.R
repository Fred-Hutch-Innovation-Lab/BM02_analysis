# Setup ----
source(here('figure_scripts/utils.R')) 
# library(readxl)


# Load data ----
## Kit costs
cost_data_ref <- read.table(here('data/cost.txt'), sep='\t', header = TRUE)  |>
  mutate(cost = readr::parse_number(cost),
         max_cells = readr::parse_number(max_cells)) |>
  
  group_by(kit_family) |>
  arrange(kit_family, max_samples, max_cells) |>
  mutate(kit_throughput_n = as.factor(row_number())) |>
  mutate(baseline_cost_per_cell = cost / max_cells)
  
## Models ----
models <- readRDS(here('rds/3p/sat_models.rds'))
models_u <- models$models_u
models_g <- models$models_g
models <- readRDS(here('rds/5p/sat_models.rds'))
models_u <- c(models_u, models$models_u)
models_g <- c(models_g, models$models_g)
rm(models)
parse_coef <- function(x){
  summary(x$model)$coefficients |>
    as.data.frame() |>
    rownames_to_column('var') |>
    select(var, Estimate)
}
model_coef_g <- lapply(models_g, parse_coef) |>
  rbindlist(idcol = 'Kit') |> 
  dcast(Kit ~ var)
model_coef_u <- lapply(models_u, parse_coef) |>
  rbindlist(idcol = 'Kit') |>
  dcast(Kit ~ var)



## Conversion efficiency
read_eff <- rbind(
  read.csv(here('data/3p/sequencing_efficiency_2.txt'), sep='\t') |>
    mutate(prop = reads_in_matrix / reads_in_fastqs) |>
    merge(metadata_3p, by='Sample') |>
    group_by(Kit) |>
    summarize(eff = mean(prop)),
  read.csv(here('data/5p/sequencing_efficiency.txt'), sep='\t') |>
    mutate(prop = reads_in_matrix / reads_in_fastqs) |>
    merge(metadata_5p, by='Sample') |>
    group_by(Kit) |>
    summarize(eff = mean(prop)))
read_eff

# Prepare plotdata ----

## RD50
# model_function <- function(nreads, a, c) {
#   # michaelis menten
#   (a * nreads) / (c + nreads)
# }
# y = (a * nreads) / (c + nreads)
# y * (c + nreads) = (a * nreads)
# y*c + y*nreads = (a * nreads)
# y*c = (a * nreads) - y*nreads
# y*c = nreads*(a - y)
# y*c / (a - y) = nreads

readfit <- function(y, a, c) {
  y*c / (a - y)
}

# req_reads_g <- lapply(model_coef_g, function(x) {
#   readfit(1000, x['a',1], x['c',1])
# })|>
#   as.data.frame() |>
#   melt(value.name = 'target_seq_depth')
# req_reads_g
# 
# req_reads_u <- lapply(model_coef_u, function(x) {
#   readfit(2000, x['a',1], x['c',1])
# }) |>
#   as.data.frame() |>
#   melt(value.name = 'target_seq_depth')
# req_reads_u

## Define vars ----

ncells = 20000
nsamples = 4
seq_cost_coef = 1.5e-06 # 1,000,000,000 reads / $1,500

# Plot ----

## Sequencing

## Table ----

plotdata <- read_eff |>
  dplyr::rename('read_conversion_eff' = eff) |>
  # rowwise() |>
  # mutate(max_umi = model_coef_u[[Kit]]['a',1],
  #        rd50 = model_coef_u[[Kit]]['c',1]) |>
  merge(model_coef_u, by='Kit') |>
  dplyr::rename('rd50' = c, 'max_umi' = a) |>
  mutate(req_reads_2000_umi = readfit(2000, max_umi, rd50)) |>
  ungroup() |>
  mutate(target_reads_per_cell = req_reads_2000_umi / read_conversion_eff) |>
  mutate(seq_cost = target_reads_per_cell * ncells * nsamples * seq_cost_coef)  |>
  mutate(Kit = factor(Kit, levels = c(
    'Flex', 'NextGEM5P', 'NextGEM3P', 'GEMX5P', 'GEMX3P',
    'Fluent_v4', 'Fluent_V', 'Parse_v2', 'Parse_v3', 'Scale'
  ))) |>
  arrange(Kit) |>
  select(Kit, max_umi, rd50, req_reads_2000_umi, read_conversion_eff, target_reads_per_cell, seq_cost) 
write_plot_data(plotdata, file = here('figure_data/seq_cost_table.txt'))
gt(plotdata) |>
  cols_label(Kit = 'Kit',
             max_umi = 'Max UMI',
             rd50 = 'rd50',
             req_reads_2000_umi = html('Mapped reads<br>needed for<br>2000 UMI'),
             read_conversion_eff = html('Read<br>conversion<br>efficiency'),
             target_reads_per_cell = html('Target<br>seq depth<br>per cell'),
             seq_cost = html('Example<br>sequencing<br>cost')
  ) |>
  fmt_currency(c('seq_cost'), decimals = 0) |>
  fmt_number(c('read_conversion_eff'), decimals = 2) |>
  fmt_number(c('max_umi', 'rd50', 'target_reads_per_cell', 'req_reads_2000_umi'), decimals = 0) |>
  tab_footnote(
    footnote = paste0('Projected cost based on estimates for ', ncells, ' cells from ', nsamples, ' samples with a median UMI recovery of ', 2000, ' and an estimated sequencing cost of $', seq_cost_coef, ' per read'),
    locations = cells_column_labels(columns = seq_cost)
  ) ->
  figures[['seq_cost_table']]
figures[['seq_cost_table']]
gt::gtsave(figures[['seq_cost_table']], here('figures/cost/seq_cost.html'))


### UMI ----
seq_cost_modeling <- expand.grid(
  model_coef_u$Kit,
  seq(100, 4000, by = 100)
) |>
  dplyr::rename('Kit' = Var1, 'depth' = Var2) |>
  dplyr::rowwise() |>
  filter(depth < model_coef_u[model_coef_u$Kit == Kit, 'a']) |>
  mutate(req_reads = readfit(depth, 
                             model_coef_u[model_coef_u$Kit == Kit, 'a'],
                             model_coef_u[model_coef_u$Kit == Kit, 'c'])) |>
  ungroup() |>
  left_join(read_eff, by = 'Kit') |>
  mutate(target_reads = req_reads / eff) |>
  mutate(cost_per_cell = req_reads * seq_cost_coef) |>
  mutate(example_cost = cost_per_cell * 80000)

ggplot(seq_cost_modeling, aes(x=depth, y=cost_per_cell, color = Kit,
                              linetype = case_when(
                                Kit %in% kit_order_3p ~ "WT only",
                                Kit %in% kit_order_5p ~ "TCR"))) + 
  geom_line() +
  scale_color_manual(values = unlist(color_palette$kits), labels = label_function, breaks = kit_order_all) +
  scale_linetype_manual(values = c("TCR" = "dashed", "WT only" = "solid")) +
  scale_y_continuous(labels = scales::label_currency(),
                     sec.axis = sec_axis(~ . * 80000,
                                         labels = scales::label_currency(), 
                                         name = 'Cost for 80k cells')) +
  labs(x='UMI per cell', y='Cost per cell', color = 'Kit', linetype = 'Kit type') +
  coord_cartesian(ylim=c(0,.1)) ->
  figures[['seq_cost_curves_umi']]
my_plot_save(figures[['seq_cost_curves_umi']],
             here('figures/cost/seq_cost_curves_umi.svg'),
             device ='svglite' ,
             width = 7.5, height = 5)


### Gene ----
seq_cost_modeling <- expand.grid(
  model_coef_g$Kit,
  seq(100, 2500, by = 10)
) |>
  dplyr::rename('Kit' = Var1, 'depth' = Var2) |>
  dplyr::rowwise() |>
  filter(depth < model_coef_g[model_coef_g$Kit == Kit, 'a']) |>
  mutate(req_reads = readfit(depth, 
                             model_coef_g[model_coef_g$Kit == Kit, 'a'],
                             model_coef_g[model_coef_g$Kit == Kit, 'c']))|>
  ungroup() |>
  left_join(read_eff, by = 'Kit') |>
  mutate(target_reads = req_reads / eff) |>
  mutate(cost_per_cell = req_reads * seq_cost_coef) |>
  mutate(example_cost = cost_per_cell * 80000)

ggplot(seq_cost_modeling,
       aes(x=depth, y=cost_per_cell, color = Kit,
           linetype = case_when(
             Kit %in% kit_order_3p ~ "WT only",
             Kit %in% kit_order_5p ~ "TCR"))) + 
  geom_line() +
  scale_color_manual(values = unlist(color_palette$kits), labels = label_function,
                     breaks = kit_order_all) +
  scale_linetype_manual(values = c("TCR" = "dashed", "WT only" = "solid")) +
  # scale_linetype_manual(c(setNames(rep('solid', length(kit_order_3p)), kit_order_3p),
  #                         setNames(rep('dashed', length(kit_order_5p)), kit_order_5p))) +
  scale_y_continuous(labels = scales::label_currency(),
                     sec.axis = sec_axis(~ . * 80000,
                                         labels = scales::label_currency(), 
                                         name = 'Cost for 80k cells')) +
  labs(x='Genes per cell', y='Cost per cell', color = 'Kit', linetype = 'Kit type') +
  coord_cartesian(ylim=c(0,.1)) ->
  figures[['seq_cost_curves_gene']]
my_plot_save(figures[['seq_cost_curves_gene']], 
             here('figures/cost/seq_cost_curves_gene.svg'),
             device ='svglite' ,
             width = 7.5, height = 5)

