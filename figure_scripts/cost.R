# Setup ----
source(here('figure_scripts/utils.R')) 
# library(readxl)
library(ggpattern)
library(immunarch)


# Load data ----

## Cell counts ----

cell_loading_data_3p <- read.table(here('figure_data/3p/cell_recovery/cell_recovery_counts.txt'), header=TRUE)
cell_loading_data_5p <- read.table(here('figure_data/5p/wt/cell_recovery/tcell_recovery_counts.txt'), header=TRUE)
cell_loading_data <- merge(cell_loading_data_3p, cell_loading_data_5p, all=TRUE, by=colnames(cell_loading_data_3p)) |> 
  group_by(Kit) |>
  summarize(recovery = sum(after_doublet_filtering) / sum(cells_loaded)) |>
  dplyr::rename(Assay = Kit) |>
  mutate(kit_family = case_when(
    grepl('Fluent', Assay) ~ 'Fluent',
    Assay == 'Flex' ~ 'Flex',
    grepl('NextGEM', Assay) ~ 'NextGEM',
    grepl('GEMX', Assay) ~ 'GEMX',
    Assay == 'Parse_v3' ~ 'Parse_WT',
    Assay == 'Parse_v2' ~ 'Parse_TCR',
    Assay == 'Scale' ~ 'Scale'
  )) 

## Expected recovery ----

expected_recovery_3p <- read.csv(here('data/3p/loaded_cells.csv')) |>
  as.data.table() |>
  melt(id='Sample') |>
  merge(metadata_3p, by='Sample') |>
  group_by(Kit) |>
  filter(variable == 'target_cells_fraction') |>
  summarize(expected_cells_fraction = mean(value)) |>
  mutate(kit_family = case_when(
    grepl('Fluent', Kit) ~ 'Fluent',
    Kit == 'Flex' ~ 'Flex',
    grepl('NextGEM', Kit) ~ 'NextGEM',
    grepl('GEMX', Kit) ~ 'GEMX',
    Kit == 'Parse_v3' ~ 'Parse_WT',
    Kit == 'Parse_v2' ~ 'Parse_TCR',
    Kit == 'Scale' ~ 'Scale'
  ))

expected_recovery_5p <- read.csv(here('data/5p/loaded_cells.csv')) |>
  as.data.table() |>
  melt(id='Sample') |>
  merge(metadata_5p, by='Sample') |>
  group_by(Kit) |>
  filter(variable == 'target_cells_fraction') |>
  summarize(expected_cells_fraction = mean(value)) |>
  mutate(kit_family = case_when(
    grepl('Fluent', Kit) ~ 'Fluent',
    Kit == 'Flex' ~ 'Flex',
    grepl('NextGEM', Kit) ~ 'NextGEM',
    grepl('GEMX', Kit) ~ 'GEMX',
    Kit == 'Parse_v3' ~ 'Parse_WT',
    Kit == 'Parse_v2' ~ 'Parse_TCR',
    Kit == 'Scale' ~ 'Scale'
  ))

expected_recovery <- merge(expected_recovery_3p, expected_recovery_5p, all=TRUE, by=colnames(expected_recovery_3p)) 

cell_loading_data <- merge(cell_loading_data, expected_recovery, by.x=c('Assay', 'kit_family'), by.y=c('Kit', 'kit_family')) |> 
  mutate(cost_ratio = expected_cells_fraction / recovery)


## Kit costs ----
cost_data_ref <- read.table(here('data/cost.txt'), sep='\t', header = TRUE)  |>
  mutate(kit = trimws(kit), 
         cost = readr::parse_number(cost),
         max_cells = readr::parse_number(max_cells)) |>
  
  group_by(kit_family) |>
  arrange(kit_family, max_samples, max_cells) |>
  mutate(kit_throughput_n = as.factor(row_number())) |>
  mutate(baseline_cost_per_cell = cost / max_cells)

cost_data_per_kit <- cost_data_ref |>
  filter(kit %in% c(
    '10X FLEX 4 x 4',
    '10X Next GEM 4rxn',
    '10X GEM-X  4rxn',
    'Fluent (both)',
    'Parse v3- 12 rxns',
    'Parse TCR-WT 4 rxns',
    'Scale ET'
  )) |>
  merge(cell_loading_data, by='kit_family', all=TRUE) |>
  mutate(observed_cost = baseline_cost_per_cell  * cost_ratio) |>
  mutate(kit_family = case_when(grepl('Parse', kit_family) ~ 'Parse', .default = kit_family)) |>
  mutate(kit_family = factor(kit_family, levels = c('Flex', 'NextGEM', 'GEMX', 'Fluent', 'Parse', 'Scale')),
         Assay = factor(Assay, levels = kit_order_all))

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



## Conversion efficiency ----
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

## Clonotypes ----

fig_objs <- readRDS(here('rds/5p/tcr/tcr_paired.Rds'))
clonotype_counts <- data.table::rbindlist(lapply(fig_objs, function(x) {
  x$data |>
    group_by(chain) |>
    summarize(n=n())
}), idcol = 'Sample') |>
  mutate(chain_group = case_when(
    chain == 'TRA' ~ 'partial',
    chain == 'TRB' ~ 'partial',
    chain == 'TRA;TRA' ~ 'partial',
    chain == 'TRB;TRB' ~ 'partial',
    chain == 'TRA;TRB' ~ 'complete',
    .default = 'complex'
  )) |>
  mutate(chain_group = factor(chain_group, levels=c('complete', 'complex', 'partial'))) |>
  group_by(Sample, chain_group) |>
  summarize(n = sum(n)) |>
  mutate(total = cumsum(n))
  
# clonotype_counts$prop <- clonotype_counts$n / unlist(lapply(clonotype_counts$Sample, return_total_cells))


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

## Per clonotype ----
plotdata <- merge(clonotype_counts, metadata_5p, by='Sample') |>
  left_join(cost_data_per_kit, join_by('Kit' =='Assay')) |>
  mutate(cost_per_clonotype_group = cost / total) |>
  select(Sample, chain_group, n, total, Kit, Individual, cost_per_clonotype_group)
ggplot(plotdata, aes(x=chain_group, y=cost_per_clonotype_group, fill=Kit)) +
  geom_boxplot() +
  geom_point(aes(shape = Individual)) +
  facet_wrap(~ Kit, labeller = label_function(mode='TCR')) +
  scale_x_discrete(labels = str_to_title) + 
  scale_fill_manual(values = color_palette$kits, labels = label_function(mode='TCR')) +
  scale_y_continuous(labels = scales::label_currency()) +
  labs(x='Clone type', y='Cost per clone', caption = "Complete = 1 TRA AND 1 TRB\nComplex = 1+ TRA AND 1+ TRB\nPartial = 1+ TRA XOR 1+ TRB") ->
  figures[['price_per_clonotype']]
write_plot_data(plotdata, file = here('figure_data/cost_per_clonotype.txt'))
my_plot_save(figures[['price_per_clonotype']],
             here('figures/cost/cost_per_clonotype.svg'),
             device ='svglite' ,
             width = 12, height = 5)

## Per cell/reagents ----
plotdata <- cost_data_per_kit |>
  select(kit_family, Assay, baseline_cost_per_cell, observed_cost) |>
  as.data.table() 

plotdata |>
data.table::melt(id.vars = c('kit_family', 'Assay')) |>
  mutate(variable = factor(variable, levels = c('baseline_cost_per_cell', 'observed_cost'))) |>
ggplot(aes(x=Assay, y=value, 
             fill=Assay, group=variable )) +
  geom_col_pattern(aes(pattern=variable), position='identity', alpha=1, pattern_key_scale_factor = 0.1) +
  facet_grid(~ kit_family, scales='free_x') +
  scale_x_discrete(labels = label_function(mode='WT')) + 
  scale_fill_manual(values = color_palette$kits, labels = label_function(mode='WT'), breaks = kit_order_all) +
  scale_pattern_manual(values=c('stripe', 'none'), labels = c('Expected', 'Observed'))  +
  scale_y_continuous(labels = scales::label_currency()) +
  guides(fill = guide_legend(override.aes = list(pattern = 'none')),
         pattern = guide_legend(override.aes = list(fill = 'white', color = 'black'))) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x='Kit', y='Cost per cell', fill = 'Kit', pattern = 'Cost group') ->
  figures[['price_per_cell']]

plotdata |>
  # select(Assay, kit, cost, max_cells, baseline_cost_per_cell, expected_cells_fraction, recovery, cost_ratio, observed_cost) |>
  write_plot_data(file = here('figure_data/cost_per_cell_table.txt'))
my_plot_save(figures[['price_per_cell']],
             here('figures/cost/cost_per_cell.svg'),
             device ='svglite' ,
             width = 10, height = 5)

## Sequencing ----

## Table ----

# plotdata <- read_eff |>
#   dplyr::rename('read_conversion_eff' = eff) |>
  # rowwise() |>
  # mutate(max_umi = model_coef_u[[Kit]]['a',1],
  #        rd50 = model_coef_u[[Kit]]['c',1]) |>
  # merge(model_coef_u, by='Kit') |>
plotdata <- cost_data_per_kit |>
  merge(model_coef_u, by.x = 'Assay', by.y='Kit') |>
  dplyr::rename('rd50' = c, 'max_umi' = a) |>
  mutate(req_reads_2000_umi = readfit(2000, max_umi, rd50),
         observed_cells = round(max_cells / cost_ratio)) |>
  ungroup() |>
  # mutate(target_reads_per_cell = req_reads_2000_umi / read_conversion_eff) |>
  mutate(seq_cost_80k_cells = req_reads_2000_umi * ncells * nsamples * seq_cost_coef,
         full_seq_cost = req_reads_2000_umi * observed_cells * seq_cost_coef)  |>
  mutate(Assay = factor(Assay, levels = c(
    'Flex', 'NextGEM5P', 'NextGEM3P', 'GEMX5P', 'GEMX3P', 'Fluent_v4', 'Fluent_V', 'Parse_v2', 'Parse_v3', 'Scale'
  ))) |>
  arrange(Assay) |>
  select(Assay, max_umi, rd50, req_reads_2000_umi, seq_cost_80k_cells, observed_cells, full_seq_cost) 
write_plot_data(plotdata, file = here('figure_data/seq_cost_table.txt'))
gt(plotdata) |>
  cols_label(Assay = 'Kit',
             max_umi = 'Max UMI',
             rd50 = 'rd50',
             req_reads_2000_umi = html('Reads per cell<br>needed for<br>2000 UMI'),
             # read_conversion_eff = html('Read<br>conversion<br>efficiency'),
             # target_reads_per_cell = html('Target<br>seq depth<br>per cell'),
             seq_cost_80k_cells = html('Normalized<br>sequencing<br>cost (80k cells)'),
             observed_cells = html('Expected<br>cell recovery'),
             full_seq_cost = html('Sequencing<br>cost of<br>whole kit')
  ) |>
  fmt_currency(c('seq_cost_80k_cells', 'full_seq_cost'), decimals = 0) |>
  # fmt_number(c('read_conversion_eff'), decimals = 2) |>
  fmt_number(c('max_umi', 'rd50', 'req_reads_2000_umi', 'observed_cells'), decimals = 0)|>
  tab_footnote(
    footnote = paste0('Projected cost based on an estimated sequencing cost of $', seq_cost_coef, ' per read'),
    locations = cells_column_labels(columns = c(full_seq_cost, seq_cost_80k_cells))
  ) ->
  figures[['seq_cost_table']]
figures[['seq_cost_table']]
gt::gtsave(figures[['seq_cost_table']], here('figures/cost/seq_cost.html'))


### UMI ----
seq_cost_modeling <- expand.grid(
  model_coef_u$Kit,
  seq(100, 7500, by = 100)
) |>
  dplyr::rename('Kit' = Var1, 'depth' = Var2) |>
  left_join(model_coef_u, by = 'Kit') |>
  filter(depth < a) |>
  mutate(req_reads = readfit(depth, a, c)) |>
  left_join(read_eff, by = 'Kit') |>
  mutate(cost_per_cell = req_reads * seq_cost_coef) |>
  mutate(example_cost = cost_per_cell * 80000)

ggplot(seq_cost_modeling, aes(x=depth, y=cost_per_cell, color = Kit,
                              linetype = case_when(
                                Kit %in% kit_order_3p ~ "3'",
                                Kit %in% kit_order_5p ~ "5' (WT only)"))) + 
  geom_line() +
  lims(x=c(0,7500)) +
  scale_color_manual(values = unlist(color_palette$kits), labels = label_function(mode='WT'), breaks = kit_order_all) +
  scale_linetype_manual(values = c("3'" = "solid", "5' (WT only)" = "dashed")) +
  scale_y_continuous(labels = scales::label_currency(),
                     sec.axis = sec_axis(~ . * 80000,
                                         labels = scales::label_currency(), 
                                         name = 'Cost for 80k cells'), limits = c(0, .1)) +
  labs(x='UMI per cell', y='Sequencing cost per cell', color = 'Kit', linetype = 'Kit type') ->#+
  # coord_cartesian(ylim=c(0,.1)) ->
  figures[['seq_cost_curves_umi']]
my_plot_save(figures[['seq_cost_curves_umi']],
             here('figures/cost/seq_cost_curves_umi.svg'),
             device ='svglite' ,
             width = 8, height = 5)


### Gene ----
seq_cost_modeling <- expand.grid(
  model_coef_g$Kit,
  seq(100, 3000, by = 10)
) |>
  dplyr::rename('Kit' = Var1, 'depth' = Var2) |>
  left_join(model_coef_g, by = 'Kit') |>
  filter(depth < a) |>
  mutate(req_reads = readfit(depth, a, c)) |>
  left_join(read_eff, by = 'Kit') |>
  mutate(cost_per_cell = req_reads * seq_cost_coef) |>
  mutate(example_cost = cost_per_cell * 80000)

ggplot(seq_cost_modeling,
       aes(x=depth, y=cost_per_cell, color = Kit,
           linetype = case_when(
             Kit %in% kit_order_3p ~ "3'",
             Kit %in% kit_order_5p ~ "5' (WT only)"))) + 
  geom_line() +
  lims(x=c(0,3000)) +
  scale_color_manual(values = unlist(color_palette$kits), labels = label_function(mode='WT'),
                     breaks = kit_order_all) +
  scale_linetype_manual(values = c("3'" = "solid", "5' (WT only)" = "dashed")) +
  # scale_linetype_manual(c(setNames(rep('solid', length(kit_order_3p)), kit_order_3p),
  #                         setNames(rep('dashed', length(kit_order_5p)), kit_order_5p))) +
  scale_y_continuous(labels = scales::label_currency(),
                     sec.axis = sec_axis(~ . * 80000,
                                         labels = scales::label_currency(), 
                                         name = 'Cost for 80k cells'),limits = c(0,.1)) +
  labs(x='Genes per cell', y='Sequencing cost per cell', color = 'Kit', linetype = 'Kit type') ->#+
  # coord_cartesian(ylim=c(0,.1)) ->
  figures[['seq_cost_curves_gene']]
my_plot_save(figures[['seq_cost_curves_gene']], 
             here('figures/cost/seq_cost_curves_gene.svg'),
             device ='svglite' ,
             width = 8, height = 5)

