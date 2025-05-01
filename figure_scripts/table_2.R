library(tidyverse)
library(gt)
library(here)

source(here('figure_scripts/utils.R'))

plotdata <- read.csv(here('figure_data/table_2.txt'), sep='\t') |> 
  column_to_rownames('Metric') |>
  separate_wider_delim(cols = everything(), delim = "+/-", names_sep=c('_M'), too_few = 'align_start', cols_remove = FALSE) |>
  mutate(across(everything(), str_trim)) |>
  rename_with(~ gsub('_M1', '_mean', .x), dplyr::ends_with('_M1')) |>
  rename_with(~ gsub('_M2', '_se', .x), dplyr::ends_with('_M2')) |>
  mutate(across(dplyr::ends_with('_mean'), parse_number)) |>
  mutate(across(dplyr::ends_with('_se'), parse_number)) 

head(plotdata)

#|>
  mutate(across(c(Cell.Target.Per.Sample,  #Read.1 ,Index.1, Index.2, Read.2, 
           Recommended.Cell.Loading, Instrument.cost, Kit.cost, Kit.Cost.per.sample..n.4., Kit.cost.per.cell),
         ~ as.numeric(gsub("[$,]", "", .))))

figures[['table_1']] <- gt(plotdata) |>
  tab_header(
    title = html('<b>Table 1</b><br>Protocol summaries'),
    subtitle = NA,
  ) |>
  tab_footnote(
    'Does not incuded indices',
    locations= cells_column_labels(c(Kit.Cost.per.sample..n.4., Kit.cost.per.cell))
  ) |> 
  ## Grouping rows
  tab_row_group(
    id = "5'",
    label = "Benchmark 2",
    rows = 8:10
  ) |>
  tab_row_group(
    id = "3'",
    label = 'Benchrmark 1',
    rows = 1:7
  )  |>
  row_group_order(c("3'", "5'")) |>
  ## Column labels
  cols_label(
    Protocol = 'Protocol',
    Kit.Abbreviation = 'Kit ID',
    Assay.Type = 'Assay type',
    Sample.Type = 'Sample type',
    Sample.Process = 'Sample process',
    Detection.Method = 'Detection method',
    Read.1 = 'Read 1',
    Index.1 = 'Index 1',
    Index.2 = 'Index 2',
    Read.2 = 'Read 2',
    Cell.Target.Per.Sample = 'Cell target per sample',
    Recommended.Cell.Loading = 'Recommended cell loading',
    Instrument.cost = 'Instrument cost',
    Kit.cost = 'Kit cost', 
    Kit.Cost.per.sample..n.4. = 'Kit cost per sample (n=4)',
    Kit.cost.per.cell = 'Kit cost per cell'
  ) |>
  ## Formatting columns
  fmt_currency(c(Kit.cost.per.cell), decimals = 2) |>
  fmt_currency(c(Instrument.cost, Kit.cost, Kit.Cost.per.sample..n.4.), decimals = 0) |>
  fmt_integer(c(Cell.Target.Per.Sample, Recommended.Cell.Loading)) |>
  ## Spanners
  tab_spanner(
    "Cost",
    columns = c(
      # Cell.Target.Per.Sample, 
      # Recommended.Cell.Loading,
      Instrument.cost, 
      Kit.cost, 
      Kit.Cost.per.sample..n.4.,
      Kit.cost.per.cell
    )) |>
  tab_spanner(
    'Read length (nt)',
    columns = c(
      Read.1,
      Index.1,
      Index.2,
      Read.2
    )
  )

gtsave(figures[['table_1']], here('figures/tables/table_1.html'))
# my_plot_save(figures[['table_1']], here('figures/tables/table_1.svg'), width = 7, height = 13)
