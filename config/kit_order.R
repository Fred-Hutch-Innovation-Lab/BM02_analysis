kit_order_3p <- c(
  'Flex',
  'NextGEM3P',
  'GEMX3P',
  'Fluent_v4',
  'Fluent_V',
  'Parse_v3',
  'Scale'
)
kit_order_5p <- c(
  'NextGEM5P',
  'GEMX5P',
  'Parse_v2'
)
label_function <- function(kits) {
  lapply(kits, function(kit) {
    case_when(
      kit == 'NextGEM3P' ~ "NextGEM 3'",
      kit == 'GEMX3P' ~ "GEMX 3'",
      kit == 'Fluent_v4' ~ "Fluent v4",
      kit == 'Fluent_V' ~ "Fluent V",
      kit == 'Parse_v3' ~ "Parse v3",
      kit == 'Parse_v2' ~ "Parse v2",
      kit == 'NextGEM5P' ~ "NextGEM 5'",
      kit == 'GEMX5P' ~ "GEMX 5'",
      .default = kit
    )
  }) |>
    unlist() |>
    unname()
}