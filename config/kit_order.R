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
kit_order_all <- c(
  'Flex',
  'NextGEM5P',
  'NextGEM3P',
  'GEMX5P',
  'GEMX3P',
  'Fluent_v4',
  'Fluent_V',
  'Parse_v2',
  'Parse_v3',
  'Scale'
)
label_function <- function(kits) {
  lapply(kits, function(kit) {
    case_when(
      grepl('Flex', kit) ~ "Flex",
      grepl('NextGEM3P', kit) ~ "NextGEM 3'",
      grepl('GEMX3P', kit) ~ "GEMX 3'",
      grepl('Fluent_v4', kit) ~ "Fluent v4",
      grepl('Fluent_V', kit) ~ "Fluent V",
      grepl('Parse_v3', kit) ~ "Parse v3",
      grepl('Parse_v2', kit) ~ "Parse TCR",
      grepl('Scale', kit) ~ "Scale",
      grepl('NextGEM5P', kit) ~ "NextGEM 5'",
      grepl('GEMX5P', kit) ~ "GEMX 5'",
      .default = kit
    )
  }) |>
    unlist() |>
    unname()
}
