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

label_function <- function(mode = 'both', wrap = FALSE, width = 20) {
  fun <- function(labels) {
    lapply(labels, function(value){
      x <- case_when(
        grepl('Flex', value) ~ "Flex",
        grepl('NextGEM3P', value) ~ "NextGEM 3'",
        grepl('GEMX3P', value) ~ "GEMX 3'",
        grepl('Fluent_v4', value) ~ "Fluent v4",
        grepl('Fluent_V', value) ~ "Fluent V",
        grepl('Parse_v3', value) ~ "Parse v3",
        grepl('Parse_v2', value) ~ case_when(
          mode == 'TCR' ~ "Parse v2 TCR",
          mode == 'WT' ~ "Parse v2 WT",
          mode == 'both' ~ "Parse v2 WT+TCR",
          mode == 'clean' ~ 'Parse v2'
        ),
        grepl('Scale', value) ~ "Scale",
        grepl('NextGEM5P', value) ~ case_when(
          mode == 'TCR' ~ "NextGEM 5' TCR",
          mode == 'WT' ~ "NextGEM 5' WT",
          mode == 'both' ~ "NextGEM 5' WT+TCR",
          mode == 'clean' ~ "NextGEM 5'"
        ),
        grepl('GEMX5P', value) ~ case_when(
          mode == 'TCR' ~ "GEMX 5' TCR",
          mode == 'WT' ~ "GEMX 5' WT",
          mode == 'both' ~ "GEMX 5' WT+TCR",
          mode == 'clean' ~ "GEMX 5'"
        ),
        .default = value
      )
      if (wrap) {
        x <- str_wrap(x, width = width)
      }
      
      return(x)
    })
  }
  structure(fun, class = 'labeller')
}
