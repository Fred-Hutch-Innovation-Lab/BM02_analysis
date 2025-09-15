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
  'NextGEM3P',
  'NextGEM5P',
  'GEMX3P',
  'GEMX5P',
  'Fluent_v4',
  'Fluent_V',
  'Parse_v3',
  'Parse_v2',
  'Scale'
)

label_function <- function(mode = 'both', wrap = FALSE, width = 20) {
  fun <- function(labels) {
    lapply(labels, function(value){
      x <- case_when(
        grepl('Flex', value) ~ "NextGEM Flex",
        grepl('NextGEM3P', value) ~ "NextGEM 3'",
        grepl('GEMX3P', value) ~ "GEM-X 3'",
        grepl('Fluent_v4', value) ~ "PIPseq v4",
        grepl('Fluent_V', value) ~ "PIPseq V",
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
          mode == 'TCR' ~ "GEM-X 5' TCR",
          mode == 'WT' ~ "GEM-X 5' WT",
          mode == 'both' ~ "GEM-X 5' WT+TCR",
          mode == 'clean' ~ "GEM-X 5'"
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
