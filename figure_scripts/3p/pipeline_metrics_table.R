library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(readxl)
library(kableExtra)
library(formattable) ## formatting percentages


# main <- function() {
#   metadata <- read.csv(here('config/3p/metadata.csv'), stringsAsFactors = TRUE)
#   metrics <- read_xlsx(here('data/3p/pipeline_summary_statistics/downsampled_data.xlsx'), skip = 1) %>%
#     filter(!METRICS %in% c('Last updated')) %>%
#     # mutate(across(-c('METRICS'), as.numeric)) %>%
#     column_to_rownames('METRICS') %>%
#     t() %>% 
#     as.data.frame() %>%
#     rownames_to_column('Sample') %>%
#     as_tibble() %>%
#     mutate(across(-c('Sample'), ~ accounting(.x, format='d'))) %>%
#     mutate(across(c('Sequencing saturation', matches('%')), percent))
#   metrics <- merge(metrics, metadata, by='Sample') %>%
#     mutate(Individual.Replicate = paste0(Individual, Replicate))
  
  
#   cols <- c('Kit', 'Individual.Replicate', #'Individual.Replicate',
#             'Sequencing saturation', '% reads w/ valid cell barcodes',
#             '% Reads mapping to transcriptome', 'Median Transcripts/cell',
#             'Median genes/cell') #, 'Total genes detected')
  
#   plotdata <- metrics %>%
#     select(all_of(cols)) %>% 
#     arrange(Kit, Individual.Replicate)
  
#   qc_metrics_table <- plotdata %>%
#     select(-Kit) %>%
#     knitr::kable(row.names=FALSE, col.names = c('', cols[3:length(cols)])) %>%
#     group_rows('Flex', 1,4) %>%
#     group_rows('Fluent', 5,8) %>%
#     group_rows('GEMX3P', 9,12) %>%
#     group_rows('NextGEM3P', 13,16) %>%
#     group_rows('Parse_V3', 17,20) %>%
#     group_rows('Scale', 21,24) 
#   return(qc_metrics_table)
# }

# qc_metrics_table <- main()
# # library(webshot2)   ## Save kable to png
# # as_image(qc_metrics_table, width=7, height=5, file=here('figures/metrics_table.png'))
