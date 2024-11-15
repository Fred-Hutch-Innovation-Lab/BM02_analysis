library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)
rbindlist <- data.table::rbindlist

config <- read_yaml(here('config/config.yaml'))
data_dirs <- list(
  "B2-FLEX_10k"="B2-FLEX_10k/full_data",
  "B2-FluentV4_25k"="B2-FluentV4_25k/OutPut_25K",
  "B2-FluentV5_25k"="B2-FluentV5_25k",
  "B2-GEMX3P_25k"="B2-GEMX3P_25k",
  "B2-NextGEM3P_25k"="B2-NextGEM3P_25k",
  "B2-Parse_25k"="B2-Parse_25k/OutPut_25K_V.1.2.1/COMBINED/",
  "B2-Scale_25k"="B2-Scale_25k",
  "B2-BD_WTA_25k"="B2-BD_WTA_25k"
)
samplesheets <- lapply(names(data_dirs), function(x) {
  read.table(here(file.path('config', paste0(x, '_samplesheet.txt'))), header = TRUE)
})
names(samplesheets) <- names(data_dirs)
metadata <- data.table::rbindlist(samplesheets, fill=TRUE)

write.table(metadata, here('config/metadata.csv'), quote=FALSE, row.names = FALSE, sep = ',')
            