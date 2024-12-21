library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)
rbindlist <- data.table::rbindlist

config <- read_yaml(here('config/tcr/config.yaml'))
data_dirs <- list(
  "gemx_5p"="gemx_5p/downsampled",
  "nextgem_5p"="nextgem_5p/downsampled",
  "parse_v2"="parse_v2"
)
samplesheets <- lapply(names(data_dirs), function(x) {
  read.table(here(file.path('config/tcr', paste0(x, '_samplesheet.txt'))), header = TRUE)
})
names(samplesheets) <- names(data_dirs)
metadata <- data.table::rbindlist(samplesheets, fill=TRUE)

write.table(metadata, here('config/tcr/metadata.csv'), quote=FALSE, row.names = FALSE, sep = ',')
            