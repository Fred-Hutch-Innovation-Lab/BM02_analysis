library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)
rbindlist <- data.table::rbindlist

config <- read_yaml(here('config/3p/config.yaml'))
data_dirs <- list(
  "flex"="flex/full_data_runs",
  "fluent_v4"="fluent_v4/downsampled_runs",
  "fluent_V"="fluent_V/downsampled_runs",
  "gemx_3p"="gemx_3p/downsampled_runs",
  "nextgem_3p"="nextgem_3p/downsampled_runs",
  "parse_v3"="parse_v3/downsampled_runs/COMBINED/",
  "scale"="scale/downsampled_runs/ScaleRNA_out_w_cell_calling"#,
  # "B2-BD_WTA_25k"="B2-BD_WTA_25k"
)
samplesheets <- lapply(names(data_dirs), function(x) {
  read.table(here(file.path('config/3p/', paste0(x, '_samplesheet.txt'))), header = TRUE)
})
names(samplesheets) <- names(data_dirs)
metadata <- data.table::rbindlist(samplesheets, fill=TRUE)

write.table(metadata, here('config/3p/metadata.csv'), quote=FALSE, row.names = FALSE, sep = ',')
            