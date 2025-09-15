library(Seurat)
library(SeuratDisk)
library(dplyr)

h5_convert <- function(obj, name, overwrite = TRUE, ...) {
  ## Convert to v3 assay
  obj[["RNA"]] <- GetAssayData(obj[["RNA"]]) |>
    CreateAssayObject()
  
  SaveH5Seurat(obj, filename = paste0(name, '.h5seurat'), overwrite = overwrite, ...) 
  Convert(paste0(name, '.h5seurat'), dest = 'h5ad', overwrite = overwrite)
}

objs_5p <- readRDS(here('rds/5p/wt/05_merged_objs_post_clustering.rds'))
lapply(names(objs_5p), function(obj) {
  saveRDS(objs_5p[[obj]], file = here('rds/publish/', paste0('bm02-', obj, '-merged.rds')), compress = TRUE)
})
rm(objs_5p)
gc()

objs_3p <- readRDS(here('rds/3p/07_post_module_scoring.rds'))
lapply(names(objs_3p), function(obj) {
  saveRDS(objs_3p[[obj]], file = here('rds/publish/', paste0('bm01-', obj, '-merged.rds')), compress = TRUE)
})
lapply(names(objs_3p), function(obj) {
  x <- DietSeurat(objs_3p[[obj]],
             assays = 'RNA',
             layers = c('data', 'counts'), 
             dimreducs = c('PCA_allgenes', 'UMAP_allgenes'))
  x@meta.data <- x@meta.data |>
    mutate(across(where(is.factor), as.character))
  h5_convert(x, here('rds/publish/', paste0('bm01-', obj)))
  # saveRDS(objs_5p[[obj]], file = here('rds/publish/', paste0('bm02-', obj, '-merged.rds')), compress = TRUE)
})
