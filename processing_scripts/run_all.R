library(here)

# for (file in c(
#   # '01-soup_channels.rds',
#   # # '01_raw_objs.rds',
#   # '02-objs_post_cell_filtering.rds',
#   # '03-objs_post_gene_filtering.rds',
#   # "04_objs_post_annotation.rds"
# )) {
#   objs <- readRDS(here(file.path('rds', file)))
#   F1A <- objs$Scale_F5A
#   F1B <- objs$Scale_F5B
#   F5A <- objs$Scale_F1A
#   F5B <- objs$Scale_F1B
#   
#   objs$Scale_F5A <- F5A
#   objs$Scale_F5B <- F5B
#   objs$Scale_F1A <- F1A
#   objs$Scale_F1B <- F1B
#   saveRDS(objs, here(file.path('rds', file)), compress = FALSE)
# }
# 
# for (file in c(
#   '05_merged_objs_post_clustering.rds',
#   '07_post_module_scoring.rds'
# )) {
#   objs <- readRDS(here(file.path('rds', file)))
#   
#   objs$Scale$orig.ident <- plyr::mapvalues(
#     objs$Scale$orig.ident,
#     from = c('Scale_F1A', 'Scale_F1B', 'Scale_F5A', 'Scale_F5B'),
#     to = c('Scale_F5A', 'Scale_F5B', 'Scale_F1A', 'Scale_F1B')
#   )
#   saveRDS(objs, here(file.path('rds', file)), compress = FALSE)
# }




rmarkdown::render(here('processing_scripts/01-ambient_RNA_check.runfile.Rmd'))
rmarkdown::render(here('processing_scripts/02-cell_filtering.runfile.Rmd'))
rmarkdown::render(here('processing_scripts/03-gene_analysis.runfile.Rmd'))
rmarkdown::render(here('processing_scripts/04-reference_annotation.runfile.Rmd'))
rmarkdown::render(here('processing_scripts/05-clustering.runfile.Rmd'))

