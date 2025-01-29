skip <- c('run_all.R', 'pipeline_metrics_table.R')
files <- dir(here('figure_scripts/3p'), pattern = '\\.R$')
files <- files[!files %in% skip]
for (file in files) {
  tryCatch(source(here('figure_scripts/3p', file)), 
           error = function(e) {print(paste0('Error running file ', file))})
}
