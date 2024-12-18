
skip <- c('run_all.R')
files <- dir(here('figure_scripts'), pattern = '\\.R$')
files <- files[!files %in% skip]
for (file in files) {
  tryCatch(source(here('figure_scripts', file)), 
           error = function(e) {print(paste0('Error running file ', file))})
}
