library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(Seurat)
library(ggrastr)    ## rasterize layers of ggplot
library(patchwork)  ## arrange plots
library(svglite)
library(data.table)
library(DT)     ## Datatables
library(gt)
# library(extrafont)
library(reshape2)  ## melt
set.seed(33)

systemfonts::register_font('Arial', plain=here('config/Arial.ttf'))
# ttf_import(here('config'), pattern = '.ttf')
# loadfonts(device = 'all', quiet=TRUE)
source(here('config/kit_order.R'))
source(here('config/color_palette.R'))

metadata_3p <- read.csv(here('config/3p/metadata.csv')) |>
  mutate(Kit = factor(Kit, levels = kit_order_3p))

metadata_5p <- read.csv(here('config/5p/metadata.csv')) |>
  mutate(Kit = factor(Kit, levels = kit_order_5p))

if (!exists(quote(figures))) {
  figures <- list()
} else {
  message('Figures object exists in environment, adding figures to existing obj. Be careful of overwrites.')
}

theme_fhil <- function(){ 
  font <- "Arial"   #assign font family up front
  
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      text = element_text(size=16, family=font)
    )
}
theme_set(theme_fhil())

my_plot_save <- function(image, path, width, height, device='svglite'){
  if (!dir.exists(dirname(path))) {
    dir.create(dirname(path), recursive=TRUE)
  }
  svglite::svglite(filename = path,
                   width = unit(width, 'in'), height = unit(height, 'in'))
  plot(image)
  dev.off()
}

write_plot_data <- function(plotdata, file, sep='\t', row.names = FALSE, quote = FALSE, ...){
  if (!dir.exists(dirname(file))) {
    dir.create(dirname(file), recursive=TRUE)
  }
  write.table(plotdata, quote = quote, sep = sep, row.names = row.names, file=file, ...)
}
