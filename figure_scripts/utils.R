library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(Seurat)
library(ggrastr)    ## rasterize layers of ggplot
library(patchwork)  ## arrange plots
library(svglite)
library(reshape2)  ## melt
set.seed(33)

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
      text = element_text(size=16)
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
