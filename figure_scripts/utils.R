library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(Seurat)
library(ggrastr)    ## rasterize layers of ggplot
library(patchwork)  ## arrange plots
library(svglite)
library(data.table)
library(DT)     ## Datatables
library(gt)
library(janitor)
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
  if (device == 'svglite') {
    svglite::svglite(filename = path,
                     width = unit(width, 'in'), height = unit(height, 'in'))
    plot(image)
    dev.off()
  } else if (device == 'pdf'){
    pdf(file = path, width = unit(width, 'in'), height = unit(height, 'in'))
    plot(image)
    dev.off()
  }
}

prettyprint_vector <- function(console_output, line_width = 80) {
  # Combine input into a single string if it's a vector of lines
  if (length(console_output) > 1) {
    console_output <- paste(console_output, collapse = " ")
  }
  
  # Remove [index] markers
  cleaned <- gsub("\\[\\d+\\]", "", console_output)
  
  # Split by whitespace, preserving quoted strings (e.g., character vectors)
  elements <- strsplit(cleaned, "(?<!\\\\)\\s+(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", perl = TRUE)[[1]]
  elements <- trimws(elements)
  elements <- elements[nzchar(elements)]  # Drop empty strings
  
  # Collapse and wrap
  full <- paste(elements, collapse = ", ")
  formatted <- strwrap(full, width = line_width)
  cat(paste0("c(\n  ", paste(formatted, collapse = "\n  "), "\n)"))
}

write_plot_data <- function(plotdata, file, sep='\t', row.names = FALSE, quote = FALSE, ...){
  if (!dir.exists(dirname(file))) {
    dir.create(dirname(file), recursive=TRUE)
  }
  
  # Convert list columns to character vectors by concatenating elements
  plotdata <- as.data.frame(plotdata) |>
    mutate(across(where(is.list), 
                  ~map_chr(., function(x) {
                    if (is.null(x)) return(NA)
                    paste(as.character(x), collapse = ";")
                  })))
  
  write.table(plotdata, quote = quote, sep = sep, row.names = row.names, file=file, ...)
}
