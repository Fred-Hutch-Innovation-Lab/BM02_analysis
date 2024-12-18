## Colors optimized for color-blind visibility
## https://www.nature.com/articles/nmeth.1618

create_color_palette <- function() {
  palette <- list(
    kits=list(
      NextGEM3P = "#E69F00", #rgb(230, 159,   0, maxColorValue = 255),  # Orange 
      Fluent_V = "#56B4E9", #rgb(86,  180, 233, maxColorValue = 255),  # Sky blue
      Fluent_v4 = "#0072B2", #rgb(0,   114, 178, maxColorValue = 255),  # Blue
      GEMX3P =    "#F0E442", #rgb(240, 228,  66, maxColorValue = 255),  # Yellow
      Parse_v2 = 'darkgreen',
      Parse_v3 =  "#009E73", #rgb(0,   158, 115, maxColorValue = 255),  # Bluish green
      Scale =     "#CC79A7", #rgb(204, 121, 167, maxColorValue = 255),  # pink
      Flex =      "#D55E00",  #rgb(213,  94,   0, maxColorValue = 255)   # Vermillion
      BD = 'black'
  ), cell_colors = list(
  #   coarse=c(     
  #   "T" = "blue",
  #   "B" = "red",
  #   "Monocyte" = "purple",
  #   "NK" = "#4DAF4A",
  #   "Megakaryocyte" = "#A65628",
  #   "Dendritic" = "orange",
  #   "pDC" = "gold",
  #   'Unknown' = 'grey'
  # ), fine = c(
    'T' = "blue",
    'CD8+ T' = "#377EB8",
    'CD4+ T' = "#80B1D3",
    'B naive' = "#FB8072",
    'B memory' = "#E41A1C",
    'Monocyte' = "purple",
    'Classical monocyte' = "#BC80BD",
    'Non-classical monocyte' = "#984EA3",
    'NK' = "#4DAF4A",
    'Megakaryocyte' = "#A65628",
    'Dendritic' =  "orange", 
    'pDC' = "gold", 
    'Granulocyte' = "#CCEBC5",
    'Erythrocyte' = "#FCCDE5",
    'Unknown' = 'grey'
  ),
  samples = list(
    "F1A" = "#117733",
    "F1B" = "#88CCEE",
    "F5A" = "#CC6677",
    "F5B" = "#AA4499"
  ))
  return(palette)
}

 extras <- c("#CCEBC5",  "#FB8072",  "#FCCDE5", "#4DAF4A")

color_palette <- create_color_palette()
