## Colors optimized for color-blind visibility
## https://www.nature.com/articles/nmeth.1618
color_palette <- function() {
  palette <- list(
    kits=list(
      NextGEM5P = "#FF7700",
      NextGEM3P = "#f39c00", 
      Fluent_V  = "#56b4e9",
      Fluent_v4 = "#0072b2", 
      GEMX3P =    "#b03800", 
      GEMX5P    = "#E41A1C",
      Parse_v2  = '#3add77',
      Parse_v3 =  "#009e73", 
      Scale =     '#ff8696', 
      Flex =      "#f0c64d" 
  ), cell_colors = list(
    'T'                      = "blue",
    'CD8+ T'                 = "#377EB8",
    'CD4+ T'                 = "#80B1D3",
    'B naive'                = "#FB8072",
    'B memory'               = "#E41A1C",
    'B'                      = 'red',
    'Monocyte'               = "purple",
    'Classical monocyte'     = "#BC80BD",
    'Non-classical monocyte' = "#984EA3",
    'NK'                     = "#4DAF4A",
    'Megakaryocyte'          = "#A65628",
    'Dendritic'              =  "orange",  
    'pDC'                    = "gold", 
    'Granulocyte'            = "#CCEBC5",
    'Erythrocyte'            = "#FCCDE5",
    'Unknown'                = 'grey'
  ),
  samples = list(
    "F1A" = "#117733",
    "F1B" = "#88CCEE",
    "F5A" = "#CC6677",
    "F5B" = "#AA4499"
  ))
  return(palette)
}

 # extras <- c("#CCEBC5", "#4DAF4A")

# "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"
color_palette <- color_palette()
