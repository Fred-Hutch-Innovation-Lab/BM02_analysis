## Colors optimized for color-blind visibility
## https://www.nature.com/articles/nmeth.1618

create_color_palette <- function() {
  palette <- c(
    Fluent_v5 =     "#E69F00", #rgb(230, 159,   0, maxColorValue = 255),  # Orange 
    NextGEM3P = "#56B4E9", #rgb(86,  180, 233, maxColorValue = 255),  # Sky blue
    GEMX3P =    "#0072B2", #rgb(0,   114, 178, maxColorValue = 255),  # Blue
    Fluent_v4 =    "#F0E442", #rgb(240, 228,  66, maxColorValue = 255),  # Yellow
    Parse_V3 =  "#009E73", #rgb(0,   158, 115, maxColorValue = 255),  # Bluish green
    Scale =     "#CC79A7", #rgb(204, 121, 167, maxColorValue = 255),  # pink
    Flex =      "#D55E00",  #rgb(213,  94,   0, maxColorValue = 255)   # Vermillion
    BD = 'black',
  )
  return(palette)
}

color_palette <- create_color_palette()
