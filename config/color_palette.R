## Colors optimized for color-blind visibility
## https://www.nature.com/articles/nmeth.1618

create_color_palette <- function() {
  palette <- c(
    other = "#E69F00",      #rgb(230, 159,   0, maxColorValue = 255),  # Orange 
    Fluent = "#56B4E9",    #rgb(86,  180, 233, maxColorValue = 255),  # Sky blue
    GEMX3P = "#009E73",    #rgb(0,   158, 115, maxColorValue = 255),  # Bluish green
    NextGEM3P = 'black', #"#F0E442", #rgb(240, 228,  66, maxColorValue = 255),  # Yellow
    Parse_V3 = "#0072B2",  #rgb(0,   114, 178, maxColorValue = 255),  # Blue
    Scale = "#CC79A7",     # rgb(204, 121, 167, maxColorValue = 255)
    Flex = "#D55E00"      #rgb(213,  94,   0, maxColorValue = 255)   # Vermillion
  )
  return(palette)
}

# Call the function and store the color palette in a variable
color_palette <- create_color_palette()
