library(readxl)

convert_xlsx_to_csv <- function(input_file, sheet = 1, output_file = NULL) {
  data <- read_excel(input_file, sheet = sheet)
  
  if (is.null(output_file)) {
    output_file <- sub("\\.xlsx$", ".csv", input_file)
  }
  
  write.csv(data, file = output_file, row.names = FALSE)
  
  message("Conversion complete: ", output_file)
}

for (file in list.files(here('data/saturation_curves'),pattern = '.xlsx$', full.names = TRUE)) {
  convert_xlsx_to_csv(file) 
}
