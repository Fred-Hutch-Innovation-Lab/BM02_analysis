library(tidyverse)
library(gt)
library(here)
library(scales)

source(here('figure_scripts/utils.R'))

plotdata <- read.csv(here('figure_data/table_2.csv')) |> 
  column_to_rownames('Metric')
  # separate_wider_delim(cols = everything(), delim = "±", names_sep=c('_M'), too_few = 'align_start', cols_remove = FALSE) |>

plotdata2 <- plotdata |>
  mutate(across(everything(), str_trim)) |>
  mutate(across(everything(), \(x) str_remove(x, pattern = "±.+"))) |>
  mutate(across(everything(), \(x) str_remove(x, pattern = "\\(.+"))) |>
  mutate(across(everything(), parse_number)) 
  
  # mutate(across(everything(), ~ gsub('(.+)(±.+)*', '\\1', .x)))# |>
  rename_with(~ gsub('_M1', '_mean', .x), dplyr::ends_with('_M1')) |>
  rename_with(~ gsub('_M2', '_se', .x), dplyr::ends_with('_M2')) |>
  mutate(across(dplyr::ends_with('_mean'), parse_number)) |>
  mutate(across(dplyr::ends_with('_se'), parse_number)) 

# #|>
#   mutate(across(c(Cell.Target.Per.Sample,  #Read.1 ,Index.1, Index.2, Read.2, 
#            Recommended.Cell.Loading, Instrument.cost, Kit.cost, Kit.Cost.per.sample..n.4., Kit.cost.per.cell),
#          ~ as.numeric(gsub("[$,]", "", .))))
parse_values <- function(values) {
  values <- str_remove(values, pattern = "±.+")
  values <- str_remove(values, pattern = "\\(.+")
  parse_number(values)
}

range_calc <- function(values, high, low, mid){
  if (!missing(med)) {
    mid = mean(values, na.rm=TRUE)
  } 
  sd <- sqrt(sum((values[!is.na(values)] - mid) ** 2)) / length(values[!is.na(values)])
  if (!missing(low)) {
    low = min(mid - sd*2, min(values, na.rm = TRUE))
  }
}

get_color_palette <- function(invert = FALSE, 
                              midpoint = c("median", "manual"), 
                              manual_mid = NULL, 
                              range_method = c("minmax", "mad"), 
                              sd_range = 4,
                              low_color = "#b03800",
                              mid_color = "lightyellow", 
                              high_color = "#009e73") {
  midpoint <- match.arg(midpoint)
  range_method <- match.arg(range_method)
  
  # Return a function to be called on the vector (e.g., by gt::data_color)
  function(x) {
    # x <- suppressWarnings(as.numeric(x))  # Ensure numeric and remove non-numeric silently
    x <- parse_values(x)
    x_na_omit <- na.omit(x)
    if (!is.numeric(x_na_omit)) stop("get_color_palette can only be applied to numeric data")
    
    # Compute midpoint
    mid_val <- if (midpoint == "median") {
      median(x_na_omit)
    } else {
      if (is.null(manual_mid)) stop("Manual midpoint must be provided when midpoint = 'manual'")
      manual_mid
    }
    
    # Compute range
    if (range_method == "minmax") {
      low_val <- min(x_na_omit)
      high_val <- max(x_na_omit)
    } else {
      mad_val <- mad(x_na_omit, constant = 1)
      median_val <- median(x_na_omit)
      low_val <- median_val - sd_range * mad_val
      high_val <- median_val + sd_range * mad_val
    }
    
    # Invert colors if requested
    if (invert) {
      tmp <- low_color
      low_color <- high_color
      high_color <- tmp
    }
    
    # Build the color palette
    color_scale <- col_numeric(
      palette = c(low_color, mid_color, high_color),
      domain = c(low_val, high_val),
      na.color = 'lightgrey'
    )
    
    x_squished <- pmax(pmin(x, high_val), low_val)
    # Return color-mapped vector
    color_scale(x_squished)
  }
}

  
  # Build gt table
plotdata[c('Flex', 'NextGEM', 'GEM-X', 'Fluent v4', 'Fluent V', 'Parse v3', 'Scale'),] |>
  select(-c(cost_per_recovered_paired_clone, proportion_of_cells_with_paired_clones)) |>
  gt(rownames_to_stub = TRUE) |>
  # tab_header(
  #   title = html('<b>Table 2</b><br>Result summaries'),
  #   subtitle = NA,
  # ) |>
  ## labeling
  cols_label(
    high_quality_cell_recovery = "High quality cell recovery",
    cost_per_cell_recovered = "Cost per high quality cell",
    reads_mapped_to_transcriptome = "Reads mapped to transcriptome", 
    reads_mapped_and_in_called_cells = "Reads contributing to counts matrix",
    maximum_gene_recovery_estimate = "Estimatated maximum gene recovery",
    rd50_gene = "rd50 gene",
    maximum_umi_recovery_estimate = "Estimated maximum transcript recovery",
    rd50_umi = "rd50 transcripts",
    median_genes_per_cell = "Median genes per cell",
    median_umi_per_cell = "Median transcripts per cell",
    median_ribosomal_umi_proportion = "Median ribosomal proportion",
    median_mitochondrial_umi_proportion = "Median mitochondrial proportion",
    informative_gene_proportion = "Reads mapping to informative genes",
    earth_mover_distance = "Earth mover distance",
    percent_residual_variance_in_gene_expression = "Residual variance in gene expression",
    successfully_detected_enriched_pathways = "Successfully detected enriched pathways"
    # cost_per_recovered_paired_clone = "Cost per recovered clone pair",
    # proportion_of_cells_with_paired_clones = "Cells with paired clones"
  ) |>
  ## Format cols
  fmt_currency(c(cost_per_cell_recovered), decimals = 2) |>
  fmt_integer(c(maximum_gene_recovery_estimate, rd50_gene, maximum_umi_recovery_estimate, rd50_umi)) |>
  ## coloring
  data_color(
    columns = c(
      rd50_gene,
      rd50_umi,
      median_ribosomal_umi_proportion,
      median_mitochondrial_umi_proportion,
      percent_residual_variance_in_gene_expression,
      earth_mover_distance
      # cost_per_recovered_paired_clone
    ),
    fn = get_color_palette(midpoint='median', range_method='mad', invert=TRUE)
  ) |>
  data_color(
    columns = c(
      high_quality_cell_recovery,
      reads_mapped_to_transcriptome,
      reads_mapped_and_in_called_cells,
      maximum_gene_recovery_estimate,
      maximum_umi_recovery_estimate,
      median_genes_per_cell,
      median_umi_per_cell,
      informative_gene_proportion,
      successfully_detected_enriched_pathways
      # proportion_of_cells_with_paired_clones
    ),
    fn = get_color_palette(midpoint='median', range_method='mad', invert=FALSE)
  ) |> data_color(
    columns = c(
      cost_per_cell_recovered
    ),
    fn = get_color_palette(midpoint='median', range_method='mad', sd_range = 5, invert=TRUE)
  ) |>
  ## Spanners
  tab_spanner("Cell recovery",columns = c(high_quality_cell_recovery, cost_per_cell_recovered)) |>
  tab_spanner("Sequencing efficiency",columns = c(
    reads_mapped_to_transcriptome, reads_mapped_and_in_called_cells, maximum_gene_recovery_estimate,
    rd50_gene, maximum_umi_recovery_estimate, rd50_umi
  )) |>
  tab_spanner("Feature recovery",columns = c(
    median_genes_per_cell, median_umi_per_cell,
    median_ribosomal_umi_proportion, median_mitochondrial_umi_proportion,
    informative_gene_proportion
  )) |>
  tab_spanner("Replicability",columns = c(
    earth_mover_distance, percent_residual_variance_in_gene_expression,
  )) |> 
  tab_style(
    style = cell_borders(sides = "left", color = "black", weight = px(6), style = "solid"),
    locations = cells_body(
      columns = c(reads_mapped_to_transcriptome, median_genes_per_cell, earth_mover_distance))
  ) -> 
  figures[['table_2']]
figures[['table_2']]

gtsave(figures[['table_2']], here('figures/tables/table_2.html'))
# my_plot_save(figures[['table_1']], here('figures/tables/table_1.svg'), width = 7, height = 13)
