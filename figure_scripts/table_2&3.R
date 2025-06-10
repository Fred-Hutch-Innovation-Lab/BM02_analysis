library(tidyverse)
library(gt)
library(here)
library(scales)

source(here('figure_scripts/utils.R'))

plotdata <- read.csv(here('figure_data/table_2.csv')) |>
  # column_to_rownames('Metric') |>
  # tibble() #|>
  separate(earth_mover_distance, sep = '-', into = c('emd_low', 'emd_high')) |>
  separate_wider_regex(percent_residual_variance_in_gene_expression, c(
    percent_residual_variance_in_gene_expression = '^\\d+',
    '\\s*\\(',
    percent_residual_variance_in_gene_expression_iqr1 = '\\d+',
    '-',
    percent_residual_variance_in_gene_expression_iqr2 = '\\d+',
    '.+'
  ), too_few = "align_start") |>
  separate_wider_delim(cols = everything(), delim = '±', names_sep='_M', too_few = 'align_start', cols_remove=TRUE) |>
  dplyr::rename_with(~ gsub('_M1', '', .x), dplyr::ends_with('_M1')) |>
  dplyr::rename_with(~ gsub('_M2', '_se', .x), dplyr::ends_with('_M2')) |>
  column_to_rownames('Metric') |>
  mutate(across(everything(), parse_number)) |>
  mutate(projected_umi_rd50 = round(maximum_umi_recovery_estimate / 2),
         projected_gene_rd50 = round(maximum_gene_recovery_estimate / 2)) 
  # separate_wider_delim(cols = everything(), delim = "±", names_sep=c('_M'), too_few = 'align_start', cols_remove = FALSE) |>
head(plotdata)
  
parse_values <- function(values) {
  # Handle vectors by applying the function to each element
  if (length(values) > 1) {
    return(sapply(values, parse_values))
  }
  
  # Handle single value
  if (grepl('-', values)) {
    # For range values, return the midpoint
    range_vals <- str_split(values, '-')[[1]]
    return(mean(sapply(range_vals, parse_number)))
  }
  
  # Remove ± and () parts and parse
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

# BM1 ----
plotdata[c('Flex', "NextGEM 3'", "GEM-X 3'", 'Fluent v4', 'Fluent V', 'Parse v3', 'Scale'),] |>
  gt(rownames_to_stub = TRUE) |>
  ## Hide columns ----
  cols_hide(columns = c(
    rd50_gene,rd50_umi,
    cost_per_recovered_paired_clone, proportion_of_cells_with_paired_clones, 
    median_mitochondrial_umi_proportion, median_ribosomal_umi_proportion,
    # emd_high, emd_low,
    ends_with('_se'), ends_with('_iqr1'), ends_with('_iqr2')
  )) |>
  ## Merge cols ----
  cols_merge_range(
    col_begin = emd_low,
    col_end = emd_high
  ) |>
  # tab_header(
  #   title = html('<b>Table 2</b><br>Result summaries'),
  #   subtitle = NA,
  # ) |>
  ## labeling ----
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
    informative_gene_proportion = "Transcripts mapping to informative genes",
    emd_low = "Earth mover distance",
    percent_residual_variance_in_gene_expression = "Residual variance in gene expression",
    successfully_detected_enriched_pathways = "Successfully detected enriched pathways"
    # cost_per_recovered_paired_clone = "Cost per recovered clone pair",
    # proportion_of_cells_with_paired_clones = "Cells with paired clones"
  ) |>
  ## Format cols ----
  fmt_currency(c(cost_per_cell_recovered), decimals = 2) |>
  fmt_percent(scale_values=FALSE, decimals = 0, c(
    high_quality_cell_recovery, 
    reads_mapped_to_transcriptome, reads_mapped_and_in_called_cells,
    informative_gene_proportion, percent_residual_variance_in_gene_expression
    )) |>
  fmt_integer(c(
    maximum_gene_recovery_estimate, rd50_gene, projected_gene_rd50,
    maximum_umi_recovery_estimate, rd50_umi, projected_umi_rd50,
    median_genes_per_cell, median_umi_per_cell, 
    )) |>
  fmt_engineering(drop_trailing_zeros = FALSE, drop_trailing_dec_mark=TRUE, decimals = 1, c( #, exp_style = 'e1'
    maximum_gene_recovery_estimate, rd50_gene, projected_gene_rd50,
    maximum_umi_recovery_estimate, rd50_umi, projected_umi_rd50,
    median_genes_per_cell, median_umi_per_cell, 
  )) |>
  cols_width(
    everything() ~ px(120)
  ) |>
  ## Merge cols ----
  cols_merge(
    columns = c(rd50_gene, projected_gene_rd50),
    pattern = "{1} reads ({2} genes)"
  ) |>
  cols_merge(
    columns = c(rd50_umi, projected_umi_rd50),
    pattern = "{1} reads ({2} UMIs)"
  ) |>
  ## Add text to cols ----
  text_transform(
    fn = \(x) {paste0(x, "/6")},
    locations = cells_body(columns =successfully_detected_enriched_pathways)
  ) |>
  ## coloring ----
  data_color(
    columns = c(
      cost_per_cell_recovered,
      rd50_gene,
      rd50_umi,
      median_ribosomal_umi_proportion,
      median_mitochondrial_umi_proportion,
      percent_residual_variance_in_gene_expression,
      emd_low
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
  # ) |> data_color(
  #   columns = c(
  #     cost_per_cell_recovered
  #   ),
  #   fn = get_color_palette(midpoint='median', range_method='mad', sd_range = 5, invert=TRUE)
  ) |>
  ## Spanners ----
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
    emd_low, percent_residual_variance_in_gene_expression,
  )) |> 
  tab_style(
    style = cell_borders(sides = "left", color = "black", weight = px(8), style = "solid"),
    locations = cells_body(
      columns = c(reads_mapped_to_transcriptome, median_genes_per_cell, emd_low))
  ) |>
  ## Title
  tab_header(
    title = "Table 2: Summary metrics for benchmark 1"
  ) |> 
  opt_align_table_header(align = c("left")) -> 
  figures[['table_2']]
figures[['table_2']]

gtsave(figures[['table_2']], here('figures/tables/table_2.html'))
# my_plot_save(figures[['table_1']], here('figures/tables/table_1.svg'), width = 7, height = 13)

# BM2 ----
plotdata[c("NextGEM 5'", "GEM-X 5'", 'Parse v2'),] |>
  gt(rownames_to_stub = TRUE) |>
  ## Hide columns ----
cols_hide(columns = c(
  rd50_gene,rd50_umi,
  # cost_per_recovered_paired_clone, proportion_of_cells_with_paired_clones, 
  median_mitochondrial_umi_proportion, median_ribosomal_umi_proportion,
  emd_high, emd_low, percent_residual_variance_in_gene_expression, successfully_detected_enriched_pathways,
  ends_with('_se'), ends_with('_iqr1'), ends_with('_iqr2')
)) |>
  ## Merge cols ----
  # tab_header(
  #   title = html('<b>Table 2</b><br>Result summaries'),
  #   subtitle = NA,
  # ) |>
  ## labeling ----
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
  informative_gene_proportion = "Transcripts mapping to informative genes",
  cost_per_recovered_paired_clone = "Cost per recovered clone pair",
  proportion_of_cells_with_paired_clones = "Cells with paired clones"
) |>
  ## Format cols ----
fmt_currency(c(cost_per_cell_recovered, cost_per_recovered_paired_clone), decimals = 2) |>
  fmt_percent(scale_values=FALSE, decimals = 0, c(
    high_quality_cell_recovery, 
    reads_mapped_to_transcriptome, reads_mapped_and_in_called_cells,
    informative_gene_proportion, proportion_of_cells_with_paired_clones
  )) |>
  fmt_integer(c(
    maximum_gene_recovery_estimate, rd50_gene, projected_gene_rd50,
    maximum_umi_recovery_estimate, rd50_umi, projected_umi_rd50,
    median_genes_per_cell, median_umi_per_cell, 
  )) |>
  fmt_engineering(drop_trailing_zeros = FALSE, drop_trailing_dec_mark=TRUE, decimals = 1, c( #, exp_style = 'e1'
    maximum_gene_recovery_estimate, rd50_gene, projected_gene_rd50,
    maximum_umi_recovery_estimate, rd50_umi, projected_umi_rd50,
    median_genes_per_cell, median_umi_per_cell, 
  )) |>
  cols_width(
    everything() ~ px(120)
  ) |>
  ## Merge cols ----
cols_merge(
  columns = c(rd50_gene, projected_gene_rd50),
  pattern = "{1} reads ({2} genes)"
) |>
  cols_merge(
    columns = c(rd50_umi, projected_umi_rd50),
    pattern = "{1} reads ({2} UMIs)"
  ) |>
  ## coloring ----
data_color(
  columns = c(
    cost_per_cell_recovered,
    cost_per_recovered_paired_clone,
    rd50_gene,
    rd50_umi,
    median_ribosomal_umi_proportion,
    median_mitochondrial_umi_proportion
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
      proportion_of_cells_with_paired_clones
    ),
    fn = get_color_palette(midpoint='median', range_method='mad', invert=FALSE)
    # ) |> data_color(
    #   columns = c(
    #     cost_per_cell_recovered
    #   ),
    #   fn = get_color_palette(midpoint='median', range_method='mad', sd_range = 5, invert=TRUE)
  ) |>
  ## Spanners ----
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
  tab_spanner(
    'TCR recovery', columns = c(cost_per_recovered_paired_clone, proportion_of_cells_with_paired_clones)
  ) |>
  tab_style(
    style = cell_borders(sides = "left", color = "black", weight = px(8), style = "solid"),
    locations = cells_body(
      columns = c(reads_mapped_to_transcriptome, median_genes_per_cell, cost_per_recovered_paired_clone))
  ) |>
  ## Title
  tab_header(
    title = "Table 3: Summary metrics for benchmark 2"
  ) |> 
  opt_align_table_header(align = c("left")) -> 
  figures[['table_3']]
figures[['table_3']]

gtsave(figures[['table_3']], here('figures/tables/table_3.html'))
# my_plot_save(figures[['table_1']], here('figures/tables/table_1.svg'), width = 7, height = 13)
