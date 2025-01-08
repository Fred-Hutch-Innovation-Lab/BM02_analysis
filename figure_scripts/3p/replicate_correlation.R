# Setup ----
source(here('figure_scripts/utils.R'))
library(DT)     ## Datatables
library(DESeq2)  ## Pseudobulk analysis
rbindlist <- data.table::rbindlist

# Load data ----
obj <- readRDS(here('rds/3p/pseudobulk_obj.Rds'))
coarse_celltypes <- data.frame(
  fine = c("CD4+ T", "CD8+ T", 'T',
           "NK",
           "Classical monocyte", "Non-classical monocyte", "Monocyte",
           "Dendritic",
           "B naive", "B memory", 'B',
           "Megakaryocyte",
           "pDC", 
           "Granulocyte", 
           "Erythrocyte",
           "Unknown"),
  coarse =  c("T", "T", 'T',
              "NK",
              "Monocyte", "Monocyte", "Monocyte",
              "Dendritic",
              "B", "B", 'B',
              "Megakaryocyte",
              "pDC", 
              "Granulocyte", 
              "Erythrocyte",
              "Unknown")
)

# Prepare plotdata ----

sampleDists <- dist(t(assays(obj)$vst))

obj_cor <- function(obj, kit, celltype, ngenes = 200, method = 'pearson') {
  plotdata <- subset(obj, select=colData(obj)$Kit==kit)
  plotdata <- subset(plotdata, select=colData(plotdata)$celltype==celltype)
  rv <- rowVars(plotdata@assays@data$vst) 
  var_feat <- names(head(sort(rv, decreasing=TRUE), 200))
  plotdata <- subset(plotdata, rownames(plotdata) %in% var_feat)
  plotdata <- assays(plotdata)$vst |>
    cor(method = 'pearson')
  return(plotdata)
}

# Plot ----
for (kit in unique(metadata_3p$Kit)) {
  plotdata <- tibble()
  # for (celltype in unique(colData(obj)$celltype)) {
    plotdata <- obj_cor(obj, kit=kit, ngenes=200, method='pearson') |>
      melt() |> 
      separate_wider_delim(Var1, '_', too_many = 'merge',
                           names = c('Sample1', 'celltype1'), cols_remove = TRUE) |>
      separate_wider_delim(Var2, '_', too_many = 'merge',
                           names = c('Sample2', 'celltype2'), cols_remove = TRUE) |>
      separate_wider_regex(Sample1, patterns = c(
        Kit1 = '^.+', '-', Individual1 = '[^-]+$'
      ), cols_remove = TRUE) |>
      separate_wider_regex(Sample2, patterns = c(
        Kit2 = '^.+', '-', Individual2 = '[^-]+$'
      ), cols_remove = TRUE) |>
      mutate(Kit1 = gsub('-', '_', Kit1),
             Kit2 = gsub('-', '_', Kit2)) |>
      # filter(celltype1 != 'Unknown') |>
      filter(celltype1 == celltype2) |>
      mutate(celltype1 = factor(celltype1, levels=coarse_celltypes$fine),
             celltype2 = factor(celltype2, levels=coarse_celltypes$fine)) |>
      arrange(celltype1, celltype2) |>
      unique()
    # plotdata <- rbind(plotdata1, plotdata)
  # }
  
  
  ggplot(plotdata,
         aes(x=Individual1,
             y=Individual2,
             label=round(value,2),
             fill=value)) +
    geom_tile() +
    geom_text() +
    scale_fill_gradient2(low = "white",
                         mid="#FFFFCC",
                         high = "#FF0000", 
                         midpoint=0.5) +
    labs(x = 'Sample 1', y = 'Sample 2', fill = 'Corr', title = label_function(kit)) +
    facet_wrap(~ celltype1, drop = FALSE, ncol=4, nrow = 4) ->
    figures[[paste0('corr_mat_pearson', kit)]]
  my_plot_save(image = figures[[paste0('corr_mat_pearson', kit)]],
               path = here('figures/3p/sample_correlation/', paste0(kit, '_sample_correlations.svg')),
               width = 12, height = 12)
}
