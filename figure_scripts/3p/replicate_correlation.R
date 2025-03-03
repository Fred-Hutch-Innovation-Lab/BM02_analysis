# Setup ----
source(here('figure_scripts/utils.R'))
library(DESeq2)  ## Pseudobulk analysis

# Load data ----
obj <- readRDS(here('rds/3p/pseudobulk_obj.Rds'))
cellcounts <- read.csv(here('rds/3p/pb_cell_counts.txt'), sep='\t') |> 
  mutate(ind = case_when(
    grepl('F5', Sample) ~ 'F5',
    grepl('F1', Sample) ~ 'F1'
  )) |>
  group_by(Kit, ind, Celltype) |>
  summarize(Count = mean(Count)) |>
  mutate(prop = Count / sum(Count))
label_order <- c(
  "B naive", "B memory",
  "CD4+ T", "CD8+ T", 
  "NK",
  "Classical monocyte", "Non-classical monocyte",
  "Dendritic", "pDC",
  "Megakaryocyte", 'Erythrocyte','Granulocyte'
)

# Functions ----

corr_dotplot <- function(corrdata) {
  corrdata |>
    ggplot(aes(x=Kit, y=Celltype)) +
    geom_point(aes(size=prop, color=y)) +
    # geom_text(aes(label=round(y, 2)), size=2) +
    scale_size_continuous(range=c(1,12), breaks = c(0.05, 0.2, 0.35)) +
    scale_color_gradient2(high = "#4DAF4A",
                          mid="grey",
                          low = "#FF0000", 
                          midpoint=0.5, oob=scales::squish) +
    scale_x_discrete(labels = label_function) +
    # theme_bw() +
    labs(x = 'Kit', y = 'Celltype',
         color = 'Average\nreplicate\ndistance', 
         size = 'Proportion\nof cells\nfrom kit')
}

# Prepare plotdata ----
sampleDists <- dist(t(assay(assays(obj)$vsd))) |>
  as.matrix() |>
  as.data.table(keep.rownames = TRUE) |>
  # rownames_to_column('Sample1') |>
  melt() |>
  separate_wider_delim(rn, '_', too_many = 'merge',
                       names = c('Sample1', 'celltype1'), cols_remove = TRUE) |>
  separate_wider_delim(variable, '_', too_many = 'merge',
                       names = c('Sample2', 'celltype2'), cols_remove = TRUE) |>
  separate_wider_regex(Sample1, patterns = c(
    Kit1 = '^.+', '-', Individual1 = '[^-]+$'
  ), cols_remove = TRUE) |>
  separate_wider_regex(Sample2, patterns = c(
    Kit2 = '^.+', '-', Individual2 = '[^-]+$'
  ), cols_remove = TRUE) |>
  mutate(Kit1 = gsub('-', '_', Kit1),
         Kit2 = gsub('-', '_', Kit2)) |>
  filter(celltype1 == celltype2) |>
  filter(Kit1 == Kit2)  |>
  filter((Individual1 == 'F1A' & Individual2 == 'F1B') | 
           (Individual1 == 'F5A' & Individual2 == 'F5B')) |>
  dplyr::rename(Kit = 'Kit1', 'Celltype' = celltype1) |>
  mutate(Kit = factor(Kit, levels = kit_order_3p),
         Celltype = factor(Celltype, levels = label_order)) |>
  mutate(ind = case_when(
    grepl('F5', Individual1) ~ 'F5',
    grepl('F1', Individual1) ~ 'F1'
  )) |>
  merge(cellcounts, by=c('Kit', 'Celltype', 'ind')) 

# Plot ----
plotdata <- sampleDists  |>
  # mutate(value = value / max(value)) |> 
  group_by(Kit, Celltype) |>
  summarize(y=round(mean(value)), prop = mean(prop)) #, prop = Count / sum(Count)


corr_dotplot(plotdata) +
  geom_hline(yintercept = 9.5, lty='dashed') +
  scale_color_gradient2(low = "#4DAF4A",
                        mid="grey",
                        high = "#FF0000", midpoint = 75, limits=c(0,150), oob=scales::squish,
                        breaks = c(0,75,150), labels = c('0', '75', '150+')) ->
  figures[['sample_distance_dotplot']]
write_plot_data(plotdata, here('figure_data/3p/replicate_correlation/rep_dist_dotplot.txt'))
my_plot_save(image = figures[['sample_distance_dotplot']], 
             path = here('figures/3p/replicate_correlation/rep_dist_dotplot.svg'), 
             width = 10, height = 5.5)

plotdata <- sampleDists |>
  filter(!Celltype %in% c('Erythrocyte', 'Granulocyte', 'Megakaryocyte')) |>
  merge(
    (cellcounts |> 
       group_by(Celltype) |> 
       summarize(Count2 = sum(Count)) |> 
       filter(!Celltype %in% c('Erythrocyte', 'Granulocyte', 'Megakaryocyte')) |>
       mutate(prop2 = Count2 / sum(Count2))
     ), .by = 'Celltype', all.x=TRUE
    ) |>
  # group_by(Celltype) |>
  # mutate(Count2 = sum(Count)) |>
  # ungroup() |>
  # mutate(prop2 = Count2 / sum(Count2)) |>
  # group_by(Sample) |>
  # mutate(prop = Count / sum(Count)) |>
  mutate(z = value * prop) |>
  group_by(Kit, ind) |>
  summarize(z = sum(z)) |>
  # mutate(z = z / max(z)) |>
  # janitor::adorn_rounding(digits = 3) |>
  as.data.table() #|>
  # dcast(Kit ~ Individual, value.var = 'z')

plotdata |>
  ggplot(aes(x=ind, y=z, fill=Kit)) +
  geom_col() +
  facet_wrap(~ Kit, labeller = labeller(Kit = label_function), nrow=1) +
  scale_fill_manual(values = color_palette$kits, labels = label_function) +
  labs(x='Sample', y='Weighted distance') ->
figures[['pb_dist_table']]
write_plot_data(plotdata, here('figure_data/3p/replicate_correlation/rep_dist_barchart.txt'))
my_plot_save(image = figures[['pb_dist_table']], 
             path = here('figures/3p/replicate_correlation/rep_dist_barchart.svg'), 
             width = 10, height = 5)

# plotdata |> 
#   gt::gt() |>
#   tab_header('Composition-weighted pseudobulk replicate distance') |>
#   # cols_hide(c(low, high, avg)) |>
#   # cols_label(Kit = 'Kit', text = 'Sample distance') |>
#   data_color(columns = c("F1", 'F5'), method = 'numeric', palette = 'RdYlGn', reverse = TRUE) ->
#   figures[['pb_dist_table']]
# write_plot_data(plotdata, here('figure_data/3p/replicate_correlation/rep_dist_barchart.txt'))
# my_plot_save(image = figures[['pb_dist_table']], 
#              path = here('figures/3p/replicate_correlation/rep_dist_barchart.svg'), 
#              width = 2.75, height = 3.2)


