# Setup ----
source(here('figure_scripts/utils.R')) 
# library(transport)
# library(waddR)

# fig_objs <- readRDS(here('rds/3p/07_post_module_scoring.rds'))
emd_results_p <- read.csv(here('rds/3p/emd_bootstraps.txt'), sep = '\t')

# Prepare plotdata ----
plotdata <- 
  emd_results_p %>%
  group_by(Kit, Individual) %>%
  summarize(
    mean = mean(EMD),
    conf_low = t.test(EMD)$conf.int[1],
    conf_high = t.test(EMD)$conf.int[2],
    .groups = "drop"
  ) |> arrange(mean) |>
  mutate(Kit = factor(Kit, levels = kit_order_3p),
         ci = mean - conf_low) |>
  arrange(Kit, Individual)

# Calculate EMD
write_plot_data(plotdata, here('figure_data/3p/earthmover.txt'))

# Plot ----
plotdata <-
  emd_summary |>
  select(-c(conf_low, conf_high)) |>
  # data.table::as.data.table() |> 
  # knitr::kable() |>
  janitor::adorn_rounding(digits = 3) |>
  mutate(EMD = paste0(mean, ' ± ', ci)) |>
  group_by(Kit) |>
  mutate(tmp = mean(mean)) |>
  arrange(tmp, Individual) |>
  select(-c(tmp, ci)) |>
  as.data.table() #|>
# plotdata |> 
#   dcast(Kit ~ Individual)  |>
#   gt::gt() |>
#   data_color(
#     columns = c(F1, F5),
#     method='numeric',
#     # colors = scales::col_numeric(
#     #   palette = c("green", "yellow", "red"), # Customize the color scale
#     #   domain = NULL # Dynamic domain will be calculated later
#     # ),
#     apply_to = "fill",
#     fn = function(x) {
#       # Extract the mean dynamically from the cell values
#       mean_values <- as.numeric(sub(" ±.*", "", x))
#       
#       # Define the color scale using the extracted means
#       color_scale <- scales::col_numeric(
#         palette = c("#006837", "#FFFFBF", "#A50026"),
#         domain = c(0, .75),#range(mean_values, na.rm = TRUE)
#         na.color='red'
#       )
#       # Apply the color scale to the extracted means
#       color_scale(mean_values)
#     }
#   ) ->
#   figures[['emd_table']]

plotdata |>
  ggplot(aes(x = Individual, y=mean, fill=Kit)) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width=0.4, colour="black", alpha=0.9) +
  scale_fill_manual(values = color_palette$kits, labels = label_function) +
  facet_wrap(~ Kit, labeller = labeller( Kit = label_function), nrow=1) +
  labs(x='Individual', y='EMD') ->
  figures[['emd_barchart']]
figures[['emd_barchart']]

my_plot_save(figures[['emd_barchart']],
             here('figures/3p/emd.svg'),
             device ='svglite' ,
             width = 10, height = 5)
