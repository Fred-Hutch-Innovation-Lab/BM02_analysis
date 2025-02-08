# Setup ----
source(here('figure_scripts/utils.R')) 
library(readxl)

# Load data ----
cost_data_kit <- read_csv(here('data/cost.txt'), sep='\t')  |>
  mutate(cost = readr::parse_number(cost)) |>
    rowwise() %>%
    mutate(
      Extended = list(data.frame(
        samples = 1:max_samples#,
      ))
    ) %>%
    unnest(cols = c(Extended)) |>
  mutate(cost_per_sample = ifelse(all_in, cost / samples, cost / max_samples),
         cells_per_sample = ifelse(all_in, max_cells / samples, max_cells / max_samples)) #|>
  select(-c(Samples, Cells, Cost)) |>
  as.data.table()
# 
# cost_data_5p <- cost_data_kit |>
#   filter(kit %in% c(''))

# cost_data <- expand.grid(samples = 1:16, Kit = unique(cost_data_kit$Kit)) %>%
#   left_join(cost_data_kit, by = "Kit", relationship = "many-to-many") %>%
#   mutate(
#     Rxns = as.numeric(Rxns), 
#     num_items_needed = ceiling(samples / Rxns),
#     total_cost_item = num_items_needed * Cost
#   ) %>%
#   group_by(Kit, samples) %>%
#   summarize(total_cost = sum(total_cost_item, na.rm = TRUE), .groups = "drop") |>
#   left_join(as.data.frame(cells_per_kit) |> rownames_to_column('Kit'), by='Kit') |>
#   mutate(Kit = factor(Kit, levels=rev(kit_order)))


cost_data_kit |>
  filter(!kit_family %in% c('Parse_WT', 'Parse_TCR', 'Scale')) |>
  ggplot(aes(x=samples, color=kit, group=kit)) +
  geom_line(aes(y=cost_per_sample, linewidth=cells_per_sample), alpha=0.8, lty='solid') +
  # geom_step(aes(y=cells_per_sample), direction='hv', alpha=0.8, lty='dashed') +
  theme_bw() +
  facet_wrap(~ kit_family, nrow=1, scales='fixed') +
  # scale_linewidth_continuous(range=c(1,4), breaks = c(10000, 20000, 30000)) +
  # scale_color_manual(values = unlist(color_palette$kits), labels = label_function) +
  labs(x='Samples', y='Cost per sample', color = 'Kit') #+ 
  # guides(color = guide_legend(override.aes = list(linewidth = 4)))

ggsave(filename = 'cost_curves.png', path = here('figures'), device = 'png',
       width = unit(5, 'in'), height = unit(5, 'in'))

cost_data_kit |>
  ggplot(aes(x=samples, color=Kit, group=Kit)) +
  geom_line(aes(y=cells_per_sample), alpha=0.8, lty='solid') +
  # geom_step(aes(y=cells_per_sample), direction='hv', alpha=0.8, lty='dashed') +
  theme_bw() +
  # scale_linewidth_continuous(range=c(1,4), breaks = c(10000, 20000, 30000)) +
  scale_color_manual(values = unlist(color_palette$kits), labels = label_function) +
  labs(x='Samples', y='Cells per sample', color = 'Kit') #+ 
# guides(color = guide_legend(override.aes = list(linewidth = 4)))