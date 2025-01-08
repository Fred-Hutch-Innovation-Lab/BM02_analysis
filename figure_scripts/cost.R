library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)    
library(readxl)
library(data.table)

source(here('config/color_palette.R'))
source(here('config/kit_order.R'))
kit_order <- append(kit_order_3p, 'Parse_v2', which(kit_order_3p == 'Parse_v3') - 1)
metadata <- read.csv(here('config/3p/metadata.csv')) %>%
  mutate(Kit = factor(Kit, levels = kit_order))


cost_data_kit <- read_tsv(here('data/cost_v2.txt'))  %>%
    rowwise() %>%
    mutate(
      Extended = list(data.frame(
        samples = 1:Samples#,
      ))
    ) %>%
    unnest(cols = c(Extended)) |>
  mutate(cost_per_sample = Cost / samples,
         cells_per_sample = Cells / samples) |>
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
  ggplot(aes(x=samples, color=Kit, group=Kit)) +
  geom_line(aes(y=cost_per_sample), alpha=0.8, lty='solid') +
  # geom_step(aes(y=cells_per_sample), direction='hv', alpha=0.8, lty='dashed') +
  theme_bw() +
  # scale_linewidth_continuous(range=c(1,4), breaks = c(10000, 20000, 30000)) +
  scale_color_manual(values = unlist(color_palette$kits), labels = label_function) +
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