library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)    
library(readxl)

source(here('config/color_palette.R'))
kit_order <- read.table(here('config/3p/kit_order.txt'))$V1
kit_order <- append(kit_order, 'Parse_v2', which(kit_order == 'Parse_v3') - 1)
# metadata <- read.csv(here('config/metadata.csv')) %>%
  # mutate(Kit = factor(Kit, levels = kit_order))
cells_per_kit <- c(
  'Flex'=10000,
  'NextGEM3P'=10000,
  'GEMX3P'=20000,
  'Fluent_v4'=20000,
  'Fluent_V'=20000,
  'Parse_v3'=25000,
  'Parse_v2'=25000,
  'Scale'=31250
)


cost_data_kit <- read_xlsx(here('data/Benchmark Cost Analysis.xlsx'), sheet = 'plotting2',
                           .name_repair = 'minimal', col_names = TRUE,
                           skip = 0)
cost_data_kit

cost_data <- expand.grid(samples = 1:16, Kit = unique(cost_data_kit$Kit)) %>%
  left_join(cost_data_kit, by = "Kit", relationship = "many-to-many") %>%
  mutate(
    Rxns = as.numeric(Rxns), 
    num_items_needed = ceiling(samples / Rxns),
    total_cost_item = num_items_needed * Cost
  ) %>%
  group_by(Kit, samples) %>%
  summarize(total_cost = sum(total_cost_item, na.rm = TRUE), .groups = "drop") |>
  left_join(as.data.frame(cells_per_kit) |> rownames_to_column('Kit'), by='Kit') |>
  mutate(Kit = factor(Kit, levels=rev(kit_order)))


cost_data |>
  ggplot(aes(group=Kit, y=total_cost, x=samples, color=Kit, linewidth = cells_per_kit)) +
  geom_step(direction='hv', alpha=0.8, lty='solid') +
  theme_bw() +
  scale_linewidth_continuous(range=c(1,4), breaks = c(10000, 20000, 30000)) +
  scale_color_manual(values = unlist(color_palette$kits)) +
  labs(x='Samples', y='Total cost', color = 'Kit', linewidth='Cells per sample') + 
  guides(color = guide_legend(override.aes = list(linewidth = 4)))

ggsave(filename = 'cost_curves.png', path = here('figures'), device = 'png',
       width = unit(5, 'in'), height = unit(5, 'in'))
