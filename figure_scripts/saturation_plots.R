library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations

source(here('config/color_palette.R'))

metadata <- read.csv(here('config/metadata.csv'), stringsAsFactors = TRUE)
sat_curves <- list()
for (file in list.files(here('data/saturation_curves'), pattern = '.csv$', full.names = TRUE)) {
  sample <- gsub('.+/B2_(.+).csv', '\\1', file)
  x <- read.csv(file)
  x$Sample <- sample
  if ('median_genes_per_cell' %in% colnames(x)) {
    lookup <- c('genes' = 'median_genes_per_cell', 'reads' = 'reads_per_cell')
    x <- rename(x, any_of(lookup)) |>
      select(genes, reads, Sample) 
  }
  # rename_with(.fn = function(col) {
    #   new_names <- col
    #   for (z in names(substring_map)) {
    #     if (grepl(pattern=z, x=col)) {
    #       new_names <- substring_map[[z]]
    #     }
    #   }
    #   return(new_names)
    # })
  sat_curves[[sample]] <- x
}
sat_curves <- data.table::rbindlist(sat_curves, use.names = TRUE)
sat_curves <- sat_curves %>%
  mutate(Sample = gsub('FBv4', 'FBv4', Sample)) %>%
  mutate(Sample = gsub('FBv5', 'FBv5', Sample)) %>%
  mutate(Sample = gsub('Parse', 'PA_V3', Sample)) %>%
  mutate(Sample = gsub('NextGEM', 'NextGEM3P', Sample)) %>%
  mutate(Sample = gsub('GEMX', 'GEMX3P', Sample)) %>%
  mutate(Sample = gsub('FL', 'Flex', Sample))

sat_curves <- merge(sat_curves, metadata, by='Sample')

ggplot(sat_curves, aes(x=as.numeric(reads), y=as.numeric(genes),
                       linetype = paste0(Individual, Replicate), color = Kit)) +
  scale_color_manual(values = color_palette) +
  geom_line() +
  theme_bw() +
  labs(x='Average reads per cell', y='Median genes per cell', linetype='Sample')
ggsave('saturation_curves.png', path = here('figures'), width = unit(7, 'in'), height = unit(5, 'in'))


# 
# samples <- c('Scale_F1A', 'FB_2_F5A', 'PA_V3_F5A', 
#              'Scale_F5A', 'FB_2_F5B', 'PA_V3_F1B')
# polygon_data <- sat_curves %>%
#   filter(Kit %in% c('Parse_V3', 'Fluent', 'Scale')) %>%
#   mutate(ID=paste0(Individual, Replicate)) %>%
#   filter(Sample %in% samples) %>%
#   mutate(Sample = factor(Sample, levels=samples)) %>%
#   arrange(Kit, Sample, case_when(
#     Sample %in% c('Scale_F1A', 'FB_2_F5A', 'PA_V3_F5A') ~ reads,
#     Sample %in% c('Scale_F5A', 'FB_2_F5B', 'PA_V3_F1B') ~ desc(reads),
#   ))
# 
# sat_curves %>%
#   filter(Kit %in% c('Parse_V3', 'Fluent', 'Scale')) %>%
#   mutate(ID=paste0(Individual, Replicate)) %>%
# ggplot(aes(x=reads, y=genes, color = Kit, linetype=ID)) +
#   scale_color_manual(values = color_palette) +
#   scale_fill_manual(values = color_palette) +
#   geom_line() +
#   geom_polygon(data = polygon_data, aes(x=reads, y=genes, group=Kit, fill=Kit), alpha=0.2, inherit.aes = FALSE)+
#   theme_bw() +
#   labs(x='Average reads per cell', y='Median genes per cell', linetype='Sample')
# ggsave('saturation_curves_blog.svg', path = here('figures/blog'), width = unit(7, 'in'), height = unit(5, 'in'))
