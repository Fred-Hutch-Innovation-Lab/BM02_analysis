library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations

source(here('config/color_palette.R'))

metadata <- read.csv(here('config/metadata.csv'), stringsAsFactors = TRUE)
sat_curves <- list()
for (file in list.files(here('data/saturation_curves'), pattern = '.csv$', full.names = TRUE)) {
  sample <- gsub('.+/B2_(.+).csv', '\\1', file)
  x <- read.csv(file)
  x$Sample <- sample
  sat_curves[[sample]] <- x
}
sat_curves <- data.table::rbindlist(sat_curves)
sat_curves <- sat_curves %>%
  mutate(Sample = gsub('FB', 'FB_2', Sample)) %>%
  mutate(Sample = gsub('Parse', 'PA_V3', Sample)) %>%
  mutate(Sample = gsub('NextGEM', 'NextGEM3P', Sample)) %>%
  mutate(Sample = gsub('GEMX', 'GEMX3P', Sample)) %>%
  mutate(Sample = gsub('FL', 'Flex', Sample))

sat_curves <- merge(sat_curves, metadata, by='Sample')

ggplot(sat_curves, aes(x=reads, y=genes, linetype = paste0(Individual, Replicate), color = Kit)) +
  scale_color_manual(values = color_palette) +
  geom_line() +
  theme_bw() +
  labs(x='Average reads per cell', y='Median genes per cell', linetype='Sample')
ggsave('saturation_curves.png', path = here('figures'), width = unit(7, 'in'), height = unit(5, 'in'))
