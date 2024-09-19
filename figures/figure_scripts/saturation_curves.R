library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(reshape2)   ## DF manipulation
rbindlist <- data.table::rbindlist

set.seed(33)

metadata <- read.csv(here('config/metadata.csv'))
objs <- readRDS(here('rds/01-soup_channels.rds'))

samples <- names(objs)
reads = c(1000,1500,2000,2500,3000,4000,5000,6000, 7000, 8000, 9000, 10000)
downsampling_stats <- array(
  0,
  c(length(samples), length(reads)),
  list(objs = samples, reads)
)
for (i in 1:length(samples)) {
  for (j in 1:length(reads)) {
    read_count <- reads[j] * dim(objs[[samples[i]]]$toc)[2]
    prop <- round(read_count / sum(objs[[samples[i]]]$toc), 3)
    if (prop < 1) {
      counts <- scuttle::downsampleMatrix(objs[[samples[i]]]$toc, prop, bycol = FALSE)
      downsampling_stats[i,j] <- median(colSums(counts > 0))
    } else {
      downsampling_stats[i,j] <- NA
    }
    
  }
}
downsampling_stats <- reshape2::melt(downsampling_stats, varnames = c('sample', 'reads'))
## Add current seq depth
for (i in 1:length(samples)) {
  read_count <- round(sum(objs[[samples[i]]]$toc) / dim(objs[[samples[i]]]$toc)[2])
  downsampling_stats <- downsampling_stats %>% 
    add_row(sample=samples[i], reads=read_count, value=median(colSums(objs[[samples[i]]]$toc > 0)))
}
downsampling_stats <- filter(downsampling_stats, !is.na(value))


downsampling_stats %>%
  mutate(kit = gsub('(.+)_.*$', '\\1', sample),
         individual = gsub('.*_(F\\d).?$', '\\1', sample)) %>%
  ggplot(aes(x=reads, y=value, group=sample, color=individual)) +
  # geom_boxplot() +
  geom_point(alpha=0.5) +
  geom_line(aes(group=sample), alpha=0.5, lty='solid', inherit.aes = TRUE) +
  theme_bw() +
  facet_wrap(~ kit) +
  labs(x = 'Average transcripts per cell', y='Median genes per cell')
