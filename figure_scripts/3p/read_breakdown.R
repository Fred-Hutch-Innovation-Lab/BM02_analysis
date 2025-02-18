# Setup ----
source(here('figure_scripts/utils.R'))

# Load data ----
fig_objs <- readRDS(here('rds/3p/03-expression_summary.rds'))
gene_families <- list(
  'mitochondrial' = '^MT-',
  'ribosomal' = '^RP[SL]',
  'unidentified' = c('^AC', '^LOC', 'orf', '^ENSG')
)

# Plotting functions ----
# Prepare plotdata ----
find_group <- function(gene, patterns=gene_families) {
  matched_group <- purrr::keep(patterns, function(pat_vec) {
    any(sapply(pat_vec, function(pat) grepl(pat, gene)))
  }) |> unlist()
  if (length(matched_group) == 1) {
    return(matched_group)
  } else if (length(matched_group) > 1){
    return('unidentified')
  } else {
    return("identified_protein_coding")
  }
}
find_group_vectorized <- Vectorize(find_group)
fig_objs <- lapply(fig_objs, function(x) {
  x |>
    mutate(group = case_when(
      grepl(',', gene_biotype) ~ 'other',
      gene_biotype == 'protein_coding' ~ case_when(
        grepl('^MT-', gene) ~ 'mitochondrial',
        grepl('^RP[SL]', gene) ~ 'ribosomal',
        grepl('^AC', gene) ~ 'unidentified',
        grepl('^LOC', gene) ~ 'unidentified',
        grepl('orf', gene) ~ 'unidentified',
        grepl('^ENSG', gene) ~ 'unidentified',
        .default = 'identified_protein_coding'
      ), #find_group_vectorized(gene)
      gene_biotype == 'lncRNA' ~ 'lncRNA',
      grepl('^M[Tt]', gene_biotype) ~ 'mitochondrial',
      gene_biotype == 'rRNA' ~ 'ribosomal',
      grepl('RNA$', gene_biotype) ~ 'other',
      grepl('pseudogene', gene_biotype) ~ 'other',
      # grepl('IG', gene_biotype) ~ 'IG_gene',
      # grepl('TR', gene_biotype) ~ 'TR_gene',
      # !is.na(gene_biotype) ~ gene_biotype,
      .default = 'other'
    ))
})


# Plot ----
## Shared genes -----

plotdata <- lapply(fig_objs, function(x) {
  x |>
    group_by(group, shared) |>
    summarize(y = sum(portion_of_total_reads))
}) |>
  data.table::rbindlist(idcol='Sample') |>
  merge(metadata_3p, by='Sample') |> 
  filter(shared == TRUE) |>
  mutate(group = factor(group, levels = c(
    # 'not shared', 
    'other',
    'mitochondrial',
    'ribosomal',
    'lncRNA',
    'unidentified',
    'identified_protein_coding'
  )))
plotdata |>
  ggplot(aes(x=paste0(Individual, Replicate), y=y, fill=group)) +
  geom_col(position='stack') +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(
    'other' = '#CC79A7',
    'lncRNA' = '#0072B2',
    'mitochondrial' = '#D55E00',
    'ribosomal' = '#F0E442',
    'unidentified' = 'darkblue',
    'identified_protein_coding' = '#4DAF4A'
  )) + 
  theme(axis.text.x = element_text(angle=45, vjust=0.5),
        panel.spacing=unit(0, "lines")) + 
  facet_wrap(~ Kit, nrow=1, scales='free_x', labeller = labeller(Kit = label_function)) +
  labs(y='Percent of total reads', x='Sample', fill='Feature biotype'#, 
    # caption = paste(
    #   'Only showing reads assigned to genes found in all kits',
    #   '"Not shared" indicates a gene is not found in all kits.',
    #   'All other features are identified in all kits.',
    #   'Identified protein coding genes have symbols, non-identified genes only have location information or ENS IDs.',
    #   '"Other" indicates the gene could not be identified into a single biotype of these options.',
    # sep='\n')
    ) ->
  figures[['Expression_barchart_shared']]

my_plot_save(image = figures[['Expression_barchart_shared']], 
             path = here('figures/3p/read_utilization/shared_genes.svg'), 
             width = 12, height = 4)
plotdata |>
  mutate(sample = paste0(Individual, Replicate)) |>
  select(c(Kit, sample, group, y)) |>
  mutate(y = round(y * 100, 2)) |> 
  dcast(Kit + sample ~ group) |>
write_plot_data(here('figure_data/3p/read_utilization/shared_genes.txt'))
## Unique genes -----

plotdata <- lapply(fig_objs, function(x) {
  x |>
    group_by(group, shared) |>
    summarize(y = sum(portion_of_total_reads))
}) |>
  data.table::rbindlist(idcol='Sample') |>
  merge(metadata_3p, by='Sample') |> 
  filter(shared == FALSE) |>
  mutate(group = factor(group, levels = c(
    # 'not shared', 
    'other',
    'mitochondrial',
    'ribosomal',
    'lncRNA',
    'unidentified',
    'identified_protein_coding'
  )))
plotdata |>
  ggplot(aes(x=paste0(Individual, Replicate), y=y, fill=group)) +
  geom_col(position='stack') +
  scale_fill_manual(values = c(
    'other' = '#CC79A7',
    'lncRNA' = '#0072B2',
    'mitochondrial' = '#D55E00',
    'ribosomal' = '#F0E442',
    'unidentified' = 'darkblue',
    'identified_protein_coding' = '#4DAF4A'
  ), labels = c(
    'Other', 'Mitochondrial', 'Ribosomal', 'lncRNA', 'Unidentified', 'Identified protein coding'
  )) + 
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle=45, vjust=0.5),
        panel.spacing=unit(0, "lines")) + 
  facet_wrap(~ Kit, nrow=1, scales='free_x', labeller = labeller(Kit = label_function)) +
  labs(y='Percent of total reads', x='Sample', fill='Feature biotype'#, 
# caption = paste(
#   'Only showing genes not found in all kits.',
#   'These genes are not necessarily unique to a kit, just not found in all kits',
#   '"Not shared" indicates a gene is not found in all kits.',
#   'All other features are identified in all kits.',
#   'Identified protein coding genes have symbols, non-identified genes only have location information or ENS IDs.',
#   '"Other" indicates the gene could not be identified into a single biotype of these options.',
#   sep='\n')
) ->
  figures[['Expression_barchart_unique']]
my_plot_save(image = figures[['Expression_barchart_unique']], 
             path = here('figures/3p/read_utilization/unique_genes.svg'), 
             width = 12, height = 4)

plotdata |>
  mutate(sample = paste0(Individual, Replicate)) |>
  select(c(Kit, sample, group, y)) |>
  mutate(y = round(y * 100, 2)) |> 
  dcast(Kit + sample ~ group) |>
  write_plot_data(here('figure_data/3p/read_utilization/unique_genes.txt'))

## Joint ----

plotdata <- lapply(fig_objs, function(x) {
  x |>
    group_by(group, shared) |>
    summarize(y = sum(portion_of_total_reads))
}) |>
  data.table::rbindlist(idcol='Sample') |>
  merge(metadata_3p, by='Sample') |>
  mutate(group = ifelse(shared, group, 'not shared')) |>
  mutate(group = factor(group, levels = c(
    'not shared', 
    'other',
    # 'ambiguous',
    'mitochondrial',
    'ribosomal',
    'lncRNA',
    'unidentified',
    'identified_protein_coding'
  )))
plotdata |>
  ggplot(aes(x=paste0(Individual, Replicate), y=y, fill=group)) +
  geom_col(position='stack') +
  scale_fill_manual(values = c(
    # 'ambiguous' = 'grey',
    'other' = '#CC79A7',
    'lncRNA' = '#0072B2',
    'mitochondrial' = '#D55E00',
    'ribosomal' = '#F0E442',
    'unidentified' = 'darkblue',
    'identified_protein_coding' = '#4DAF4A',
    'not shared' = 'black'
  ), labels = c(
    'Not shared', 'Other', 'Mitochondrial', 'Ribosomal', 'lncRNA', 'Unidentified', 'Identified protein coding'
  )) + 
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle=45, vjust=0.5),
        panel.spacing=unit(0, "lines")) + 
  facet_wrap(~ Kit, nrow=1, scales='free_x', labeller = labeller(Kit = label_function)) +
  labs(y='Percent of total reads', x='Sample', fill='Feature biotype'#, 
       # caption = paste('"Not shared" indicates a gene is not found in all kits.',
       #                  'All other features are identified in all kits.',
       #                  'Identified protein coding genes have symbols.',
       #                  'Non-identified genes only have location information or ENS IDs.',
       #                 '"Other" indicates the gene could not be identified into a single biotype of these options.',
       #                  sep='\n')
       ) ->
  figures[['Expression_barchart_all']]
my_plot_save(image = figures[['Expression_barchart_all']], 
             path = here('figures/3p/read_utilization/all_genes.svg'), 
             width = 12, height = 4)

plotdata |>
  mutate(sample = paste0(Individual, Replicate)) |>
  select(c(Kit, sample, group, y)) |>
  mutate(y = round(y * 100, 5)) |> 
  group_by(Kit, sample, group) |>
  summarize(y = sum(y)) |>
  dcast(Kit + sample ~ group, value.var = 'y') |>
  write_plot_data(here('figure_data/3p/read_utilization/all_genes.txt'))
