library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(yaml)       ## parses config yaml

# Load data ----
fig_objs <- readRDS(here('rds/wt/03-expression_summary.rds'))
kit_order <- read.table(here('config/kit_order.txt'))$V1
metadata <- read.csv(here('config/wt/metadata.csv')) |>
  mutate(kit = factor(kit, levels = kit_order))
  
# # Aggregate expression at sample level -----
# detected_genes_sample_level <- lapply(fig_objs, function(obj) {
#   genes <- obj@assays$RNA@features@.Data 
#   genes <- genes %>% 
#     as.data.frame() %>%
#     filter(counts == TRUE) %>% 
#     rownames()
#   detected_genes <- rowSums(obj@assays$RNA@layers$counts > 0) > 10 
#   detected_genes <- genes[detected_genes]
#   detected_genes
# })
# 
# # Define gene families
# human = useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl", verbose = TRUE)
# info <- getBM(attributes = c(#"ensembl_gene_id",
#   "external_gene_name",
#   "gene_biotype",
#   "transcript_length",
#   # "percentage_gene_gc_content",
#   "chromosome_name"),
#   mart = human) |>
#   filter(!is.na(external_gene_name) & external_gene_name != '') |>  
#   group_by(external_gene_name) |>
#   summarize(med_transcript_length = median(transcript_length),
#             # mean_transcript_length = mean(transcript_length), 
#             # ensembl_gene_id = paste0(ensembl_gene_id, collapse=', '), 
#             gene_biotype = paste0(unique(gene_biotype), collapse=', '))
# gene_families <- list(
#   'mitochondrial' = '^MT-',
#   'ribosomal' = '^RP[SL]',
#   'protein_coding_wout_gene_name' = c('^AC', '^LOC', 'orf', '^ENSG')
#   # 'ZNF' = '^ZNF',
#   # 'long intergenic non-coding' = '^LINC',
#   # 'histones' = c('^HIST', '^H\\d')
# )
# 
# 
# # Append group info ----
# find_group <- function(gene, patterns=gene_families) {
#   matched_group <- purrr::keep(patterns, function(pat_vec) {
#     any(sapply(pat_vec, function(pat) grepl(pat, gene)))
#   })
#   if (length(matched_group) > 0) {
#     return(paste0(names(matched_group), collapse = ','))
#   } else {
#     return("protein_coding_other")
#   }
# }
# find_group_vectorized <- Vectorize(find_group)
# gene_exp_summary <- lapply(gene_exp_summary, function(x) {
#   merge(x=x, y=info, by.x='row.names', by.y='external_gene_name', all.x=TRUE, all.y=FALSE) |>
#     dplyr::rename('gene' = 'Row.names') |>
#     mutate(shared = gene %in% shared_genes)
# })

# Map gene biotypes ----
fig_objs <- lapply(fig_objs, function(x) {
  x |>
    mutate(group = case_when(
      grepl(',', gene_biotype) ~ 'ambiguous',
      gene_biotype == 'protein_coding' ~ find_group_vectorized(gene),
      gene_biotype == 'lncRNA' ~ 'lncRNA',
      grepl('^M[Tt]', gene_biotype) ~ 'mitochondrial',
      gene_biotype == 'rRNA' ~ 'ribosomal',
      grepl('RNA$', gene_biotype) ~ 'other_RNA',
      grepl('pseudogene', gene_biotype) ~ 'pseudogene',
      # grepl('IG', gene_biotype) ~ 'IG_gene',
      # grepl('TR', gene_biotype) ~ 'TR_gene',
      # !is.na(gene_biotype) ~ gene_biotype,
      .default = 'Non protein coding'
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
  merge(metadata, by='Sample') |> 
  filter(shared == TRUE) |>
  mutate(group = factor(group, levels = c(
    'not shared', 
    'other',
    'ambiguous',
    'mitochondrial',
    'ribosomal',
    'lncRNA',
    'unidentified_protein_coding',
    'other_protein_coding'
  )))
plotdata |>
  ggplot(aes(x=paste0(Individual, Replicate), y=y, fill=group)) +
  geom_col(position='stack') +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  scale_fill_manual(values = c(
    'ambiguous' = 'grey',
    'lncRNA' = '#0072B2',
    'mitochondrial' = '#D55E00',
    'ribosomal' = '#F0E442',
    'unidentified_protein_coding' = 'darkblue',
    'other_protein_coding' = '#4DAF4A',
    'other' = '#CC79A7',
    'not shared' = 'black'
  )) + 
  theme(axis.text.x = element_text(angle=45, vjust=0.5),
        panel.spacing=unit(0, "lines")) + 
  facet_wrap(~ Kit, nrow=1, scales='free_x') +
  labs(y='Percent of total reads', x='Sample', fill='Feature biotype', 
       caption = 'Only showing reads assigned to genes found in all kits') ->
  figures_list[['Expression_barchart_shared']]
figures_list[['Expression_barchart_shared']]
ggsave(here('figures/wt/read_utilization/shared_genes.png'), 
       figures_list[['Expression_barchart_shared']],
       width = unit(9, 'in'), height = unit(4, 'in'), device = 'png')


## Unique genes -----

plotdata <- lapply(fig_objs, function(x) {
  x |>
    group_by(group, shared) |>
    summarize(y = sum(portion_of_total_reads))
}) |>
  data.table::rbindlist(idcol='Sample') |>
  merge(metadata, by='Sample') |> 
  filter(shared == FALSE) |>
  mutate(group = factor(group, levels = c(
    'not shared', 
    'other',
    'ambiguous',
    'mitochondrial',
    'ribosomal',
    'lncRNA',
    'unidentified_protein_coding',
    'other_protein_coding'
  )))
plotdata |>
  ggplot(aes(x=paste0(Individual, Replicate), y=y, fill=group)) +
  geom_col(position='stack') +
  theme_bw() +
  scale_fill_manual(values = c(
    'ambiguous' = 'grey',
    'lncRNA' = '#0072B2',
    'mitochondrial' = '#D55E00',
    'ribosomal' = '#F0E442',
    'unidentified_protein_coding' = 'darkblue',
    'other_protein_coding' = '#4DAF4A',
    'other' = '#CC79A7',
    'not shared' = 'black'
  )) + 
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle=45, vjust=0.5),
        panel.spacing=unit(0, "lines")) + 
  facet_wrap(~ Kit, nrow=1, scales='free_x') +
  labs(y='Percent of total reads', x='Sample', fill='Gene family', 
       caption = 'Only showing genes not found in all kits.\nThese genes are not necessarily unique to a kit, just not found in all kits') ->
  figures_list[['Expression_barchart_unique']]
figures_list[['Expression_barchart_unique']]
ggsave(here('figures/wt/read_utilization/unique_genes.png'), 
       figures_list[['Expression_barchart_unique']],
       width = unit(9, 'in'), height = unit(4, 'in'))

## Joint ----

plotdata <- lapply(fig_objs, function(x) {
  x |>
    group_by(group, shared) |>
    summarize(y = sum(portion_of_total_reads))
}) |>
  data.table::rbindlist(idcol='Sample') |>
  merge(metadata, by='Sample') |>
  mutate(group = ifelse(shared, group, 'not shared')) |>
  mutate(group = factor(group, levels = c(
    'not shared', 
    'other',
    'ambiguous',
    'mitochondrial',
    'ribosomal',
    'lncRNA',
    'unidentified_protein_coding',
    'other_protein_coding'
  )))
plotdata |>
  ggplot(aes(x=paste0(Individual, Replicate), y=y, fill=group)) +
  geom_col(position='stack') +
  scale_fill_manual(values = c(
    'ambiguous' = 'grey',
    'lncRNA' = '#0072B2',
    'mitochondrial' = '#D55E00',
    'ribosomal' = '#F0E442',
    'unidentified_protein_coding' = 'darkblue',
    'other_protein_coding' = '#4DAF4A',
    'other' = '#CC79A7',
    'not shared' = 'black'
  )) + 
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust=0.5),
        panel.spacing=unit(0, "lines")) + 
  facet_wrap(~ Kit, nrow=1, scales='free_x') +
  labs(y='Percent of total reads', x='Sample', fill='Feature biotype', 
       caption = 'Only showing reads assigned to genes found in all kits') ->
  figures_list[['Expression_barchart_all']]
figures_list[['Expression_barchart_all']]
ggsave(here('figures/wt/read_utilization/all_genes.png'), 
       figures_list[['Expression_barchart_all']],
       width = unit(9, 'in'), height = unit(4, 'in'), device = 'png')
