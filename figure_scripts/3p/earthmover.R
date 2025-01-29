library(transport)
source(here('figure_scripts/utils.R'))
library(waddR)

fig_objs <- readRDS(here('rds/3p/07_post_module_scoring.rds'))

kits <- unique(metadata_3p$Kit)
plotdata <- data.frame(row.names = kits, F1 = rep(NA, length(kits)), F5=rep(NA, length(kits)))
for (kit in kits){
  for (ind in c('F1', 'F5')) {
    umap1 <- subset(fig_objs[[kit]], orig.ident == paste0(kit, '_', ind, 'A')) 
    # umap1 <- umap1@reductions$UMAP_allgenes@cell.embeddings
    
    umap2 <- subset(fig_objs[[kit]], orig.ident == paste0(kit, '_', ind, 'B'))
    # umap2 <- umap2@reductions$UMAP_allgenes@cell.embeddings
    
    # if (nrow(umap1) > nrow(umap2)) {
    #   umap1 <- umap1[sample(1:nrow(umap1), nrow(umap2), replace = FALSE),]
    # } else if (nrow(umap2) > nrow(umap1)) {
    #   umap2 <- umap2[sample(1:nrow(umap2), nrow(umap1), replace = FALSE),]
    # }
    # 
    # umap1 <- pp(umap1)
    # umap2 <- pp(umap2)
    # emd_result <- transport::wasserstein(umap1, umap2, prob = FALSE)
    result <- wasserstein.sc(umap1, umap2, permnum=100, method='TS') ## increase permnum for real runs
    plotdata[kit, ind] <- emd_result
  }
}

# Calculate EMD
write_plot_data(plotdata, here('figure_data/3p/earthmover.txt'))

# emd_result <- transport::wasserstein(weights1, weights2, cost_matrix)
print(emd_result)