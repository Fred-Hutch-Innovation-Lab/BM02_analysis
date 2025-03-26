# Setup ----
source(here('figure_scripts/utils.R'))
library(fgsea)
library(msigdbr)

# Load data ----
gsea_results <- readRDS(here('rds/3p/gsea_results.rds'))

# Prepare plotdata ----
pathways_of_interest <- c('TCR', 'MONOCYTE', 'BCR') #'TCRA','TCYTOTOXIC', 'CTL', 'TID', 'LYM', 'LYMPHOCYTE',, 'BLYMPHOCYTE'
plotdata <- rbindlist(gsea_results, idcol='contrast') |>
  mutate(kit = case_when(
    grepl('Flex', contrast) ~ 'Flex',
    grepl('NextGEM3P', contrast) ~ 'NextGEM3P',
    grepl('GEMX3P', contrast) ~ 'GEMX3P',
    grepl('Fluent_v4', contrast) ~ 'Fluent_v4',
    grepl('Fluent_V', contrast) ~ 'Fluent_V',
    grepl('Parse_v3', contrast) ~ 'Parse_v3',
    grepl('Scale', contrast) ~ 'Scale',
  ), comp = case_when(
    grepl('B.vs.T', contrast) ~ 'B cells vs T cells',
    grepl('T.vs.M', contrast) ~ 'T cells vs Monocytes',
    grepl('M.vs.B', contrast) ~ 'Monocytes vs B cells',
    grepl('B.vs.all', contrast) ~ 'B',
    grepl('T.vs.all', contrast) ~ 'T',
    grepl('M.vs.all', contrast) ~ 'M',
    .default = NA
  )) |>
  filter(!is.na(comp), name %in% pathways_of_interest) |>
  mutate(name = factor(name, levels = pathways_of_interest),
         sig = case_when(
           # padj < 0.0005 ~ '***',
           # padj < 0.005 ~ '**',
           padj < 0.05 ~ '*',
           .default = NA
         )) |>
  mutate(kit = factor(kit, kit_order_3p)) |>
  filter(comp %in% c('B cells vs T cells', 'T cells vs Monocytes', 'Monocytes vs B cells'))

# Plot ----
ggplot(plotdata, aes(x=kit, color=NES, size=-log(padj), label = sig, y = name)) +
  geom_point() +
  facet_wrap(~ comp, nrow=1, labeller = label_wrap_gen()) +
  scale_x_discrete(labels = label_function()) + 
  scale_color_gradient2(low="blue",
                       mid="white",
                       high="red",
                       midpoint=0,
                       breaks=c(-2,-1,0,1,2),
                       limits=c(min(plotdata$NES,-2),
                              max(plotdata$NES,2))) +
  theme(panel.grid.major = element_line(colour = "grey92"),
        panel.grid.minor = element_line(colour = "grey92"),
        panel.grid.major.y = element_line(colour = "grey92"),
        panel.grid.minor.y = element_line(colour = "grey92")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(), 
        # panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank()) +
  labs(x="Kit", y='BioCarta Pathway',
       color = "Normalized\nenrichment\nscore",
       size="Adjusted p-value",
       caption = '* indicates adjusted p-value < 0.05') +
  geom_text(na.rm = TRUE, color = 'white', size = 5) + 
  scale_size_area(max_size = 10, breaks=c(1,2,3), labels = c(1, 0.1, 0.01)) ->
  figures[['gsea_dotplot']]

# Save plot
my_plot_save(figures[['gsea_dotplot']], 
             path = here('figures/3p/gsea_dotplot.svg'),
             width = 13.5,
             height = 4.5)


write_plot_data(plotdata, file = here('figure_data/3p/gsea_dotplot.txt')) 
