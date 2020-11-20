rm(list = ls())

setwd("~/Desktop/distinct project/EXPERIMENTAL data/Kang (muscData)/")
load("results/res_cell_clust_7Samples.RData")

RES_ALL = list()
for(i in 1:3){
  RES_ALL[[i]] = c( res_distinct[[i]],
                    RES_PB[[i]] )
  names(RES_ALL[[i]]) = c("distinct.cpm", "distinct.log2-cpm", "distinct.vstresiduals",
                          "edgeR.counts", "limma-voom.counts", "edgeR.cpm",           
                          "limma-trend.log2-cpm", "limma-trend.vstresiduals")
}

str(RES_ALL[[1]])
str(RES_ALL[[2]])
str(RES_ALL[[3]])

##### FILTER RESULTS (AT LEAST 10 NON-ZERO GENES TO ANALYZE A CLUSTER-GENE):
gene_clust_keep = c()

clust_levels = levels(clust)

sce_counts = as.matrix(sce_counts)
rownames(sce_counts) = seq_len(nrow(sce_counts))

for(i in 1:length(clust_levels)){
  sel = clust == clust_levels[i]
  sel[is.na(sel)] = FALSE # if clust is NA, set to FALSE (not clust_levels[i])
  
  # select genes with at least 20 counts in each cluster:
  sel_genes = rowSums( sce_counts[, sel] > 0 ) >= 20 
  
  genes_to_keep = rownames(sce_counts)[sel_genes]
  
  # if at least 1 gene:
  if(!is.null(genes_to_keep)){
    gene_clust_keep = c(gene_clust_keep, paste(clust_levels[i], genes_to_keep, sep = "_") )
  }
}
head(gene_clust_keep); tail(gene_clust_keep)

# FILTER RESULTS (REMOVE LOWLY ABUNDANT CLUSTERS-GENES):
RES_ALL_filt = lapply(RES_ALL, 
                      function(RES){
                        lapply(RES, 
                               function(res){
                                 res$clust_gene = paste(res$cluster_id, res$gene, sep = "_")
                                 sel = res$clust_gene %in% gene_clust_keep
                                 res[sel,]
                               })
                      })
sapply(RES_ALL_filt[[1]], nrow)
sapply(RES_ALL_filt[[2]], nrow)
sapply(RES_ALL_filt[[3]], nrow)

sapply(RES_ALL[[1]], nrow)
sapply(RES_ALL[[2]], nrow)
sapply(RES_ALL[[3]], nrow)

###################################################################################################
# Make nice null plot:
###################################################################################################
suppressMessages({
  library(data.table)
  library(dplyr)
  library(iCOBRA)
  library(ggplot2)
  library(purrr)
  library(reshape2)
})

p_vals_filt = do.call(rbind, 
                      do.call(rbind, args = lapply(1:length(RES_ALL_filt), function(X){
                        lapply(1:length(RES_ALL_filt[[X]]),
                               function(i){
                                 data.frame(p_val = RES_ALL_filt[[X]][[i]]$p_val, 
                                            p_adj.loc = RES_ALL_filt[[X]][[i]]$p_adj.loc,
                                            p_adj.glb = RES_ALL_filt[[X]][[i]]$p_adj.glb, 
                                            replicate = X, 
                                            method = names(RES_ALL_filt[[X]])[i])
                               })
                      })
                      ))
head(p_vals_filt)
tail(p_vals_filt)
dim(p_vals_filt)

levels(factor(p_vals_filt$method))


.prettify <- function(theme = NULL, ...) {
  if (is.null(theme)) theme <- "classic"
  base <- paste0("theme_", theme)
  base <- getFromNamespace(base, "ggplot2")
  base(base_size = 8) + theme(
    aspect.ratio = 1,
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.2, color = "lightgrey"),
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text = element_text(color = "black"),
    legend.key.size = unit(2, "mm"),
    strip.background = element_rect(fill = NA),
    plot.margin = unit(rep(1, 4), "mm"),
    panel.spacing = unit(0, "mm"),
    legend.margin = margin(0,0,1,0,"mm"),
    ...)}


p_vals_filt$replicate = as.factor(p_vals_filt$replicate)

p <- ggplot(data = p_vals_filt, aes(x = p_val, y = ..ndensity.., 
                                    col = method, fill = method, lty = replicate),
            colour = colours) +
  facet_wrap(~ method, ncol = 4) + 
  geom_density(adjust = 1, size = 0.3, alpha = 0.1) +
  scale_alpha_manual(values = c("TRUE" = 0.1, "FALSE" = 0.4)) +
  guides(col = FALSE,
         lty = guide_legend(ncol = 1, order = 2),
         fill = guide_legend(ncol = 1, order = 1,
                             override.aes = list(alpha = 1, col = NA))) +
  scale_x_continuous("p-value", breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
  scale_y_continuous("normalized density", breaks = c(0, 1), expand = c(0, 0.1)) +
  .prettify("bw") + theme(aspect.ratio = 1/2,
                          # legend.box.just = "left",
                          panel.grid = element_blank(),
                          panel.spacing = unit(1, "mm"),
                          panel.border = element_rect(color = "grey"),
                          strip.text = element_text(size = 5),
                          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# match colours to methods:
methods = unique(p_vals_filt$method)

library(RColorBrewer);
all_colours = c(
  brewer.pal(4, "Reds")[4:2],    # distinct cpm, logcounts, vstresid
  brewer.pal(3, "BuPu")[3:2],    # edgeR 2 methods
  brewer.pal(4, "Greens")[4:2]   # limma-trend + # limma-voom 2 methods
)

all_methods = c(
  "distinct.cpm", # distinct.cpm25FALSE, or distinct.cpm25FALSE1000,
  "distinct.log2-cpm", # distinct.logcounts25FALSE, or distinct.logcounts25FALSE1000,
  "distinct.vstresiduals", # distinct.vstresiduals25FALSE, or distinct.vstresiduals25FALSE1000,
  "edgeR.counts",
  "edgeR.cpm",
  "limma-voom.counts",
  "limma-trend.log2-cpm",
  "limma-trend.vstresiduals"
)

# all_methods = sort(all_methods)
match = match( methods,  all_methods )

colours = c(all_colours[match])


p <- ggplot(data = p_vals_filt, aes(x = p_val, y = ..ndensity.., 
                                    col = method, fill = method, lty = replicate),
            colour = colours) +
  facet_wrap(~ method, ncol = 4) + 
  geom_density(adjust = 1, size = 0.3, alpha = 0.1) + 
  scale_alpha_manual(values = c("TRUE" = 0.1, "FALSE" = 0.4)) +
  guides(col = FALSE,
         lty = guide_legend(ncol = 1, order = 2),
         fill = guide_legend(ncol = 2, order = 1,
                             override.aes = list(alpha = 1, col = NA))) +
  scale_x_continuous("p-value", breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
  scale_y_continuous("density", breaks = c(0, 1), expand = c(0, 0.1)) +
  .prettify("bw") + 
  theme(aspect.ratio = 1/2,
        #legend.box.just = "left",
        panel.grid = element_blank(),
        panel.spacing = unit(1, "mm"),
        panel.border = element_rect(color = "grey"),
        strip.text = element_text(size = 4),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(1.5)),
        axis.text.y=element_text(size=rel(1.5)),
        axis.title.y = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(1.5)),
        axis.title.x.top = element_text(size=rel(1.5)),
        legend.title=element_text(size=rel(1.5)),
        legend.text=element_text(size=rel(1.5)),
        legend.key.width=unit(0.4, "cm"),
        legend.position="bottom",
        legend.direction = "horizontal",
        legend.box="horizontal",
        strip.text.x = element_text(size = rel(2)),
        legend.margin=margin()) +
  scale_colour_manual(values = colours) +
  scale_fill_manual(values = colours)
p

p_Kang = p
save(p_Kang, file = "p_Kang.RData")

ggsave(filename = paste0("null-Kang.pdf"),
       plot = last_plot(),
       device = "pdf",
       path = "~/Desktop/distinct project/distinct Article/v1/images/NULL",
       width = 8,
       height = 4,
       units = "in",
       dpi = 300,
       limitsize = TRUE)


###################################################################################################
# FP tables:
###################################################################################################
FPs = t(sapply( split(p_vals_filt$p_val, p_vals_filt$method), function(x){
  c(mean(x< 0.1), mean(x< 0.05), mean(x< 0.01) )
}))
colnames(FPs) = c(0.1, 0.05, 0.01)

library(xtable)
xtable(FPs)

FPs_adj.loc = t(sapply( split(p_vals_filt$p_adj.loc, p_vals_filt$method), function(x){
  c(mean(x< 0.1), mean(x< 0.05), mean(x< 0.01) )
}))
colnames(FPs_adj.loc) = c(0.1, 0.05, 0.01)

library(xtable)
xtable(FPs_adj.loc)


FPs_adj.glb = t(sapply( split(p_vals_filt$p_adj.glb, p_vals_filt$method), function(x){
  c(mean(x< 0.1), mean(x< 0.05), mean(x< 0.01) )
}))
colnames(FPs_adj.glb) = c(0.1, 0.05, 0.01)

library(xtable)
xtable(FPs_adj.glb)



all_FP = data.frame(cbind(FPs, FPs_adj.glb, FPs_adj.loc))
colnames(all_FP) = rep(c(0.1, 0.05, 0.01),  3)
library(xtable)
xtable(all_FP)
# Supplementary Table 4
