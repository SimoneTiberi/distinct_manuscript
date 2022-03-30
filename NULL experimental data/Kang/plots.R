rm(list = ls())

setwd("~/Desktop/distinct project/EXPERIMENTAL data/Kang (muscData)/")
load("results/res_cell_clust_7Samples.RData")

RES_ALL = list()
for(i in 1:3){
  RES_ALL[[i]] = c( res_distinct[[i]],
                    RES_PB[[i]] )
  names_distinct = paste0("distinct.", names(res_distinct[[i]]))
  names_PB = names(RES_PB[[i]])
  
  # remove PB middle part name:
  names_PB = sapply(names_PB, function(x){
    a = strsplit(x, ".", fixed = TRUE)[[1]]
    paste0(a[1], ".", a[3])
  })
  
  # replace "scalecpm" with "cpm":
  names_PB = sub("scalecpm", "cpm", names_PB)
  
  names(RES_ALL[[i]]) = c(names_distinct, names_PB)
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
                                 if(is.null(res)){
                                   return(res)
                                 }
                                 res$clust_gene = paste(res$cluster_id, res$gene, sep = "_")
                                 sel = res$clust_gene %in% gene_clust_keep
                                 res[sel,]
                               })
                      })
sapply(RES_ALL_filt[[1]], nrow)
sapply(RES_ALL_filt[[2]], nrow)
sapply(RES_ALL_filt[[3]], nrow)

# same length for all methods
length(unique(sapply(RES_ALL_filt[[1]], nrow)))
length(unique(sapply(RES_ALL_filt[[2]], nrow)))
length(unique(sapply(RES_ALL_filt[[3]], nrow)))

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

# match colours to methods:
methods = sort(unique(p_vals_filt$method))

source("~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/all_methods.R")

match = match( methods,  methods_names )
colours = c(all_colours[match], "white")

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
       path = "~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/NULL experimental data/Figures",
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
% latex table generated in R 4.1.1 by xtable 1.8-4 package
% Thu Mar 17 16:49:37 2022
\begin{table}[ht]
\centering
\begin{tabular}{rrrr}
\hline
& 0.1 & 0.05 & 0.01 \\ 
\hline
distinct.basics & 0.13 & 0.08 & 0.02 \\ 
distinct.cpm & 0.13 & 0.08 & 0.02 \\ 
distinct.linnorm & 0.13 & 0.08 & 0.02 \\ 
distinct.logcounts & 0.13 & 0.08 & 0.03 \\ 
distinct.vstresiduals & 0.14 & 0.09 & 0.03 \\ 
edgeR.basics & 0.00 & 0.00 & 0.00 \\ 
edgeR.counts & 0.09 & 0.05 & 0.01 \\ 
edgeR.cpm & 0.05 & 0.02 & 0.00 \\ 
edgeR.linnorm & 0.08 & 0.04 & 0.01 \\ 
limma-trend.basics & 0.14 & 0.08 & 0.02 \\ 
limma-trend.cpm & 0.20 & 0.11 & 0.03 \\ 
limma-trend.linnorm & 0.15 & 0.08 & 0.02 \\ 
limma-trend.logcounts & 0.12 & 0.07 & 0.02 \\ 
limma-trend.vstresiduals & 0.13 & 0.07 & 0.01 \\ 
limma-voom.counts & 0.08 & 0.04 & 0.01 \\ 
\hline
\end{tabular}
\end{table}

FPs_adj.loc = t(sapply( split(p_vals_filt$p_adj.loc, p_vals_filt$method), function(x){
  c(mean(x< 0.1), mean(x< 0.05), mean(x< 0.01) )
}))
colnames(FPs_adj.loc) = c(0.1, 0.05, 0.01)

library(xtable)
xtable(FPs_adj.loc)
% latex table generated in R 4.1.1 by xtable 1.8-4 package
% Thu Mar 17 16:50:00 2022
\begin{table}[ht]
\centering
\begin{tabular}{rrrr}
\hline
& 0.1 & 0.05 & 0.01 \\ 
\hline
distinct.basics & 0.01 & 0.01 & 0.00 \\ 
distinct.cpm & 0.01 & 0.01 & 0.00 \\ 
distinct.linnorm & 0.01 & 0.01 & 0.00 \\ 
distinct.logcounts & 0.01 & 0.01 & 0.00 \\ 
distinct.vstresiduals & 0.01 & 0.01 & 0.00 \\ 
edgeR.basics & 0.00 & 0.00 & 0.00 \\ 
edgeR.counts & 0.00 & 0.00 & 0.00 \\ 
edgeR.cpm & 0.00 & 0.00 & 0.00 \\ 
edgeR.linnorm & 0.00 & 0.00 & 0.00 \\ 
limma-trend.basics & 0.00 & 0.00 & 0.00 \\ 
limma-trend.cpm & 0.02 & 0.01 & 0.00 \\ 
limma-trend.linnorm & 0.00 & 0.00 & 0.00 \\ 
limma-trend.logcounts & 0.00 & 0.00 & 0.00 \\ 
limma-trend.vstresiduals & 0.00 & 0.00 & 0.00 \\ 
limma-voom.counts & 0.00 & 0.00 & 0.00 \\ 
\hline
\end{tabular}
\end{table}

FPs_adj.glb = t(sapply( split(p_vals_filt$p_adj.glb, p_vals_filt$method), function(x){
  c(mean(x< 0.1), mean(x< 0.05), mean(x< 0.01) )
}))
colnames(FPs_adj.glb) = c(0.1, 0.05, 0.01)

library(xtable)
xtable(FPs_adj.glb)
% latex table generated in R 4.1.1 by xtable 1.8-4 package
% Thu Mar 17 16:50:10 2022
\begin{table}[ht]
\centering
\begin{tabular}{rrrr}
\hline
& 0.1 & 0.05 & 0.01 \\ 
\hline
distinct.basics & 0.01 & 0.00 & 0.00 \\ 
distinct.cpm & 0.01 & 0.00 & 0.00 \\ 
distinct.linnorm & 0.01 & 0.00 & 0.00 \\ 
distinct.logcounts & 0.01 & 0.01 & 0.00 \\ 
distinct.vstresiduals & 0.01 & 0.00 & 0.00 \\ 
edgeR.basics & 0.00 & 0.00 & 0.00 \\ 
edgeR.counts & 0.00 & 0.00 & 0.00 \\ 
edgeR.cpm & 0.00 & 0.00 & 0.00 \\ 
edgeR.linnorm & 0.00 & 0.00 & 0.00 \\ 
limma-trend.basics & 0.00 & 0.00 & 0.00 \\ 
limma-trend.cpm & 0.00 & 0.00 & 0.00 \\ 
limma-trend.linnorm & 0.00 & 0.00 & 0.00 \\ 
limma-trend.logcounts & 0.00 & 0.00 & 0.00 \\ 
limma-trend.vstresiduals & 0.00 & 0.00 & 0.00 \\ 
limma-voom.counts & 0.00 & 0.00 & 0.00 \\ 
\hline
\end{tabular}
\end{table}

all_FP = data.frame(cbind(FPs, FPs_adj.glb, FPs_adj.loc))
colnames(all_FP) = rep(c(0.1, 0.05, 0.01),  3)
library(xtable)
xtable(all_FP)
# Supplementary Table 4


\begin{table}[ht]
\centering
\begin{tabular}{rrrrrrrrrr}
\hline
& 0.1 & 0.05 & 0.01 & 0.1 & 0.05 & 0.01 & 0.1 & 0.05 & 0.01 \\ 
\hline
distinct.basics & 0.13 & 0.08 & 0.02 & 0.01 & 0.00 & 0.00 & 0.01 & 0.01 & 0.00 \\ 
distinct.cpm & 0.13 & 0.08 & 0.02 & 0.01 & 0.00 & 0.00 & 0.01 & 0.01 & 0.00 \\ 
distinct.linnorm & 0.13 & 0.08 & 0.02 & 0.01 & 0.00 & 0.00 & 0.01 & 0.01 & 0.00 \\ 
distinct.logcounts & 0.13 & 0.08 & 0.03 & 0.01 & 0.01 & 0.00 & 0.01 & 0.01 & 0.00 \\ 
distinct.vstresiduals & 0.14 & 0.09 & 0.03 & 0.01 & 0.00 & 0.00 & 0.01 & 0.01 & 0.00 \\ 
edgeR.basics & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 \\ 
edgeR.counts & 0.09 & 0.05 & 0.01 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 \\ 
edgeR.cpm & 0.05 & 0.02 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 \\ 
edgeR.linnorm & 0.08 & 0.04 & 0.01 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 \\ 
limma-trend.basics & 0.14 & 0.08 & 0.02 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 \\ 
limma-trend.cpm & 0.20 & 0.11 & 0.03 & 0.00 & 0.00 & 0.00 & 0.02 & 0.01 & 0.00 \\ 
limma-trend.linnorm & 0.15 & 0.08 & 0.02 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 \\ 
limma-trend.logcounts & 0.12 & 0.07 & 0.02 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 \\ 
limma-trend.vstresiduals & 0.13 & 0.07 & 0.01 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 \\ 
limma-voom.counts & 0.08 & 0.04 & 0.01 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 \\ 
\hline
\end{tabular}
\end{table}
