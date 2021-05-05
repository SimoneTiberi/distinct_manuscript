rm(list = ls())
setwd("~/Desktop/distinct project/SIMULATED data/KANG Real data/")

saving = TRUE

local = TRUE # locally vs globally adjusted p-values

res_names = list.files(paste0("results/output"), full.names =TRUE)
res_names
file.exists(res_names)

methods = c(
  "distinct.cpm",
  "distinct.logcounts",
  "distinct.vstresiduals",
  "edgeR.sum.counts",
  "edgeR.sum.scalecpm",
  "limma-voom.sum.counts",
  "limma-trend.mean.logcounts",
  "limma-trend.mean.vstresiduals"
)

# method names:
methods_names = c(
  "distinct.cpm", # distinct.cpm25FALSE, or distinct.cpm25FALSE1000,
  "distinct.log2-cpm", # distinct.logcounts25FALSE, or distinct.logcounts25FALSE1000,
  "distinct.vstresiduals", # distinct.vstresiduals25FALSE, or distinct.vstresiduals25FALSE1000,
  "edgeR.counts",
  "edgeR.cpm",
  "limma-voom.counts",
  "limma-trend.log2-cpm",
  "limma-trend.vstresiduals"
)

# FILTER LOWLY ABUNDANT gene-cluster combinations:
# re-load only if not loaded yet:

library(SingleCellExperiment)
sce = readRDS("raw_data/sce0_kang.rds")
sce <- sce[, sce$ind %in% c(101, 107, 1015, 1016, 1244, 1256, 1488)] # keep control samples only

sce <- sce[, sce$multiplets == "singlet"] # remove multiplets
sce <- sce[, !is.na(sce$cell)] # remove unassigned cells

sce$cluster_id = sce$cell
sce$sample_id = paste(sce$stim, sce$ind)
sce$group_id = sce$stim

sce <- sce[rowSums(counts(sce) > 1) > 20, ]

sce <- muscat::prepSCE(sce, "cluster_id", "sample_id", "group_id", TRUE) # prep. SCE for `muscat`
sce

# compute gene-cluster combinations with at least 10 non-zero cells.
gene_clust_keep = c()
clust_levels = levels(sce$cluster_id)
for(i in 1:length(clust_levels)){
  sel = sce$cluster_id == clust_levels[i]
  sel[is.na(sel)] = FALSE # if clust is NA, set to FALSE (not clust_levels[i])
  
  # select genes with at least 20 counts in the cluster:
  sel_genes = rowSums( assays(sce)$counts[,sel] > 0 ) >= 20
  
  genes_to_keep = rownames(assays(sce)$counts)[sel_genes]
  
  gene_clust_keep = c(gene_clust_keep, paste(clust_levels[i], genes_to_keep, sep = "___") )
}
head(gene_clust_keep); tail(gene_clust_keep)

RES = list()
for(j in seq_along(methods)){
  sel_method = grep(methods[j], res_names)
  if(length(sel_method) == 0){
    # is method not found:
    RES[[j]] = NULL
  }else{
    res = readRDS(res_names[sel_method])
    # PB methods crash: no metadata
    # <simpleError in .check_pbs(pb, check_by = TRUE): !is.null(ei <- metadata(pbs)$experiment_info) is not TRUE>
    
    res$clust_gene = paste(res$cluster_id, res$gene, sep = "___")
    sel = res$clust_gene %in% gene_clust_keep
    res = res[sel,]
    
    # MATCH RESULTS with gene-cluster names:
    matching = match(res$clust_gene, gene_clust_keep)
    head(matching); head(gene_clust_keep); head(res$clust_gene[matching])
    tail(matching); tail(gene_clust_keep); tail(res$clust_gene[matching])
    
    RES[[j]] = data.frame( p_val = res$p_val, 
                           p_adj.loc = res$p_adj.loc,
                           p_adj.glb = res$p_adj.glb,
                           gene_clust_keep = res$clust_gene)
  }
}

# get p-vals
p_val = sapply(RES, function(x){
  x$p_val
})
p_val = data.frame(p_val)
colnames(p_val) = methods_names

par(mfrow = c(3,3))
for(i in 1:length(methods)){
  hist(p_val[,i], main = methods[i])
}

# get adjusted p-values:
if(local){
  # get locally adjusted p-vals
  p_adj = sapply(RES, function(x){
    x$p_adj.loc
  })
}else{
  # get globally adjusted p-vals
  p_adj = sapply(RES, function(x){
    x$p_adj.glb
  })
}
p_adj = data.frame(p_adj)
colnames(p_adj) = methods_names

dim(p_val)

# remove rows with NA's
sel_NAS = rowSums(is.na(p_val)) == 0
p_val = p_val[ sel_NAS, ]
p_adj = p_adj[ sel_NAS, ]

dim(p_val)

# make UpSetR PLOT:
round(100 * colMeans(p_adj < 0.01, na.rm = TRUE))

library(UpSetR)
top_ranking = data.frame(ifelse(p_adj <= 0.05, 1, 0))

if(saving){
  dev.off()
  pdf(file = "~/Desktop/distinct project/distinct Article/v1/images/Kang/UpsetR.pdf",
      width = 7,
      height = 7,
      onefile = FALSE)
  
  upset( top_ranking, 
         nsets = 25, 
         nintersects = 25,
         keep.order = TRUE,
         order.by = "freq")
  
  dev.off()
}

upset( top_ranking, 
       nsets = 25, 
       nintersects = 25,
       keep.order = TRUE,
       order.by = "freq")

thrd = 0.05
sel = rowMeans(p_adj[,c(1:3)] < thrd) == 1 & rowMeans(p_adj[,-c(1:3)] > thrd) == 1
sum(sel)

gene_cluster = RES[[1]]$gene_clust_keep[sel_NAS]

gene_cluster = strsplit(gene_cluster, "___")

gene_cluster_vec = sapply(gene_cluster, function(x){
  paste(x[[1]], x[[2]])
})

exp = metadata(sce)$experiment_info

library(distinct)

# identify interesting profiles:
TOP_TOP = scan(what = "character", sep = ".")
.Dendritic cells EIF3K
.Dendritic cells SRSF9
.Dendritic cells MRPL41
.Dendritic cells NDUFA4
.Dendritic cells RPL24
.Dendritic cells HNRNPA0
.Dendritic cells MARCKSL1
.Dendritic cells GTF3C6
.FCGR3A+ Monocytes RPL13



TOP_TOP = TOP_TOP[TOP_TOP != ""]
TOP_TOP

# plotting function for density plot.
source("~/Desktop/distinct project/CDF plots/MY_plot_densities.R")

assays(sce)$log2_cpms = assays(sce)$logcounts

PP = list()
dev.off()
library(ggplot2)
library(ggpubr)
for(i in which(sel)) {
  cluster = gene_cluster[[i]][1]
  gene = gene_cluster[[i]][2]
  
  #if (TRUE) {
  if( paste(cluster, gene) %in% TOP_TOP){
    if (gene == "HLA-B_ENSG00000234745") {
      gene_main = "HLA-B"
    } else{
      gene_main = gene
    }
    
    print(paste(i, cluster, gene))
    
    # Density plot:
    p_1 <- plot_densities(
      x = sce,
      gene = gene,
      cluster = cluster,
      name_assays_expression = "log2_cpms"
    ) +
      geom_hline(yintercept = 0,
                 linetype = "dashed",
                 color = "grey") +
      geom_vline(xintercept = 0,
                 linetype = "dashed",
                 color = "grey") +
      theme(
        axis.text = element_text(size = rel(2)),
        axis.title = element_text(size = rel(2)),
        plot.title = element_text(size = rel(2.5)),
        legend.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        legend.key.width = unit(2, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        strip.text.x = element_text(size = rel(2)),
        legend.margin = margin()
      ) +
      labs(title = paste(cluster, "-", gene_main))

    p_3 <- plot_densities(
      x = sce,
      gene = gene,
      cluster = cluster,
      name_assays_expression = "log2_cpms",
      group_level = TRUE
    ) +
      geom_hline(yintercept = 0,
                 linetype = "dashed",
                 color = "grey") +
      geom_vline(xintercept = 0,
                 linetype = "dashed",
                 color = "grey") +
      theme(
        axis.text = element_text(size = rel(2)),
        axis.title = element_text(size = rel(2)),
        plot.title = element_text(size = rel(2)),
        legend.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        legend.key.width = unit(2, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        strip.text.x = element_text(size = rel(2)),
        legend.margin = margin()
      ) +
      labs(title = paste(cluster, "-", gene_main))
    
    p_my_densities = MY_plot_densities(
      x = sce,
      gene = gene,
      cluster = cluster,
      name_assays_expression = "log2_cpms",
      adjust = 1,
      size = 1.5
    ) +
      geom_hline(yintercept = 0,
                 linetype = "dashed",
                 color = "grey") +
      geom_vline(xintercept = 0,
                 linetype = "dashed",
                 color = "grey") +
      theme(
        axis.text = element_text(size = rel(2)),
        axis.title = element_text(size = rel(2)),
        plot.title = element_text(size = rel(2)),
        legend.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        legend.key.width = unit(2, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        strip.text.x = element_text(size = rel(2)),
        legend.margin = margin()
      ) +
      labs(title = paste(gene_main))
    
    
    # cdf plot:
    p_2 <- plot_cdfs(
      x = sce,
      gene = gene,
      cluster = cluster,
      name_assays_expression = "log2_cpms"
    ) +
      theme(
        axis.text = element_text(size = rel(2)),
        axis.title = element_text(size = rel(2)),
        plot.title = element_text(size = rel(2)),
        legend.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        legend.key.width = unit(2, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        strip.text.x = element_text(size = rel(2)),
        legend.margin = margin()
      ) +
      labs(title = paste(cluster, "-", gene_main))
    
    p_4 <- plot_cdfs(
      x = sce,
      gene = gene,
      cluster = cluster,
      name_assays_expression = "log2_cpms",
      group_level = TRUE,
      size = 1.5
    ) +
      theme(
        axis.text = element_text(size = rel(2)),
        axis.title = element_text(size = rel(2)),
        plot.title = element_text(size = rel(2)),
        legend.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        legend.key.width = unit(2, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        strip.text.x = element_text(size = rel(2)),
        legend.margin = margin()
      ) +
      labs(title = paste(cluster, "-", gene_main))
      
    
    AA = egg::ggarrange(
      plots =
        list(
          p_1 + theme(legend.position = "none") + xlab(""),
          p_2 + theme(legend.position = "none") + xlab("") + ylab("ECDF") + labs(title = ""),
          p_my_densities + theme(legend.position = "none") + labs(title = ""),
          p_4 + theme(legend.position = "none") + ylab("ECDF") + labs(title = "")
        ),
      bottom = get_legend(p_4),
      label.args = list(gp = grid::gpar(font = 2,
                                        cex = 3)),
      ncol = 2,
      nrow = 2
    )
    
    if (saving) {
      ggsave(
        filename = gsub(" ", "", paste0(cluster, "-", gene, ".pdf")),
        plot = AA,
        device = "pdf",
        path = "~/Desktop/distinct project/distinct Article/v1/images/Kang/",
        width = 16,
        height = 16,
        units = "in",
        dpi = 300,
        limitsize = TRUE
      )
    }
    
    p_1 <- p_1 + labs(title = paste(gene_main)) +
      theme( plot.title = element_text(size = rel(2)) )
    
    PP[[ which(TOP_TOP == paste(cluster, gene)) ]] = list(p_1, p_2, p_3, p_4, p_my_densities)
  }
}


for(i in seq_along(TOP_TOP)){
  PP[[i]][[5]] = PP[[i]][[5]] + theme(aspect.ratio = 0.8, # or 1 for 8 plots
                                      panel.grid = element_blank(),
                                      panel.spacing = unit(1, "mm"),
                                      panel.border = element_rect(color = "grey"),
                                      strip.text = element_text(size = 6),
                                      axis.text = element_text(size = rel(1.5)),
                                      axis.title = element_text(size=rel(1.5)),
                                      plot.title = element_text(size=rel(1.5)),
                                      legend.title=element_text(size=rel(1.5)),
                                      legend.text=element_text(size=rel(1.5)),
                                      legend.key.width=unit(1, "cm"),
                                      legend.position="bottom",
                                      legend.direction = "horizontal",
                                      legend.box="horizontal",
                                      strip.text.x = element_text(size = rel(1.5)),
                                      legend.margin=margin())
  
  PP[[i]][[1]] = PP[[i]][[1]] + theme(aspect.ratio = 0.8, # or 1 for 8 plots
                                      panel.grid = element_blank(),
                                      panel.spacing = unit(1, "mm"),
                                      panel.border = element_rect(color = "grey"),
                                      strip.text = element_text(size = 6),
                                      axis.text = element_text(size = rel(1.5)),
                                      axis.title = element_text(size=rel(1.5)),
                                      plot.title = element_text(size=rel(1.5)),
                                      legend.title=element_text(size=rel(1.5)),
                                      legend.text=element_text(size=rel(1.5)),
                                      legend.key.width=unit(1, "cm"),
                                      legend.position="bottom",
                                      legend.direction = "horizontal",
                                      legend.box="horizontal",
                                      strip.text.x = element_text(size = rel(1.5)),
                                      legend.margin=margin())
}
# ticker line for legend.
# and increase label title size.

AA = egg::ggarrange( plots = 
                       list(
                         PP[[1]][[5]] + xlab("") + theme(legend.position = "none"),
                         PP[[2]][[5]] + xlab("") + ylab("") + theme(legend.position = "none"),
                         PP[[3]][[5]] + xlab("") + ylab("") + theme(legend.position = "none"),
                         PP[[4]][[5]] + xlab("") + theme(legend.position = "none"),
                         PP[[5]][[5]] + xlab("") + ylab("") + theme(legend.position = "none"),
                         PP[[6]][[5]] + xlab("") + ylab("") + theme(legend.position = "none"),
                         PP[[7]][[5]] + theme(legend.position = "none"),
                         PP[[8]][[5]] + ylab("") + theme(legend.position = "none"),
                         PP[[9]][[5]] + ylab("") + theme(legend.position = "none")),
                     bottom = get_legend( PP[[1]][[5]] ),
                     ncol = 3, nrow = 3)
AA

if(saving){
  ggsave(filename = "9_in_1_horizontal.pdf",
         plot = AA,
         device = "pdf",
         path = "~/Desktop/distinct project/distinct Article/v1/images/Kang/",
         width = 12,
         height = 12,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}


AA = egg::ggarrange( plots = 
                       list(
                         PP[[1]][[1]] + xlab("") + theme(legend.position = "none"),
                         PP[[2]][[1]] + xlab("") + ylab("") + theme(legend.position = "none"),
                         PP[[3]][[1]] + xlab("") + ylab("") + theme(legend.position = "none"),
                         PP[[4]][[1]] + xlab("") + theme(legend.position = "none"),
                         PP[[5]][[1]] + xlab("") + ylab("") + theme(legend.position = "none"),
                         PP[[6]][[1]] + xlab("") + ylab("") + theme(legend.position = "none"),
                         PP[[7]][[1]] + theme(legend.position = "none"),
                         PP[[8]][[1]] + ylab("") + theme(legend.position = "none"),
                         PP[[9]][[1]] + ylab("") + theme(legend.position = "none")),
                     bottom = get_legend( PP[[1]][[5]] ),
                     ncol = 3, nrow = 3)
AA

if(saving){
  ggsave(filename = "9_in_1_horizontal_INDIVIDUAL.pdf",
         plot = AA,
         device = "pdf",
         path = "~/Desktop/distinct project/distinct Article/v1/images/Kang/",
         width = 12,
         height = 12,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}
