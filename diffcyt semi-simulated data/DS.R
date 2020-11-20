#################################################
# scDiff comparisons: 'less distinct' simulations
# Lukas Weber, Mar 2020
#################################################
rm(list = ls())

setwd("~/Desktop/distinct project/SIMULATED data/FULL diffcyt comparison/DS/")

# Script to run diffcyt methods on 'less distinct' simulations using datasets
# from HDCytoData package (Weber_BCR_XL_sim)

library(diffcyt)
library(SummarizedExperiment)
library(HDCytoData)



###############
# Load datasets
###############

# Load data from HDCytoData package

data_FC <- list(
  main = Weber_BCR_XL_sim_main_SE(),
  less_50pc = Weber_BCR_XL_sim_less_distinct_less_50pc_SE(), 
  less_75pc = Weber_BCR_XL_sim_less_distinct_less_75pc_SE()
)



######################
# Run diffcyt-DS-limma
######################

# names of replicates
rep_names <- names(data_FC)

# lists to store objects
out_diffcyt_DS_limma_less_distinct  <- 
  out_DS_limma <- out_DS_LMM <-
  out_DS_byCell_limma <- out_DS_byCell_LMM <-
  out_objects_diffcyt_DS_limma_less_distinct <- 
  runtime_diffcyt_DS_limma_less_distinct <- vector("list", length(rep_names))
names(out_diffcyt_DS_limma_less_distinct) <- 
  names(out_DS_limma) <-  names(out_DS_LMM) <-
  names(out_objects_diffcyt_DS_limma_less_distinct) <- 
  names(runtime_diffcyt_DS_limma_less_distinct) <- rep_names


# run for each replicate

res_distinct = res_distinct_bycell = list()

for (s in 1:length(rep_names)) {
  
  
  ##################
  # diffcyt pipeline
  ##################
  
  # --------------------
  # pre-processing steps
  # --------------------
  
  runtime_preprocessing <- system.time({
    
    # get dataset
    d_se <- data_FC[[s]]
    
    colnames(d_se)[colData(d_se)$marker_class == "type"]
    colnames(d_se)[colData(d_se)$marker_class == "state"]
    
    # transform data
    d_se <- transformData(d_se, cofactor = 5)
    
    # clustering
    # (runtime: ~5 sec with xdim = 10, ydim = 10)
    seed <- 1234
    d_se <- generateClusters(d_se, xdim = 10, ydim = 10, seed_clustering = seed)
    
    length(table(rowData(d_se)$cluster_id))  # number of clusters
    nrow(rowData(d_se))                      # number of cells
    sum(table(rowData(d_se)$cluster_id))
    min(table(rowData(d_se)$cluster_id))     # size of smallest cluster
    max(table(rowData(d_se)$cluster_id))     # size of largest cluster
    
    # calculate cluster cell counts
    d_counts <- calcCounts(d_se)
    
    dim(d_counts)
    rowData(d_counts)
    length(assays(d_counts))
    
    # calculate cluster medians
    # (note: requires updating column names due to change in formatting of
    # column names in HDCytoData package vs. diffcyt paper)
    colnames(d_se) <- gsub("\\(.*$", "", colnames(d_se))
    
    d_medians <- calcMedians(d_se)
    
    dim(d_medians)
    rowData(d_medians)
    length(assays(d_medians))
    names(assays(d_medians))
    
    # calculate medians by cluster and marker
    d_medians_by_cluster_marker <- calcMediansByClusterMarker(d_se)
    
    dim(d_medians_by_cluster_marker)
    length(assays(d_medians_by_cluster_marker))
    
    # calculate medians by sample and marker
    d_medians_by_sample_marker <- calcMediansBySampleMarker(d_se)
    
    dim(d_medians_by_sample_marker)
    length(assays(d_medians_by_sample_marker))
    
  })
  
  
  # ---------------------------------
  # store data objects (for plotting)
  # ---------------------------------
  
  out_objects_diffcyt_DS_limma_less_distinct[[s]] <- list(
    d_se = d_se, 
    d_counts = d_counts, 
    d_medians = d_medians, 
    d_medians_by_cluster_marker = d_medians_by_cluster_marker, 
    d_medians_by_sample_marker = d_medians_by_sample_marker
  )
  
  # --------------------------------------------
  # test for differential states within clusters
  # --------------------------------------------
  
  # contrast (to compare 'spike' vs. 'base')
  # note: include fixed effects for 'patient_id'
  contrast_vec <- c(0, 1, 0, 0, 0, 0, 0, 0, 0)
  
  runtime_tests <- system.time({
    
    # set up design matrix
    # note: include fixed effects for 'patient_id'
    design <- createDesignMatrix(metadata(d_se)$experiment_info, 
                                 c("group_id", "patient_id"))
    design
    
    # set up contrast matrix
    contrast <- createContrast(contrast_vec)
    contrast
    
    # run tests
    res <- testDS_limma(d_counts, d_medians, design, contrast, plot = FALSE)
    
  })
  
  # run DS_LMM method:
  contrast_vec <- c(0, 1)
  
  runtime_tests_LMM <- system.time({
    
    # set up model formula
    # note: include random effects for 'patient_id'
    formula <- createFormula(metadata(d_se)$experiment_info, 
                             cols_fixed = "group_id", cols_random = "patient_id")
    formula
    
    # set up contrast matrix
    contrast <- createContrast(contrast_vec)
    contrast
    
    # run tests
    res_LMM <- testDS_LMM(d_counts, d_medians, formula, contrast)
  })
  
  # ---------------------------------------------
  # store results at cluster level (for plotting)
  # ---------------------------------------------
  
  res_clusters <- as.data.frame(rowData(res))
  out_DS_limma[[s]] <- res_clusters
  
  res_clusters_LMM <- as.data.frame(rowData(res_LMM))
  out_DS_LMM[[s]] <- res_clusters_LMM
  
  # select DS markers only (and first 2 cols, which are removed in the test):
  sel_cols = colData(d_se)$marker_class == "state"
  
  # invert the row/column structure:
  d_SE_sub <- SummarizedExperiment(
    assays = list(exprs = t(assays(d_se)$exprs[, sel_cols]) ), 
    colData = rowData(d_se), 
    rowData = colnames(d_se)[sel_cols], 
    metadata = list(experiment_info =  metadata(d_se)$experiment_info )
  )
  
  design_distinct = model.matrix(~metadata(d_SE_sub)$experiment_info$group_id)
  rownames(design_distinct) = metadata(d_SE_sub)$experiment_info$sample_id
  
  library(distinct)
  set.seed(61217)
  print(system.time({
    res_DS_cells = distinct_test(d_SE_sub, 
                                 name_assays_expression = "exprs",
                                 name_sample = "sample_id",
                                 design = design_distinct,
                                 min_non_zero_cells = 0,
                                 n_cores = 4)
  }))
  res_distinct[[s]] = res_DS_cells
  
  # marker_id is the same as gene id!
  all(levels(res_DS_cells$gene) %in% levels(rowData(res)$marker_id))
  
  ##############################
  # Return results at cell level
  ##############################
  
  # Note: diffcyt methods return results for each cluster-marker combination. To enable
  # performance comparisons between methods at the cell level, we assign the same p-values
  # to all cells within a given cluster-marker combination.
  
  # Note: return cell-level results for marker pS6 only, since the comparative evaluations
  # (ROC curves etc) are based on pS6 in B cells only.
  
  # identify B cells (these contain the true differential signal)
  is_B_cell <- rowData(d_se)$B_cell
  
  # match cluster-level p-values for marker pS6 to individual cells
  
  n_state_markers <- length(colnames(d_se)[colData(d_se)$marker_class == "state"])
  
  stopifnot(nrow(rowData(res)) == nlevels(rowData(d_se)$cluster_id) * n_state_markers)
  stopifnot(all(levels(rowData(res)$cluster_id) == levels(rowData(d_se)$cluster_id)))
  stopifnot(all(levels(rowData(res)$cluster_id) %in% rowData(res)$cluster_id))
  
  stopifnot(nrow(rowData(res_LMM)) == nlevels(rowData(d_se)$cluster_id) * n_state_markers)
  stopifnot(all(levels(rowData(res_LMM)$cluster_id) == levels(rowData(d_se)$cluster_id)))
  stopifnot(all(levels(rowData(res_LMM)$cluster_id) %in% rowData(res_LMM)$cluster_id))
  
  # select results for pS6
  res_pS6 <- res[rowData(res)$marker_id == "pS6", ]
  res_pS6_LMM <- res_LMM[rowData(res)$marker_id == "pS6", ]
  
  res_DS_cells_pS6 = res_DS_cells[res_DS_cells$gene == "pS6", ]
  
  # match cells to clusters
  ix_match     <- match(rowData(d_se)$cluster_id, rowData(res_pS6)$cluster_id)
  ix_match_LMM <- match(rowData(d_se)$cluster_id, rowData(res_pS6_LMM)$cluster_id)
  
  ix_match_cells   <- match(rowData(d_se)$cluster_id, res_DS_cells_pS6$cluster_id)
  
  p_vals_clusters <- rowData(res_pS6)$p_val
  p_adj_clusters <- rowData(res_pS6)$p_adj
  p_vals_cells <- p_vals_clusters[ix_match]
  p_adj_cells <- p_adj_clusters[ix_match]
  
  p_vals_clusters_LMM <- rowData(res_pS6_LMM)$p_val
  p_adj_clusters_LMM <- rowData(res_pS6_LMM)$p_adj
  p_vals_cells_LMM <- p_vals_clusters_LMM[ix_match_LMM]
  p_adj_cells_LMM <- p_adj_clusters_LMM[ix_match_LMM]
  
  p_vals_DS_cells <- res_DS_cells_pS6$p_val[ix_match_cells]
  p_adj_DS_cells  <- res_DS_cells_pS6$p_adj.glb[ix_match_cells]
  
  # set up data frame with results (for marker pS6) and true B-cell status at cell level
  
  res_p_vals <- p_vals_cells
  res_p_adj <- p_adj_cells
  
  res_p_vals_LMM <- p_vals_cells_LMM
  res_p_adj_LMM <- p_adj_cells_LMM
  
  # replace NAs (due to filtering) to ensure same cells are returned for all methods
  res_p_vals[is.na(res_p_vals)] <- 1
  res_p_adj[is.na(res_p_adj)] <- 1
  
  res_p_adj_LMM[is.na(res_p_adj_LMM)] <- 1
  res_p_adj_LMM[is.na(res_p_adj_LMM)] <- 1
  
  stopifnot(length(res_p_vals) == length(res_p_adj))
  stopifnot(length(res_p_vals) == length(is_B_cell))
  
  stopifnot(length(res_p_vals_LMM) == length(res_p_adj_LMM))
  stopifnot(length(res_p_vals_LMM) == length(is_B_cell))
  
  res <- data.frame(p_val = res_p_vals, 
                    p_adj = res_p_adj, 
                    B_cell = is_B_cell)
  # store results
  out_DS_byCell_limma[[s]] <- res
  
  res_LMM <- data.frame(p_val = res_p_vals_LMM, 
                        p_adj = res_p_adj_LMM, 
                        B_cell = is_B_cell)
  # store results
  out_DS_byCell_LMM[[s]] <- res_LMM
  
  res_DS_cells <- data.frame(p_val = p_vals_DS_cells, 
                             p_adj = p_adj_DS_cells, 
                             B_cell = is_B_cell)
  # store results
  res_distinct_bycell[[s]] <- res_DS_cells
}

#same size?
dim(out_DS_limma[[1]]); 
dim(out_DS_LMM[[1]]); 
dim(res_distinct[[1]]);

# all B-cells are differential (in the spikein group)
rd = rowData(data_FC[[1]])
table(rd$B_cell, rd$spikein, rd$group_id)

# save.image("DS/results/DS_by_cluster_CellLevel.RData")

##################################################################################################
# Plots:
##################################################################################################
rm(list =ls())
setwd("~/Desktop/distinct project/SIMULATED data/FULL diffcyt comparison/DS")
load("results/DS_by_cluster_CellLevel.RData")
dev.off()

library(ggplot2)
library(RColorBrewer);
colours = c(
  brewer.pal(4, "BuPu")[4],    
  brewer.pal(4, "Greens")[4],
  brewer.pal(4, "Reds")[4],    
  "white")


gg_roc = gg_fdr = DF_list = list()
library(iCOBRA)
# make ROC and FDR plots:
for (s in 1:length(rep_names)) {
  res_limma = out_DS_byCell_limma[[s]][ order(out_DS_byCell_limma[[s]]$B_cell),]
  res_LMM = out_DS_byCell_LMM[[s]][ order(out_DS_byCell_LMM[[s]]$B_cell),]
  
  res_DS_cells = res_distinct_bycell[[s]]
  res_DS_cells = res_DS_cells[order(res_DS_cells$B_cell),]
  
  stopifnot(ncol(res_limma) == ncol(res_LMM))
  
  p_icobra = data.frame( distinct = res_DS_cells$p_val,
                         diffcyt_limma = res_limma$p_val,
                         diffcyt_LMM = res_LMM$p_val)
  
  padj_icobra = data.frame( distinct = res_DS_cells$p_adj,
                            diffcyt_limma = res_limma$p_adj,
                            diffcyt_LMM = res_LMM$p_adj)
  
  truth_cobra = data.frame(status = ifelse(res_limma$B_cell,1,0) )
  
  cobra = COBRAData(pval = p_icobra,
                    padj = padj_icobra,
                    truth = truth_cobra,
                    object_to_extend = NULL)
  
  cobraperf <- calculate_performance(cobra, binary_truth = "status")
  
  cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colours,
                                     facetted = TRUE, incloverall = FALSE,
                                     conditionalfill = FALSE)
  
  gg_roc[[s]] = plot_roc(cobraplot, linewidth = 1.5) +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          legend.text=element_text(size=rel(1.5)),
          axis.text=element_text(size=rel(2)),
          axis.title=element_text(size=rel(1.5)),
          legend.key.width=unit(2, "cm"),
          legend.title=element_text(size=rel(1.5)),
          aspect.ratio = 1,
          axis.text.x = element_text(angle = 45, hjust=0.5) )# +
  
  gg_fdr[[s]] = plot_fdrtprcurve(cobraplot, plottype = c("points", "curve"),
                                 pointsize = 4, linewidth = 1.5) +
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)) +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          legend.text=element_text(size=rel(1.5)),
          axis.text=element_text(size=rel(2)),
          axis.title=element_text(size=rel(1.5)),
          legend.key.width=unit(2, "cm"),
          legend.title=element_text(size=rel(1.5)),
          aspect.ratio = 1,
          axis.text.x = element_text(angle = 45, hjust=0.5) ) +
    geom_point(data = cobraplot@fdrtpr, 
               size = 4, alpha = 0.1, 
               colour = rep(colours[-4], 3))
  
  # make FDR table:
  DF = data.frame(method = cobraplot@fdrtpr$method, 
                  thr = as.numeric(substring(as.character(cobraplot@fdrtpr$thr), first = 4)),
                  FDR = cobraplot@fdrtpr$FDR,
                  TPR = cobraplot@fdrtpr$TPR)
  DF[order(DF$method, decreasing = TRUE),]
  
  DF_list[[s]] = DF
}
DF = do.call(cbind, DF_list)

library(xtable)
xtable( DF[,-c(5,6,9,10)], digits = 2, align = "|ll|r|rr|rr|rr|")

gg_tmp = gg_fdr[[1]] + theme(legend.text=element_text(size=rel(1.5)),
                             legend.key.width=unit(1.5, "cm"),
                             legend.title=element_text(size=rel(1.5)),
                             legend.box="horizontal", legend.position = "bottom")

# FDR AND ROC VERTICAL:
library(ggpubr)
ggarrange(gg_roc[[1]] + xlab(""),
          gg_tmp + xlab("") + ylab(""),
          gg_roc[[2]] + xlab(""),
          gg_fdr[[2]] + xlab("") + ylab(""),
          gg_roc[[3]],
          gg_fdr[[3]] + ylab(""),
          legend.grob = get_legend(gg_tmp),
          labels = c( "",
                      "main",
                      "",
                      "less 50", 
                      "", 
                      "less 75"),
          label.x = -0.1,
          font.label = list(size = 20),
          ncol = 2, nrow = 3,
          legend = "bottom", common.legend = TRUE)
# DS_SingleCellLevel.pdf
# 16 * 16

ggsave(filename = "diffcyt_SingleCell_ROC_FDR_vertical.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "~/Desktop/distinct project/distinct Article/v1/images/diffcyt",
       width = 12,
       height = 16,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

# FDR HORIZONTAL:
AA = egg::ggarrange( plots = 
                       list(gg_fdr[[1]] + xlab("") + labs(title = "main") + 
                              theme(legend.position = "none", 
                                    plot.title = element_text(hjust = 0.5, face = "bold", size=20)),
                            gg_fdr[[2]] + ylab("") + labs(title = "less 50") + 
                              theme(legend.position = "none", 
                                    plot.title = element_text(hjust = 0.5, face = "bold", size=20)),
                            gg_fdr[[3]] + xlab("") + ylab("") + labs(title = "less 75") + 
                              theme(legend.position = "none", 
                                    plot.title = element_text(hjust = 0.5, face = "bold", size=20))),
                     bottom = ggpubr::get_legend(gg_tmp),
                     ncol = 3, nrow = 1)


ggsave(filename = "diffcyt_SingleCell_FDR.pdf",
       plot = AA,
       device = "pdf",
       path = "~/Desktop/distinct project/distinct Article/v1/images/diffcyt",
       width = 12,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

# ROC HORIZONTAL:
AA = egg::ggarrange( plots = 
                       list(gg_roc[[1]] + xlab("") + labs(title = "main") + 
                              theme(legend.position = "none", 
                                    plot.title = element_text(hjust = 0.5, face = "bold", size=20)),
                            gg_roc[[2]] + ylab("") + labs(title = "less 50") + 
                              theme(legend.position = "none", 
                                    plot.title = element_text(hjust = 0.5, face = "bold", size=20)),
                            gg_roc[[3]] + xlab("") + ylab("") + labs(title = "less 75") + 
                              theme(legend.position = "none", 
                                    plot.title = element_text(hjust = 0.5, face = "bold", size=20))),
                     bottom = ggpubr::get_legend(gg_tmp),
                     ncol = 3, nrow = 1)


ggsave(filename = "diffcyt_SingleCell_ROC.pdf",
       plot = AA,
       device = "pdf",
       path = "~/Desktop/distinct project/distinct Article/v1/images/diffcyt",
       width = 12,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
