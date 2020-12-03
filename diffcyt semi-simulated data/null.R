R: "/usr/local/R/R-4.0.0/bin/R"
alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.0 /usr/local/R/R-4.0.0/bin/R'
export R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.0 ; /usr/bin/time -v /usr/local/R/R-4.0.0/bin/R -e '.libPaths()'
R

rm(list = ls())

setwd("~/Desktop/distinct project/EXPERIMENTAL data/Bodenmiller (diffcyt)")

# Script to run scDiff and diffcyt methods on null simulations from diffcyt
# paper (using Weber_BCR_XL_sim_null datasets, which are based on the original
# Bodenmiller_BCR_XL dataset)


library(diffcyt)
library(SummarizedExperiment)
library(HDCytoData)
library(distinct)


###############
# Load datasets
###############

# Load null simulation datasets from HDCytoData package

data_null_sims <- list(
  rep1 = Weber_BCR_XL_sim_null_rep1_SE(), 
  rep2 = Weber_BCR_XL_sim_null_rep2_SE(), 
  rep3 = Weber_BCR_XL_sim_null_rep3_SE()
)


######################
# Run diffcyt-DS-limma
######################

# names of replicates
rep_names <- names(data_null_sims)

# lists to store objects
out_diffcyt_DS_limma_null  <-
  out_DS_LMM <-
  out_DS_limma <- 
  out_objects_diffcyt_DS_limma_null <- 
  runtime_diffcyt_DS_limma_null <- vector("list", length(rep_names))
names(out_diffcyt_DS_limma_null) <- 
  names(out_DS_limma) <- 
  names(out_objects_diffcyt_DS_limma_null) <- 
  names(runtime_diffcyt_DS_limma_null) <- rep_names


# run for each replicate

res_scDiff_permSamples = res_distinct = list()

for (s in 1:length(rep_names)) {
  ##################
  # diffcyt pipeline
  ##################
  
  # --------------------
  # pre-processing steps
  # --------------------
  
  runtime_preprocessing <- system.time({
    
    # get dataset
    d_se <- data_null_sims[[s]]
    
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
    d_counts <- calcCounts(d_se) # summarizes counts at the cluster level: 1 number per cluster.
    
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
  
  out_objects_diffcyt_DS_limma_null[[s]] <- list(
    d_se = d_se, 
    d_counts = d_counts, 
    d_medians = d_medians, 
    d_medians_by_cluster_marker = d_medians_by_cluster_marker, 
    d_medians_by_sample_marker = d_medians_by_sample_marker
  )
  
  # --------------------------------------------
  # test for differential states within clusters
  # --------------------------------------------
  
  # contrast (to compare 'null2' vs. 'null1')
  # note: include fixed effects for 'patient_id'
  contrast_vec <- c(0, 1, 0, 0, 0, 0, 0, 0, 0)
  
  runtime_tests <- system.time({
    
    # set up design matrix
    # note: include fixed effects for 'patient_id'
    design <- createDesignMatrix(metadata(d_se)$experiment_info, cols_design = 1:2)
    design
    
    # set up contrast matrix
    contrast <- createContrast(contrast_vec)
    contrast
    
    # run tests
    res <- testDS_limma(d_counts, d_medians, design, contrast, path = ".")
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
  
  # ---------------------------------------------
  # store results at cluster level (for plotting)
  # ---------------------------------------------
  
  res_clusters <- as.data.frame(rowData(res))
  out_DS_limma[[s]] <- res_clusters
  
  res_clusters_LMM <- as.data.frame(rowData(res_LMM))
  out_DS_LMM[[s]] <- res_clusters_LMM
}

# save.image("results/null by cluster.RData")

##################################################################################################
# Plots:
##################################################################################################

rm(list =ls())

setwd("~/Desktop/distinct project/EXPERIMENTAL data/Bodenmiller (diffcyt)")
load("results/null by cluster.RData")

# Make nice null plot:

suppressMessages({
  library(data.table)
  library(dplyr)
  library(iCOBRA)
  library(ggplot2)
  library(purrr)
  library(reshape2)
})

for(i in 1:3){
  res_distinct[[i]]$replicate = i
  out_DS_limma[[i]]$replicate = i
  out_DS_LMM[[i]]$replicate = i
}

res_distinct = do.call(rbind, res_distinct)
res_DS_limma = do.call(rbind, out_DS_limma)
res_DS_LMM = do.call(rbind, out_DS_LMM)

head(res_distinct)

RES = data.frame(p_val = c(res_distinct$p_val, res_DS_limma$p_val, res_DS_LMM$p_val),
                 method = rep(c("distinct", "diffcyt_limma", "diffcyt_LMM"), each = nrow(res_distinct)),
                 replicate = c(res_distinct$replicate, res_DS_limma$replicate, res_DS_LMM$replicate) 
                 , stringsAsFactors = TRUE)
head(RES)

class(RES$p_val)
class(RES$method)

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


RES$replicate = as.factor(RES$replicate)

# scale colours of lines and inside:
library(RColorBrewer);
colours = c(
  brewer.pal(4, "BuPu")[4],    
  brewer.pal(4, "Greens")[4],
  brewer.pal(4, "Reds")[4])



p <- ggplot(data = RES, aes(x = p_val, y = ..ndensity.., 
                            col = method, fill = method, lty = replicate)) +
  facet_wrap(~ method, ncol = 4) + 
  geom_density(adjust = 0.5, size = 0.3, alpha = 0.1) +
  scale_alpha_manual(values = c("TRUE" = 0.1, "FALSE" = 0.4)) +
  guides(col = FALSE,
         lty = guide_legend(ncol = 1, order = 2),
         fill = guide_legend(ncol = 1, order = 1,
         override.aes = list(alpha = 1, col = NA))) +
  scale_x_continuous("p-value", breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
  scale_y_continuous("density", breaks = c(0, 1), expand = c(0, 0.1)) +
  .prettify("bw") + theme(aspect.ratio = 0.5,
                          panel.grid = element_blank(),
                          panel.spacing = unit(1, "mm"),
                          panel.border = element_rect(color = "grey"),
                          strip.text = element_text(size = 5),
                          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(1.5)),
                          axis.text.y=element_text(size=rel(1.5)),
                          axis.title.y = element_text(size=rel(1.5)),
                          axis.title.x = element_text(size=rel(1.5)),
                          axis.title.x.top = element_text(size=rel(1.5)),
                          legend.title=element_text(size=rel(1.5)),
                          legend.text=element_text(size=rel(1.5)),
                          legend.key.width=unit(0.5, "cm"),
                          #legend.box.just = "top",
                          legend.position="bottom",
                          legend.direction = "horizontal",
                          legend.box="horizontal",
                          strip.text.x = element_text(size = rel(1.5)),
                          legend.margin=margin()) +
  scale_colour_manual(values = colours) +
  scale_fill_manual(values = colours)
p

ggsave(filename = paste0("diffcyt-null-by-cluster.pdf"),
       plot = last_plot(),
       device = "pdf",
       path = "~/Desktop/distinct project/distinct Article/v1/images/diffcyt",
       width = 8,
       height = 3,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
# Figure 1c
