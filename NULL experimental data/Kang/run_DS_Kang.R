R: "/usr/local/R/R-4.0.0/bin/R"
alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.0 /usr/local/R/R-4.0.0/bin/R'
export R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.0 ; /usr/bin/time -v /usr/local/R/R-4.0.0/bin/R -e '.libPaths()'
R

rm(list =ls())

library(muscData)
sce = Kang18_8vs8()

# remove multiplets and unassigned cells
sce <- sce[, sce$multiplets == "singlet"] # remove multiplets
sce <- sce[, !is.na(sce$cell)] # remove unassigned cells

# REMOVE 1 samples (potential outlier): 1039
sce <- sce[, sce$ind %in% c(101, 107, 1015, 1016, 1244, 1256, 1488)] # keep control samples only

# Select Vehicle mouse only:
sel = sce$stim == "ctrl"
mean(sel)
sce = sce[, sel ]

# remove genes with < 10 non-zero cells (across ALL clusters).
sel_genes = rowSums( assays(sce)$counts>0 ) >= 10
mean(sel_genes)
sce = sce[ sel_genes , ]

sce
colData(sce)

# 8 patients, info abuot individuals stored here:
cd = colData(sce)
cd
table(cd$ind) # 8 individuals
table(cd$stim) # 2 experimental conditions: control (ctrl) and stimulated (stim)
table(cd$ind, cd$stim)
table(cd$cluster) # 10 clusters
table(cd$cell) # 8 cell types
table(cd$cluster, cd$cell) # cluster is not cell

individuals_selected = unique(cd$ind)

# make 3 random splits of samples in 2 groups:
n_groups = 3

n_samples = length(individuals_selected); n_samples
COMBN = combn(n_samples, 3)

set.seed(123456)
sel = sample.int(ncol(COMBN), n_groups)
GROUPS = lapply(1:3, 
                function(i){
                  factor(ifelse(seq_len(n_samples) %in%  COMBN[,sel[i] ], "A", "B"))
                })
GROUPS

sce_counts = assays(sce)$counts

##### CHOOSE THE CLUSTER HERE:
clust = factor(cd$cell)

library(scran)
library(scater)

sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)
assays(sce)$cpm <- calculateCPM(sce)

# store counts and remove sce object:
sce_counts = assays(sce)$counts
rownames(sce_counts) = seq_len(nrow(sce_counts))

log_CPM = assays(sce)$logcounts
rownames(log_CPM) = seq_len(nrow(log_CPM))

sce_CPM = as.matrix(assays(sce)$cpm)
rownames(sce_CPM) = seq_len(nrow(sce_CPM))

vstresiduals = sctransform::vst( as.matrix(sce_counts) )$y

rm(sce)

res_distinct = RES_PB = list()

par(mfrow = c(3,6))
for(g in 1:n_groups){
  groups = GROUPS[[g]]
  
  experiment_info_sel = data.frame(sample_id = factor(individuals_selected),
                                   group_id = groups ) # randomly separate samples in two groups (with prob 0.5)
  
  matching = match(cd$ind, experiment_info_sel$sample_id)
  head(cd$ind); head(experiment_info_sel$sample_id[matching])
  group_id = experiment_info_sel$group_id[matching]
  
  colData_sel = list(sample_id = factor(cd$ind),
                     cluster_id = clust,
                     group_id = factor(group_id))
  
  x <- SingleCellExperiment(assays = list(counts = as.matrix(sce_counts),
                                          cpm = as.matrix(sce_CPM),
                                          logcounts = as.matrix(log_CPM),
                                          vstresiduals = vstresiduals),
                            colData = colData_sel,
                            metadata = list(experiment_info = experiment_info_sel,
                                            n_cells = table(colData_sel$sample_id)))

  
  # WITHOUT as.matrix(sce_counts)) aggregateData fails!
  
  print(paste("x created - ", g))
  
  # PB WITH COUNTS:
  # PSEUDO-BULK on raw counts
  library(muscat)
  pb <- aggregateData(x,
                      assay = "counts", fun = "sum",
                      by = c("cluster_id", "sample_id"))
  assayNames(pb)
  
  RES_PB[[g]] = list()
  
  set.seed(61217)
  try({ res <- pbDS(pb, method = c("edgeR"),
                    min_cells = 1,
                    verbose = FALSE, filter = "none")
  RES = do.call(rbind, res$table[[1]])
  hist(RES$p_val)
  
  RES_PB[[g]]$edgeR.counts = RES
  }, TRUE)
  
  set.seed(61217)
  try({ res <- pbDS(pb, method = c("limma-voom"),
                    min_cells = 1,
                    verbose = FALSE, filter = "none")
  RES = do.call(rbind, res$table[[1]])
  hist(RES$p_val)
  
  RES_PB[[g]]$limma_voom.counts = RES
  }, TRUE)
  
  
  # PB WITH CPM:
  # PSEUDO-BULK on raw counts
  rm(pb)
  pb <- aggregateData(x,
                      assay = "cpm", fun = "sum",
                      by = c("cluster_id", "sample_id"))
  assayNames(pb)
  
  set.seed(61217)
  try({ res <- pbDS(pb, method = c("edgeR"),
                    min_cells = 1,
                    verbose = FALSE, filter = "none")
  RES = do.call(rbind, res$table[[1]])
  hist(RES$p_val)
  
  RES_PB[[g]]$edgeR.cpm = RES
  }, TRUE)
  
  
  # PB WITH logcounts:
  # PSEUDO-BULK on raw counts
  rm(pb)
  pb <- aggregateData(x,
                      assay = "logcounts", fun = "sum",
                      by = c("cluster_id", "sample_id"))
  assayNames(pb)
  
  set.seed(61217)
  try({ res <- pbDS(pb, method = c("limma-trend"),
                    min_cells = 1,
                    verbose = FALSE, filter = "none")
  RES = do.call(rbind, res$table[[1]])
  hist(RES$p_val)
  
  RES_PB[[g]]$limma_trend.logcounts = RES
  }, TRUE)
  
  
  # PB WITH vstresid:
  # PSEUDO-BULK on raw counts
  rm(pb)
  pb <- aggregateData(x,
                      assay = "vstresiduals", fun = "sum",
                      by = c("cluster_id", "sample_id"))
  assayNames(pb)
  
  set.seed(61217)
  try({ res <- pbDS(pb, method = c("limma-trend"),
                    min_cells = 1,
                    verbose = FALSE, filter = "none")
  RES = do.call(rbind, res$table[[1]])
  hist(RES$p_val)
  
  RES_PB[[g]]$limma_trend.vstresiduals = RES
  }, TRUE)
  
  
  print(paste("PB done - ", g))
  
  # run distinct:
  design_distinct = model.matrix(~metadata(x)$experiment_info$group_id)
  rownames(design_distinct) = metadata(x)$experiment_info$sample_id
  
  res_distinct[[g]] = list()
  
  library(distinct)
  set.seed(61217)
  print(system.time({
    res = distinct_test(x, 
                        name_assays_expression = "cpm",
                        design = design_distinct,
                        min_non_zero_cells = 10,
                        n_cores = 10)
  }))
  res_distinct[[g]]$cpm = res
  
  print(paste("distinct cpm done - ", g))
  
  set.seed(61217)
  print(system.time({
    res = distinct_test(x, 
                        name_assays_expression = "logcounts",
                        design = design_distinct,
                        min_non_zero_cells = 10,
                        n_cores = 10)
  }))
  res_distinct[[g]]$logcounts = res
  
  print(paste("distinct logcounts done - ", g))
  
  set.seed(61217)
  print(system.time({
    res = distinct_test(x, 
                        name_assays_expression = "vstresiduals",
                        design = design_distinct,
                        min_non_zero_cells = 0,
                        n_cores = 10)
  }))
  res_distinct[[g]]$vstresiduals = res  
  
  print(paste("distinct vstresiduals done - ", g))
}

save(GROUPS, clust, n_groups, res_distinct, RES_PB, sce_counts, file = paste0("results/res_cell_clust_7Samples.RData") )
# load(paste0("results/res_cell_clust_7Samples.RData"))

