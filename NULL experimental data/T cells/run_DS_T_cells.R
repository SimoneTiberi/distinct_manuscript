R: "/usr/local/R/R-4.0.0/bin/R"
alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.0 /usr/local/R/R-4.0.0/bin/R'
export R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.0 ; /usr/bin/time -v /usr/local/R/R-4.0.0/bin/R -e '.libPaths()'
R

rm(list =ls())
load("results/sce.RData")

individuals_selected = unique(sce$sample_id)
n_samples = length(individuals_selected)

# compute ALL possible combinations of 12 individuals split into 2 groups of 6:
COMBN = combn(n_samples, 6)

# randomly select 3 splits of groups into 6 samples each:
set.seed(123456)
sel = sample.int(ncol(COMBN), 3)
GROUPS = lapply(1:3, 
                function(i){
                  factor(ifelse(seq_len(n_samples) %in%  COMBN[,sel[i] ], "A", "B"))
                })
GROUPS
n_groups = length(GROUPS)

# remove genes with < 10 non-zero cells (across ALL clusters).
sel_genes = rowSums( assays(sce)$counts>0 ) >= 10
mean(sel_genes)
sce = sce[ sel_genes , ]

# 8 patients, info abuot individuals stored here:
cd = colData(sce)
cd
table(cd$sample_id) # 12 individuals
table(cd$sample_type) # 15 clusters (for cells)

##### CHOOSE THE CLUSTER HERE:
# clust = cd$major_cluster
cd
table(cd$sample_type)
table(cd$major_cluster)

# sample_id is correlated
table(cd$sample_id, cd$sample_type)
table(cd$sample_id, cd$major_cluster)

##### CHOOSE THE CLUSTER HERE:
set.seed(61217)
library(scran)
library(scater)

g <- buildSNNGraph(sce, use.dimred="PCA")
clust = factor(igraph::cluster_walktrap(g)$membership)
print(table(clust))

# store counts and remove sce object:
sce_counts = as.matrix(assays(sce)$counts)
rownames(sce_counts) = seq_len(nrow(sce_counts))

log_CPM = as.matrix(assays(sce)$logcounts)
rownames(log_CPM) = seq_len(nrow(log_CPM))

sce_CPM = as.matrix(assays(sce)$cpms)
rownames(sce_CPM) = seq_len(nrow(sce_CPM))

# vstresiduals = sctransform::vst( as.matrix(sce_counts) )$y
# sctransform::vst crashes: we use DESeq2::vst instead:
vstresiduals = as.matrix(DESeq2::vst( as.matrix(sce_counts) ))
rownames(vstresiduals) = seq_len(nrow(vstresiduals))

rm(sce)

par(mfrow = c(2,3))
res_distinct =  RES_PB = list()
for(g in 1:n_groups){
  groups = GROUPS[[g]]
  
  experiment_info_sel = data.frame(sample_id = factor(individuals_selected),
                                   group_id = groups ) # randomly separate samples in two groups (with prob 0.5)
  
  matching = match(cd$sample_id, experiment_info_sel$sample_id)
  head(cd$sample_id); head(experiment_info_sel$sample_id[matching])
  group_id = experiment_info_sel$group_id[matching]
  
  colData_sel = list(sample_id = factor(cd$sample_id),
                     cluster_id = clust,
                     group_id = factor(group_id))
  
  x <- SingleCellExperiment(assays = list(counts = sce_counts,
                                          cpm = sce_CPM,
                                          logcounts = log_CPM,
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
                        n_cores = 4)
  }))
  res_distinct[[g]]$cpm = res
  
  print(paste("distinct cpm done - ", g))
  
  set.seed(61217)
  print(system.time({
    res = distinct_test(x, 
                        name_assays_expression = "logcounts",
                        design = design_distinct,
                        min_non_zero_cells = 10,
                        n_cores = 4)
  }))
  res_distinct[[g]]$logcounts = res
  
  print(paste("distinct logcounts done - ", g))
  
  set.seed(61217)
  print(system.time({
    res = distinct_test(x, 
                        name_assays_expression = "vstresiduals",
                        design = design_distinct,
                        min_non_zero_cells = 10,
                        n_cores = 4)
  }))
  res_distinct[[g]]$vstresiduals = res  
  
  print(paste("distinct vstresiduals done - ", g))
}

save(GROUPS, clust, n_groups, res_distinct, RES_PB, sce_counts, file = paste0("results/res_automatic_clust.RData") )
# load(paste0("results/res_automatic_clust.RData"))
