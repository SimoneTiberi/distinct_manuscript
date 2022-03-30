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
sel_genes = rowSums( assays(sce)$counts>0 ) >= 20
mean(sel_genes)
sce = sce[ sel_genes , ]
rm(sel_genes)

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

sce = runPCA(sce)

set.seed(61217)
g <- buildSNNGraph(sce, use.dimred="PCA")
clust = factor(igraph::cluster_walktrap(g)$membership)
print(table(clust))

# store counts and remove sce object:
sce_counts = as.matrix(assays(sce)$counts)
rownames(sce_counts) = seq_len(nrow(sce_counts))

logcounts = as.matrix(assays(sce)$logcounts)
rownames(logcounts) = seq_len(nrow(logcounts))

cpm = as.matrix(assays(sce)$cpm)
rownames(cpm) = seq_len(nrow(cpm))

vstresiduals = as.matrix(assays(sce)$vstresiduals)
rownames(vstresiduals) = seq_len(nrow(vstresiduals))

linnorm = as.matrix(assays(sce)$linnorm)
rownames(linnorm) = seq_len(nrow(linnorm))

basics = as.matrix(assays(sce)$basics)
rownames(basics) = seq_len(nrow(basics))

rm(sce)

###############################################################################################
# define PB methods:
###############################################################################################
pb <- dplyr::bind_rows(
  expand.grid(
    stringsAsFactors = FALSE,
    assay =  "counts", fun = "sum", scale = FALSE,
    method = c("edgeR", "limma-voom"),
    treat = c(FALSE)
  ),
  expand.grid(
    stringsAsFactors = FALSE, scale = FALSE,
    assay = c("logcounts", "vstresiduals", "linnorm", "basics"),
    fun = "mean", method = c("limma-trend")
  ),
  expand.grid(
    stringsAsFactors = FALSE, scale = FALSE,
    assay = c("linnorm", "basics"),
    fun = "mean", method = c("edgeR")
  ),
  data.frame(
    stringsAsFactors = FALSE, scale = TRUE,
    assay = "cpm", fun = "sum", method = c("edgeR", "limma-trend"),
    treat = c(FALSE)
  )
)
pb$treat[is.na(pb$treat)] <- FALSE
pb$id <- with(pb, sprintf("%s%s.%s.%s%s",
                          method, ifelse(treat, "-treat", ""),
                          fun, ifelse(scale, "scale", ""), assay))

###############################################################################################
# run:
###############################################################################################
res_distinct = RES_PB = list()

library(SingleCellExperiment)
library(muscat)
library(distinct)

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
                                          cpm = cpm,
                                          logcounts = logcounts,
                                          vstresiduals = vstresiduals,
                                          linnorm = linnorm,
                                          basics = basics),
                            colData = colData_sel,
                            metadata = list(experiment_info = experiment_info_sel,
                                            n_cells = table(colData_sel$sample_id)))
  
  print(paste("sce defined - ", g))
  
  if(FALSE){
    RES_PB[[g]] = list()
    # PB methods:
    for(i in 1:nrow(pb)){
      assay = pb$assay[i]
      fun = pb$fun[i]
      scale = pb$scale[i]
      method = pb$method[i]
      
      PB <- aggregateData(x,
                          assay = assay, fun = fun,
                          scale = scale,
                          by = c("cluster_id", "sample_id"))
      
      set.seed(61217)
      try({ res <- pbDS(PB, method = method,
                        min_cells = 1,
                        verbose = FALSE, filter = "none")
      RES = do.call(rbind, res$table[[1]])
      # hist(RES$p_val)
      
      RES_PB[[g]][[i]] = RES
      }, TRUE)
    }
    names(RES_PB[[g]]) = pb$id
    
    print(paste("PB done - ", g))
  }
  
  ###############################################################################################
  # run distinct:
  ###############################################################################################
  assays = c("cpm", "logcounts", "vstresiduals", "linnorm", "basics")
  
  design_distinct = model.matrix(~metadata(x)$experiment_info$group_id)
  rownames(design_distinct) = metadata(x)$experiment_info$sample_id
  
  res_distinct[[g]] = list()
  
  library(distinct)
  for(i in 1:length(assays)){
    assay = assays[i]
    
    set.seed(61217)
    res = distinct_test(x, 
                        name_assays_expression = assay,
                        design = design_distinct,
                        min_non_zero_cells = 20,
                        n_cores = 10)
    res_distinct[[g]][[i]] = res
    
    print(paste("distinct done - ", assay))
  }
  names(res_distinct[[g]]) = assays
}

save(GROUPS, clust, n_groups, res_distinct, RES_PB, sce_counts, file = paste0("results/res_automatic_clust.RData") )