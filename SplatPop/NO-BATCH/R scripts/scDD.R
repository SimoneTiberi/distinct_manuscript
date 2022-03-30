setwd("~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/splatPOP/NO-BATCH/")
rm(list = ls())

library(SingleCellExperiment)

# sim parameters:
Loc   = c(0.2, 0.5, 1, 1.5, 0, 0, 0, 0)
Scale = c(0, 0, 0, 0, 0.2, 0.5, 1, 1.5)

n_sim = length(Loc)

for(i in 1:n_sim){
  name = paste0("results/sce-Loc",Loc[i],"-Scale",Scale[i],".RDS")
  sce = readRDS(name)
  
  assays = c("cpm", "linnorm", "basics")
  length(assays)
  
  library(doParallel)
  library(parallel)
  library(doRNG)
  
  suppressWarnings({
    cl <- makeCluster(6, setup_strategy = "sequential")
  })
  registerDoParallel(cl, 6)   
  
  set.seed(169612)
  
  library(scDD)
  res_scDD = foreach(aa = 1:length(assays),
                     .packages = c("scDD", "dplyr", "SingleCellExperiment"),
                     .errorhandling = "stop") %dorng% {
                       sel = which(names(assays(sce)) == assays[aa])
                       normcounts(sce) <- as.matrix(assays(sce)[[sel]])
                       
                       kids <- levels(sce$cluster_id)
                       cells_by_k <- split(colnames(sce), sce$cluster_id)
                       
                       res <- lapply(kids, function(k) {
                         res <- results(scDD(sce[, cells_by_k[[k]]],
                                             min.nonzero = 1, condition = "group_id",
                                             categorize = FALSE, testZeroes = FALSE,
                                             param = BiocParallel::MulticoreParam(workers = 1)))
                         data.frame(
                           gene = rownames(sce),
                           cluster_id = k,
                           p_val = res$nonzero.pvalue,
                           p_adj.loc = res$nonzero.pvalue.adj,
                           row.names = NULL,
                           stringsAsFactors = FALSE)
                       })
                       library(dplyr)
                       
                       df <- bind_rows(res)
                       df$p_adj.glb <- p.adjust(df$p_val)
                       df
                     }
  
  names(res_scDD) = assays
  
  name_res = paste0("results/res_scDD-Loc",Loc[i],"-Scale",Scale[i],".RData")
  save(res_scDD, file = name_res)
  
  print(i)
}

