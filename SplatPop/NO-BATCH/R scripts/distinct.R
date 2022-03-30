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
  
  library(distinct)

  group = metadata(sce)$experiment_info$group_id
  design = model.matrix(~group)
  # rownames of the design must indicate sample ids:
  rownames(design) = metadata(sce)$experiment_info$sample_id
  
  assays = c( "logcounts", "cpm", "vstresiduals", "linnorm", "basics")
  length(assays)
  
  library(doParallel)
  library(parallel)
  library(doRNG)
  
  suppressWarnings({
    cl <- makeCluster(6, setup_strategy = "sequential")
  })
  registerDoParallel(cl, 6)   
  
  set.seed(169612)
  res_distinct = foreach(aa = 1:length(assays),
                         .packages = c("distinct"),
                         .errorhandling = "stop") %dorng% {
                           res <- distinct_test(sce,
                                                name_assays_expression = assays[[aa]],
                                                n_cores = 1,
                                                design = design,
                                                column_to_test = 2,
                                                min_non_zero_cells = 1)
                           res
                         }
  names(res_distinct) = assays
  
  name_res = paste0("results/res_distinct-Loc",Loc[i],"-Scale",Scale[i],".RData")
  save(res_distinct, file = name_res)
  
  print(i)
}