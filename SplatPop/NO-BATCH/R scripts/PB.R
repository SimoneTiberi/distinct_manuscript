setwd("~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/splatPOP/NO-BATCH/")
rm(list = ls())

library(SingleCellExperiment)

# sim parameters:
Loc   = c(0.2, 0.5, 1, 1.5, 0, 0, 0, 0)
Scale = c(0, 0, 0, 0, 0.2, 0.5, 1, 1.5)

n_sim = length(Loc)

for(ss in 1:n_sim){
  name = paste0("results/sce-Loc",Loc[ss],"-Scale",Scale[ss],".RDS")
  sce = readRDS(name)
  
  # UPDATE PB as:
  pb <- dplyr::bind_rows(
    expand.grid(
      stringsAsFactors = FALSE,
      assay =  "counts", fun = "sum", scale = FALSE,
      method = c("edgeR", "limma-voom"),
      treat = c(FALSE)
    ),
    expand.grid(
      stringsAsFactors = FALSE, scale = FALSE,
      assay = c("logcounts", "linnorm", "basics"),
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
  
  res_PB = list()
  
  set.seed(169612)
  
  library(muscat)
  for(i in 1:nrow(pb)){
    assay = pb$assay[i]
    fun = pb$fun[i]
    scale = pb$scale[i]
    
    PB <- aggregateData(sce, assay, fun = fun, scale = scale)
    
    # define the design matrix:
    group = factor(colData(PB)$group_id)
    design = model.matrix(~group)
    rownames(design) = colnames(PB) # design rownames should have sample names.
    
    res <- pbDS( pb = PB, design = design, coef = 2,
                 filter = "none", verbose = FALSE)
    res <- dplyr::bind_rows(res$table[[1]])
    
    res_PB[[i]] = res
  }
  names(res_PB) = pb$id
  
  name_res = paste0("results/res_PB-Loc",Loc[ss],"-Scale",Scale[ss],".RData")
  save(res_PB, file = name_res)
  
  print(ss)
}
