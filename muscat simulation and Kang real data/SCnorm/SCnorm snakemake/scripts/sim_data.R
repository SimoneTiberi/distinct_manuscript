suppressMessages({
  library(dplyr)
  library(jsonlite)
  library(muscat)
  library(scater)
  library(sctransform)
  library(SingleCellExperiment) 
})

# load data & simulation parameters
sce <- readRDS(args$sce)
sim_pars <- fromJSON(args$sim_pars)
set.seed(sim_pars$seed + as.numeric(wcs$i))

if(sim_pars$p_dv > 0){
  library(purrr)
  source("scripts/sim_DV.R")
  
  sim <- simData_DV(sce,
                    paired = TRUE, lfc = as.numeric(sim_pars$lfc),
                    ng = nrow(sce), nc = sim_pars$nc,
                    force = TRUE,
                    ns = sim_pars$ns, nk = sim_pars$nk,
                    p_dd = sim_pars$p_dd, probs = sim_pars$probs, p_DV = sim_pars$p_dv)
}else{
  sim <- simData(sce, 
                 paired = TRUE, lfc = as.numeric(sim_pars$lfc),
                 ng = nrow(sce), nc = sim_pars$nc,
                 force = TRUE,
                 ns = sim_pars$ns, nk = sim_pars$nk,
                 p_dd = sim_pars$p_dd, probs = sim_pars$probs)
}

sim <- sim[rowSums(counts(sim) > 0) >= 10, ]
sim <- sim[sample(nrow(sim), min(nrow(sim), sim_pars$ng)), ]

gi <- metadata(sim)$gene_info 
gi <- dplyr::filter(gi, gene %in% rownames(sim))
metadata(sim)$gene_info <- gi

# SCnorm normalized data:
library(SCnorm)

# define clusters
clusters = levels(sim$cluster_id)

# define SCnorm normalized counts:
assays(sim)$SCnorm = assays(sim)$counts

# perform normalization within each cluster
for(cl in clusters){
  # select 1 cluster only:
  sel = sim$cluster_id == cl
  
  # select counts and sample for selected cluster
  Data = assays(sim)$counts[,sel]
  Conditions = sim$sample_id[sel]
  
  # normalize data:
  norm_subset <- SCnorm(Data = Data,
                        Conditions = Conditions,
                        useZerosToScale = TRUE,
                        FilterCellNum = 10,
                        NCores=10)
  
  
  # replace selected cells with normalized data:
  assays(sim)$SCnorm[,sel] = assays(norm_subset)$normcounts
}

saveRDS(sim, args$sim)
