suppressMessages({
    library(dplyr)
    library(jsonlite)
    library(muscat)
    library(muscat.batch)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
    library(magrittr)
})

# load data & simulation parameters
sce <- readRDS(args$sce)
sim_pars <- fromJSON(args$sim_pars)
set.seed(sim_pars$seed + as.numeric(wcs$i))

# define the design of the batch:
batch_design = matrix(c(1,1,2, 2,1,2), ncol = 2, byrow = FALSE)

# simulate from muscat.batch:::simData (a modified version of muscat::simData)
if(sim_pars$p_dv > 0){
  library(purrr)
  source("scripts/sim_DV_batch.R")
  
  sim <- simData_DV_batch(sce,
                          paired = TRUE,
                          lfc = as.numeric(sim_pars$lfc),
                          ng = nrow(sce), nc = sim_pars$nc,
                          ns = sim_pars$ns, nk = sim_pars$nk,
                          p_dd = sim_pars$p_dd, probs = sim_pars$probs,
                          rel_be_c = rep(0.4, sim_pars$nk),
                          lfc_be = 0.5, batch_design = batch_design, nb = 2,
                          p_DV = sim_pars$p_dv)
}else{
  sim <- muscat.batch:::simData(sce,
                                paired = TRUE,
                                lfc = as.numeric(sim_pars$lfc),
                                ng = nrow(sce), nc = sim_pars$nc,
                                ns = sim_pars$ns, nk = sim_pars$nk,
                                p_dd = sim_pars$p_dd, probs = sim_pars$probs,
                                rel_be_c = rep(0.4, sim_pars$nk),
                                lfc_be = 0.5, batch_design = batch_design, nb = 2)
}

# store batch in the experimental info:
MD = metadata(sim)$experiment_info
MD$batch_id = factor(unlist(batch_design))
metadata(sim)$experiment_info = MD

# store batch id for every cell:
matching = match(sim$sample_id, MD$sample_id)
batch = MD$batch_id[matching]
  
sim$batch_id = batch

# cluster_id must be factor:
sim$cluster_id = factor(sim$cluster_id)

# filter lowly abundant genes:
sim <- sim[rowSums(counts(sim) > 0) >= 10, ]
sim <- sim[sample(nrow(sim), min(nrow(sim), sim_pars$ng)), ]

gi <- metadata(sim)$gene_info 
gi <- dplyr::filter(gi, gene %in% rownames(sim))
metadata(sim)$gene_info <- gi

sim <- computeLibraryFactors(sim)
sim <- logNormCounts(sim)
assays(sim)$cpm <- calculateCPM(sim)
assays(sim)$vstresiduals <- suppressWarnings(
    vst(counts(sim), show_progress = FALSE)$y)

# add log2-vstresiduals
assays(sim)$logvstresiduals <- log2(assays(sim)$vstresiduals - min(assays(sim)$vstresiduals) + 1)

saveRDS(sim, args$sim)
