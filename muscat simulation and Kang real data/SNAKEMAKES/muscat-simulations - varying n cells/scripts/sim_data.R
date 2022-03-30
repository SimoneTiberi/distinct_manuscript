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

sim <- computeLibraryFactors(sim)
sim <- logNormCounts(sim)
assays(sim)$cpm <- calculateCPM(sim)
assays(sim)$vstresiduals <- suppressWarnings(
    vst(counts(sim), show_progress = FALSE)$y)

# basics:
library(BASiCS)

sim$BatchInfo = sim$sample_id

Chain <- BASiCS_MCMC(sim, WithSpikes = FALSE, Regression = FALSE,
                     N = 10^3,
                     Thin = 10,
                     Burn = 500)

assays(sim)$basics <- BASiCS_DenoisedCounts(sim, Chain)

# add log2-vstresiduals
assays(sim)$logvstresiduals <- log2(assays(sim)$vstresiduals - min(assays(sim)$vstresiduals) + 1)

# add Linnorm normalized data:
library(Linnorm)
Linnorm_data <- Linnorm.Norm( assays(sim)$counts )
assays(sim)$linnorm = Linnorm_data

saveRDS(sim, args$sim)
