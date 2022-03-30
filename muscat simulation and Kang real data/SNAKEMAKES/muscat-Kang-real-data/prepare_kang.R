suppressMessages({
    library(scater)
    library(muscat)
    library(SingleCellExperiment)
    library(sctransform)
    library(BASiCS)
    library(Linnorm)
})

sce = readRDS("data/raw_data/sce0_kang.rds")

sce <- sce[, sce$multiplets == "singlet"] # remove multiplets
sce <- sce[, !is.na(sce$cell)] # remove unassigned cells

# REMOVE 2 samples (potential outliers): 1039
sce <- sce[, sce$ind %in% c(101, 107, 1015, 1016, 1244, 1256, 1488)] # keep control samples only

sce$cluster_id = sce$cell
sce$sample_id = paste(sce$stim, sce$ind)
sce$group_id = sce$stim
colData(sce) = colData(sce)[, -c(1:5)]
sce <- sce[rowSums(counts(sce) > 0) >= 20, ]

sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)
assays(sce)$cpm <- calculateCPM(sce)

vstresiduals = sctransform::vst( assays(sce)$counts, show_progress = FALSE)$y
assays(sce)$vstresiduals = vstresiduals

# basics:
sim$BatchInfo = sim$sample_id

Chain <- BASiCS_MCMC(sim, WithSpikes = FALSE, Regression = FALSE,
                     N = 10^3,
                     Thin = 10,
                     Burn = 500)

assays(sim)$basics <- BASiCS_DenoisedCounts(sim, Chain)

# add Linnorm normalized data:
Linnorm_data <- Linnorm.Norm( assays(sim)$counts )
assays(sim)$linnorm = Linnorm_data

# add experimental_info via muscat::prepSCE function
sce <- prepSCE(sce, "cluster_id", "sample_id", "group_id", TRUE) # prep. SCE for `muscat`

saveRDS(sce, file = "data/raw_data/sce0_kang.rds")
