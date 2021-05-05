args <- R.utils::commandArgs(
    trailingOnly = TRUE, 
    asValues = TRUE)

suppressMessages({
    library(muscat)
    library(SingleCellExperiment)
})

sce <- readRDS(args$input_sce) # load data
reducedDims(sce) <- NULL # remove dimensionality reductions
sce <- sce[, sce$multiplets == "singlet"] # remove multiplets
sce <- sce[, !is.na(sce$cell)] # remove unassigned cells
sce <- sce[, sce$stim == "ctrl"] # keep control samples only

# REMOVE 2 samples (potential outliers): 1488 and 1039
sce <- sce[, sce$ind %in% c(101, 107, 1015, 1016, 1244, 1256, 1488)] # keep control samples only

sce$sample_id <- factor(paste0(sce$stim, sce$ind)) # construct sample IDs
sce <- prepSCE(sce, "cell", "sample_id", "stim", TRUE) # prep. SCE for `muscat`
saveRDS(sce, args$output_sce) # write SCE to .rds
