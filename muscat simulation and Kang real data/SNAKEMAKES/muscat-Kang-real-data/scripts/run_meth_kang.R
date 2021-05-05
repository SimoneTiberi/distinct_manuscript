suppressMessages({
    library(scater)
    library(muscat)
    library(dplyr)
    library(jsonlite)
    library(SingleCellExperiment)
    library(sctransform)
})

sce <- readRDS(args$sce)

sce <- sce[, sce$multiplets == "singlet"] # remove multiplets
sce <- sce[, !is.na(sce$cell)] # remove unassigned cells

# REMOVE 2 samples (potential outliers): 1039
sce <- sce[, sce$ind %in% c(101, 107, 1015, 1016, 1244, 1256, 1488)] # keep control samples only

sce$cluster_id = sce$cell
sce$sample_id = paste(sce$stim, sce$ind)
sce$group_id = sce$stim
colData(sce) = colData(sce)[, -c(1:5)]
sce <- sce[rowSums(counts(sce) > 1) > 20, ]

vstresiduals = sctransform::vst( assays(sce)$counts, show_progress = FALSE)$y
assays(sce)$vstresiduals = vstresiduals

# add experimental_info via muscat::prepSCE function
sce <- prepSCE(sce, "cluster_id", "sample_id", "group_id", TRUE) # prep. SCE for `muscat`

nk <- length(kids <- levels(sce$cluster_id))
meth_pars <- as.list(fromJSON(args$meth_pars))

# increase number of threads used by AD & mixed model methods
if (grepl("AD|MM|MAST", wcs$mid)) 
    meth_pars$n_threads <- as.numeric(args$n_threads)

# run method & write results to .rds
source(fun <- args$fun)
fun <- gsub("(.R)", "", basename(fun))
res <- get(fun)(sce, meth_pars, ds_only = FALSE)$tbl

# assure all gene-cluster combinations are presents in results table
if (!inherits(res, "error")) {
    res <- mutate_at(res, c("gene", "cluster_id"), as.character)
    res <- as.data.frame(rowData(sce)) %>% 
        mutate(gene = rownames(.)) %>% 
        replicate(n = nk, simplify = FALSE) %>% 
        bind_rows %>% 
        mutate(cluster_id = rep(kids, each = nrow(sce))) %>% 
        left_join(res, by = c("gene", "cluster_id")) %>% 
        mutate(mid = wcs$mid)
}

saveRDS(res, args$res)
