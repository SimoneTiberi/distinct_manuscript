suppressMessages({
    library(dplyr)
    library(scater)
    library(scDD)
    library(sctransform)
    library(SingleCellExperiment)
})

apply_permscdd <- function(sce, pars, ds_only = TRUE) {
    t <- system.time({
        if (ds_only) {
            normcounts(sce) <- switch(pars$assay,
                basics = assay(sce, "basics"),
                linnorm = assay(sce, "linnorm"),
                cpm = assay(sce, "cpm"),
                logcounts = 2^logcounts(sce)-1,
                vstresiduals = exp(assay(sce, "vstresiduals")))
        } else {
            normcounts(sce) <- switch(pars$assay, 
                logcounts = normcounts(logNormCounts(computeLibraryFactors(sce), log = FALSE)),
                vstresiduals = exp(vst(counts(sce), show_progress = FALSE)$y))
        }

	# only run for the 200 cell per cluster-sample combination:
	cond = ncol(sce) == 18 * 200
	if( cond ){
          res <- tryCatch(run_permscdd(sce), error = function(e) e)
	}else{
		genes = rownames(sce)
		n_genes = length(genes)

		clust = levels(sce$cluster_id)
		n_clust = length(clust)

		N = n_genes * n_clust

		res = data.frame(gene = rep(genes, n_clust),
                  cluster_id = rep(clust, each = n_genes),
                  p_val = rep(NA, N ), 
                  p_adj.loc = rep(NA, N ), 
                  p_adj.glb = rep(NA, N ))
	}
    })[[3]]
    list(rt = t, tbl = res)
}

run_permscdd <- function(sce) {
    kids <- levels(sce$cluster_id)
    cells_by_k <- split(colnames(sce), sce$cluster_id)
    if (is(normcounts(sce), "dgCMatrix"))
        normcounts(sce) <- as.matrix(normcounts(sce))
    suppressMessages(
        res <- lapply(kids, function(k) {
            res <- results(scDD(sce[, cells_by_k[[k]]], 
		permutations = 100,
                min.nonzero = 1, condition = "group_id",
                categorize = TRUE, testZeroes = FALSE,
                param = BiocParallel::MulticoreParam(workers = 3)))
            data.frame(
                gene = rownames(sce), 
                cluster_id = k,
                p_val = res$nonzero.pvalue, 
                p_adj.loc = res$nonzero.pvalue.adj,
                row.names = NULL, 
                stringsAsFactors = FALSE)
        }) 
    )
    df <- bind_rows(res)
    df$p_adj.glb <- p.adjust(df$p_val)
    return(df)
}

