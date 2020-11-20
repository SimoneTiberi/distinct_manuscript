suppressMessages({
    library(muscat)
    library(SingleCellExperiment)
})

apply_mm <- function(sce, pars, ds_only = TRUE) {
    pars <- pars[names(pars) != "id"]
    if (is.null(pars$n_threads))
        pars$n_threads <- 1
    t <- system.time({
        pars[sapply(pars, `==`, "")] <- NULL
	
	# add patient id as covariate:
	sce$patient_id = as.numeric(factor(substring(sce$sample_id, 1, 7)))
	sce$batch_id = as.numeric(colData(sce)$batch_id)
	pars$covs = c("batch_id", "patient_id")

        res <- tryCatch(
            do.call(mmDS, c(list(sce, verbose = FALSE), pars)),
            error = function(e) e)
        if (!inherits(res, "error"))
            res <- dplyr::bind_rows(res)
    })[[3]]
    list(rt = t, tbl = res)
}
