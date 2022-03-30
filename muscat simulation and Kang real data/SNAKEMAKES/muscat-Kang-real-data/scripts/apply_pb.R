suppressMessages({
    library(limma)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
})

apply_pb <- function(sce, pars, ds_only = TRUE) {
    t1 <- system.time({
        a <- pars$assay
        pb <- aggregateData(sce, a, fun = pars$fun, scale = pars$scale)
    })[[3]]
    t2 <- system.time({
	# define the design matrix:
	group = metadata(sce)$experiment_info$group_id
	# get sample names:
	patient_id = factor(substring(metadata(sce)$experiment_info$sample_id, 6, 100))
	design = model.matrix(~group + patient_id)

      rownames(design) = colnames(pb) # design rownames should have sample names.
      res <- tryCatch(
        do.call(pbDS, c(
          # test the second column of the design matrix (group):
          list(pb = pb, design = design, coef = 2,
               filter = "none", verbose = FALSE),
          pars[names(pars) %in% names(formals(pbDS))])),
        error = function(e) e)
      if (!inherits(res, "error"))
        res <- dplyr::bind_rows(res$table[[1]])
    })[[3]]
    list(rt = c(t1, t2), tbl = res)
}
