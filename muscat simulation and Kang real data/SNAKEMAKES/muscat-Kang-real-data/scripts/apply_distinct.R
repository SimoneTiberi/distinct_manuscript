suppressMessages({
  library(muscat)
  library(SingleCellExperiment)
  library(distinct)
})

apply_distinct <- function(sce, pars, ds_only = TRUE) {
  t <- system.time({
    a <- pars$assay
    cores = as.numeric(pars$cores)

    # create the design of the study:
    samples = metadata(sce)$experiment_info$sample_id

    patient_id = factor(substring(metadata(sce)$experiment_info$sample_id, 6, 100))

    group = metadata(sce)$experiment_info$group_id
    design = model.matrix(~group + patient_id)
    # rownames of the design must indicate sample ids:
    rownames(design) = samples

    # by default, using counts, then try all the other ones too:
    res <- tryCatch(distinct_test(sce, 
                                  name_assays_expression = a,
				  n_cores = cores,
				  design = design,
				  column_to_test = 2,
                                  min_non_zero_cells = 1), error = function(e) e)
    
    res$gene = as.character(res$gene)
    res$cluster_id = as.character(res$cluster_id)
  })[[3]]
  list(rt = t, tbl = res)
}
