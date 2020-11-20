suppressMessages({
  library(muscat)
  library(SingleCellExperiment)
  library(distinct)
})

apply_distinct <- function(sce, pars, ds_only = TRUE) {
  t <- system.time({
    a <- pars$assay # assay (counts, cpm, logcounts or vstresiduals)
    cores = as.numeric(pars$cores)
    
    if (!ds_only) 
      assay(sce, a) <- switch(a, 
                              counts = counts(sce),
                              cpm = calculateCPM(counts(sce)),
                              logcounts = normalizeCounts(computeLibraryFactors(sce)),
                              vstresiduals = vst(counts(sce), show_progress = FALSE)$y)
    
    # create the design of the study:
    samples = metadata(sce)$experiment_info$sample_id
    group = metadata(sce)$experiment_info$group_id
    batch = factor(metadata(sce)$experiment_info$batch_id)
    patient_id = sapply(metadata(sce)$experiment_info$sample_id,
                       function(x){ strsplit(x, ".", fixed = TRUE)[[1]][[1]]})

    design = model.matrix(~group + patient_id + batch)
    # rownames of the design must indicate sample ids:
    rownames(design) = samples

    # by default, using counts, then try all the other ones too:
    res <- tryCatch(distinct_test(sce, 
                                  name_assays_expression = a,
                                  design = design,
				  n_cores = cores,
                                  min_non_zero_cells = 1), error = function(e) e)
    
    res$gene = as.character(res$gene)
    res$cluster_id = as.character(res$cluster_id)
  })[[3]]
  list(rt = t, tbl = res)
}

