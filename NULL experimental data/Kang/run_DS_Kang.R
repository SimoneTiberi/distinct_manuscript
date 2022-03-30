R: "/usr/local/R/R-4.0.0/bin/R"
alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.0 /usr/local/R/R-4.0.0/bin/R'
export R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.0 ; /usr/bin/time -v /usr/local/R/R-4.0.0/bin/R -e '.libPaths()'
R

setwd("~/Desktop/distinct project/EXPERIMENTAL data/Kang (muscData)")

rm(list =ls())

library(muscData)
sce = Kang18_8vs8()

# remove multiplets and unassigned cells
sce <- sce[, sce$multiplets == "singlet"] # remove multiplets
sce <- sce[, !is.na(sce$cell)] # remove unassigned cells

# REMOVE 1 samples (potential outlier): 1039
sce <- sce[, sce$ind %in% c(101, 107, 1015, 1016, 1244, 1256, 1488)] # keep control samples only

# Select Vehicle mouse only:
sel = sce$stim == "ctrl"
mean(sel)
sce = sce[, sel ]

# remove genes with < 10 non-zero cells (across ALL clusters).
sel_genes = rowSums( assays(sce)$counts>0 ) >= 20
mean(sel_genes)
sce = sce[ sel_genes , ]

sce
colData(sce)

# 8 patients, info abuot individuals stored here:
cd = colData(sce)
cd
table(cd$ind) # 8 individuals
table(cd$stim) # 2 experimental conditions: control (ctrl) and stimulated (stim)
table(cd$ind, cd$stim)
table(cd$cluster) # 10 clusters
table(cd$cell) # 8 cell types
table(cd$cluster, cd$cell) # cluster is not cell

individuals_selected = unique(cd$ind)

# make 3 random splits of samples in 2 groups:
n_groups = 3

n_samples = length(individuals_selected); n_samples
COMBN = combn(n_samples, 3)

set.seed(123456)
sel = sample.int(ncol(COMBN), n_groups)
GROUPS = lapply(1:3, 
                function(i){
                  factor(ifelse(seq_len(n_samples) %in%  COMBN[,sel[i] ], "A", "B"))
                })
GROUPS

sce_counts = assays(sce)$counts

##### CHOOSE THE CLUSTER HERE:
clust = factor(cd$cell)

library(scran)
library(scater)

sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)
assays(sce)$cpm <- calculateCPM(sce)

# Linnorm:
library(Linnorm)
Linnorm_data <- Linnorm.Norm( assays(sce)$counts )
assays(sce)$linnorm = Linnorm_data

# vst:
vstresiduals = sctransform::vst( as.matrix(sce_counts) )$y
assays(sce)$vstresiduals = vstresiduals

# basics:
library(BASiCS)

sce$BatchInfo = sce$ind

Chain <- BASiCS_MCMC(sce, WithSpikes = FALSE, Regression = FALSE,
                     N = 10^3,
                     Thin = 10,
                     Burn = 500)

assays(sce)$basics <- BASiCS_DenoisedCounts(sce, Chain)

#save(sce, file = "results/sce.RData")
#load("results/sce.RData")

# store norm counts and remove sce object:
logcounts = as.matrix(assays(sce)$logcounts)
rownames(logcounts) = seq_len(nrow(logcounts))

cpm = as.matrix(assays(sce)$cpm)
rownames(cpm) = seq_len(nrow(cpm))

vstresiduals = as.matrix(assays(sce)$vstresiduals)
rownames(vstresiduals) = seq_len(nrow(vstresiduals))

linnorm = as.matrix(assays(sce)$linnorm)
rownames(linnorm) = seq_len(nrow(linnorm))

basics = as.matrix(assays(sce)$basics)
rownames(basics) = seq_len(nrow(basics))

rm(sce)

###############################################################################################
# define PB methods:
###############################################################################################
pb <- dplyr::bind_rows(
  expand.grid(
    stringsAsFactors = FALSE,
    assay =  "counts", fun = "sum", scale = FALSE,
    method = c("edgeR", "limma-voom"),
    treat = c(FALSE)
  ),
  expand.grid(
    stringsAsFactors = FALSE, scale = FALSE,
    assay = c("logcounts", "vstresiduals", "linnorm", "basics"),
    fun = "mean", method = c("limma-trend")
  ),
  expand.grid(
    stringsAsFactors = FALSE, scale = FALSE,
    assay = c("linnorm", "basics"),
    fun = "mean", method = c("edgeR")
  ),
  data.frame(
    stringsAsFactors = FALSE, scale = TRUE,
    assay = "cpm", fun = "sum", method = c("edgeR", "limma-trend"),
    treat = c(FALSE)
  )
)
pb$treat[is.na(pb$treat)] <- FALSE
pb$id <- with(pb, sprintf("%s%s.%s.%s%s",
                          method, ifelse(treat, "-treat", ""),
                          fun, ifelse(scale, "scale", ""), assay))

###############################################################################################
# run:
###############################################################################################

res_distinct = RES_PB = list()

library(SingleCellExperiment)

for(g in 1:n_groups){
  groups = GROUPS[[g]]
  
  experiment_info_sel = data.frame(sample_id = factor(individuals_selected),
                                   group_id = groups ) # randomly separate samples in two groups (with prob 0.5)
  
  matching = match(cd$ind, experiment_info_sel$sample_id)
  head(cd$ind); head(experiment_info_sel$sample_id[matching])
  group_id = experiment_info_sel$group_id[matching]
  
  colData_sel = list(sample_id = factor(cd$ind),
                     cluster_id = clust,
                     group_id = factor(group_id))
  
  x <- SingleCellExperiment(assays = list(counts = sce_counts,
                                          cpm = cpm,
                                          logcounts = logcounts,
                                          vstresiduals = vstresiduals,
                                          linnorm = linnorm,
                                          basics = basics),
                            colData = colData_sel,
                            metadata = list(experiment_info = experiment_info_sel,
                                            n_cells = table(colData_sel$sample_id)))
  
  print(paste("sce defined - ", g))
  
  # WITHOUT as.matrix(sce_counts)) aggregateData fails!
  
  RES_PB[[g]] = list()
  # PB methods:
  library(muscat)
  
  for(i in 1:nrow(pb)){
    assay = pb$assay[i]
    fun = pb$fun[i]
    scale = pb$scale[i]
    method = pb$method[i]
    
    PB <- aggregateData(x,
                        assay = assay, fun = fun,
                        scale = scale,
                        by = c("cluster_id", "sample_id"))
    
    set.seed(61217)
    try({ res <- pbDS(PB, method = method,
                      min_cells = 1,
                      verbose = FALSE, filter = "none")
    RES = do.call(rbind, res$table[[1]])
    
    RES_PB[[g]][[i]] = RES
    }, TRUE)
  }
  names(RES_PB[[g]]) = pb$id
  
  print(paste("PB done - ", g))
  
  ###############################################################################################
  # run distinct:
  ###############################################################################################
  assays = c("cpm", "logcounts", "vstresiduals", "linnorm", "basics")
  
  design_distinct = model.matrix(~metadata(x)$experiment_info$group_id)
  rownames(design_distinct) = metadata(x)$experiment_info$sample_id
  
  res_distinct[[g]] = list()
  
  library(distinct)
  for(i in 1:length(assays)){
    assay = assays[i]
    
    set.seed(61217)
    res = distinct_test(x, 
                        name_assays_expression = assay,
                        design = design_distinct,
                        min_non_zero_cells = 10,
                        n_cores = 8)
    res_distinct[[g]][[i]] = res
    
    print(paste("distinct done - ", assay))
  }
  names(res_distinct[[g]]) = assays
}

save(GROUPS, clust, n_groups, res_distinct, RES_PB, sce_counts, file = paste0("results/res_cell_clust_7Samples.RData") )
