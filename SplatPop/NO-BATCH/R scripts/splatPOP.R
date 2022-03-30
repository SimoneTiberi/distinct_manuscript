setwd("~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/splatPOP/NO-BATCH/")

suppressPackageStartupMessages({
  library("splatter")
  library("scater")
  library("ggplot2")
  library("VariantAnnotation")
})

library(data.table)
gff <- read.table("../gff-vcf/gencode.v19.annotation.gff3", sep="\t")
dim(gff)
gff <- subset(gff, V3 == "gene")
table(gff[,9])
gene_type = strsplit(gff[,9], ";")
gene_type = sapply(gene_type, function(x){ # 4-th element contains gene type
  x[4]
})
sel_genes = gene_type == "gene_type=protein_coding"
table(sel_genes)
gff = gff[sel_genes,]
dim(gff)
# 20 k genes, ok

# sim parameters:
Loc   = c(0.2, 0.5, 1, 1.5, 0, 0, 0, 0)
Scale = c(0, 0, 0, 0, 0.2, 0.5, 1, 1.5)

cbind(Loc, Scale)

n_sim = length(Loc)

for(i in 1:n_sim){
  # batchCells = n_cells per sample, but across both groups!
  params.group <- newSplatPopParams(eqtl.n = 0, # no eqtl                            
                                    batchCells = 200, # 200 cells per sample (across both groups)
                                    de.facLoc = Loc[i], # DE Location
                                    de.facScale = Scale[i], # DE Scale
                                    de.prob = 0.05,                               
                                    group.prob = c(0.5, 0.5))
  
  # specify n samples per group
  vcf = mockVCF(n.samples = 4)
  
  set.seed(61217)
  sim <- splatPopSimulate(gff = gff, vcf = vcf, params = params.group)
  colnames(sim) = 1:ncol(sim)
  sim
  
  # remove unnecessary objects from sim:
  assays(sim)$BatchCellMeans = NULL
  assays(sim)$BaseCellMeans = NULL
  assays(sim)$BCV = NULL
  assays(sim)$CellMeans = NULL
  assays(sim)$TrueCounts = NULL
  
  # define DE genes:
  mean(rowData(sim)$GroupDE.Group2 != rowData(sim)$GroupDE.Group1)
  rowData(sim)$DE = rowData(sim)$GroupDE.Group2 != rowData(sim)$GroupDE.Group1
  
  sim <- logNormCounts(sim)
  sim <- runUMAP(sim)
  
  plotUMAP(sim, colour_by = "Group", shape_by = "Sample")
  plotUMAP(sim, colour_by = "Batch", shape_by = "Sample")
  plotUMAP(sim, colour_by = "Sample", shape_by = "Group")
  plotUMAP(sim, colour_by = "Sample", shape_by = "Batch")
  
  # at least 10 non-zero cells per gene:
  sim <- sim[rowSums(counts(sim) > 0) >= 10, ]
  
  sim <- computeLibraryFactors(sim)
  sim <- logNormCounts(sim)
  assays(sim)$cpm <- calculateCPM(sim)
  
  # add Linnorm normalized data:
  library(Linnorm)
  Linnorm_data <- Linnorm.Norm( assays(sim)$counts )
  assays(sim)$linnorm = Linnorm_data
  
  # vstresiduals
  library(sctransform)
  assays(sim)$vstresiduals <- suppressWarnings(
    sctransform::vst(counts(sim), show_progress = FALSE)$y)
  
  # basics:
  library(BASiCS)
  
  sim$BatchInfo = sim$Sample
  
  set.seed(61217)
  Chain <- BASiCS_MCMC(sim, WithSpikes = FALSE, Regression = FALSE,
                       N = 10^3,
                       Threads = 8,
                       Thin = 10,
                       Burn = 500)
  
  assays(sim)$basics <- BASiCS_DenoisedCounts(sim, Chain)
  
  library(muscat)
  sim$cluster_id = "1"
  sim$sample_id = paste(sim$Sample, sim$Group, sep = "-")
  sim$group_id = sim$Group
  
  sce <- prepSCE(sim, "cluster_id", "sample_id", "group_id", drop = TRUE) # prep. SCE for `muscat`
  
  name = paste0("results/sce-Loc",Loc[i],"-Scale",Scale[i],".RDS")
  saveRDS(sce, file = name)
}