# the raw sequencing data have been deposited in the European Genome-phenome Archive database
# under study accession id EGAS00001002791 and dataset accession id EGAD0000100391031, 
# which are available in FASTQ file format upon request and approval

rm(list =ls())
x = data.table::fread("data/GSE108989_CRC.TCell.S11138.count.txt.gz", header = TRUE)
class(x); dim(x)
# 23,460 11,140
head(colnames(x)); head(rownames(x))
x[1:10, 1:10]

metadata = data.table::fread("metadata/GSE108989_family.xml/GSE108989-tbl-1.txt", header = FALSE)
dim(metadata)
#  11138
metadata
# cell id
# patient id
# major Cluster
# sample Type

# samples:
table(metadata$V2)
P0123 P0215 P0309 P0411 P0413 P0701 P0825 P0909 P1012 P1207 P1212 P1228 
1238   754   758   764   765   852  1253  1105  1065   210  1164  1210 

# samples, major Cluster
table( metadata$V3 )
table( metadata$V2, metadata$V3  )

# samples, sample type:
table( metadata$V4 )
table( metadata$V2, metadata$V4  )


# Store counts, gene ID, sample ID
counts = x[,-c(1,2)]
dim(counts)
# [1] 23459 11138

genes = as.character(x$geneID)
length(genes)
# 23459

matching = match(colnames(counts), metadata$V1)
colnames(counts)[1:10]; metadata$V1[matching][1:10]

samples = metadata$V2[matching]
major_cluster = metadata$V3[matching]
sample_type   = metadata$V4[matching]

table(samples)
table(major_cluster)
table(sample_type)

# save(counts, CPM, samples, genes, major_cluster, sample_type, file = "results/data.RData")

CPM <- calculateCPM(counts)

# PCA/tSNE/UMAP plots based on CPMs:
library(SingleCellExperiment)
sce <- SingleCellExperiment(assays = list(counts = counts,
                                        cpms = CPM),
                          colData = list(sample_id = factor(samples),
                                         major_cluster = factor(major_cluster),
                                         sample_type = factor(sample_type)),
                          metadata = list(experiment_info = data.frame( sample_id = unique(samples)),
                                          n_cells = table(samples)))

library(scater)
sce = computeLibraryFactors(sce)
sce = logNormCounts(sce)

save(sce, file = "results/sce.RData")
