rm(list = ls())
library('SampleQC')
library(muscData)
library(muscat)

sce = Kang18_8vs8()

sce <- sce[, sce$multiplets == "singlet"] # remove multiplets
sce <- sce[, !is.na(sce$cell)] # remove unassigned cells

sce$cluster_id = sce$cell
sce$sample_id = paste(sce$stim, sce$ind, sep = "_")
sce$group_id = sce$stim
# colData(sce) = colData(sce)[, -c(1:5)]
sce <- sce[rowSums(counts(sce) > 1) > 20, ]

library(scater)
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)
assays(sce)$cpm <- calculateCPM(sce)

# add experimental_info via muscat::prepSCE function
sce <- prepSCE(sce, "cluster_id", "sample_id", "group_id", FALSE) # prep. SCE for `muscat`

sce$cluster_id = sce$ind
pb <- aggregateData(sce,
                    assay = "cpm", fun = "sum",
                    by = c("cluster_id", "sample_id"))

library(edgeR)
# plot PB MDS of samples:
pbMDS(pb)

ggsave(filename = "MDS_samples.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "~/Desktop/distinct project/distinct Article/v1/images/Kang/",
       width = 6,
       height = 6,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

# corr plot:
names = names(assays(pb))
colnames = unlist(strsplit(colnames(assays(pb)[[1]]), "_"))[seq(2, 32, 2)]

PB_matrix = assays(pb)[[1]]
for(x in 1:8){
  name = names[x]
  ids = which(name == colnames)
  PB_matrix[,ids] = assays(pb)[[x]][,ids]
}
colnames(PB_matrix) = colnames(assays(pb)[[1]])

library(edgeR);

### female vs male:
group_id = factor(substr(colnames(PB_matrix), 1, 4))
cols = ifelse(group_id == "ctrl", "red", "blue")

sample_id = factor(substr(colnames(PB_matrix), 6, 100))
# 2D MDS plot based on Transcript-level counts:
d_tr = DGEList(counts = PB_matrix, samples = sample_id);
d_tr = calcNormFactors(d_tr)
m3 = plotMDS(d_tr, labels = sample_id, col = cols, main = "MDS")

dev.off() 
# corr plot:
pdf("~/Desktop/distinct project/distinct Article/v1/images/Kang/corrplot_samples.pdf",
    width = 6, height = 6) 

corrplot::corrplot(cor(PB_matrix), method="color" )

dev.off() 
