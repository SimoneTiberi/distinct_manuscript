rm(list =ls())
library('SampleQC')
library(muscData)
library(muscat)

sce = Kang18_8vs8()

sce <- sce[, sce$multiplets == "singlet"] # remove multiplets
sce <- sce[, !is.na(sce$cell)] # remove unassigned cells

sce$cluster_id = sce$cell
sce$sample_id = paste(sce$ind, sce$stim, sep = "_")
sce$group_id = sce$stim
sce$condition = sce$group_id
colData(sce) = colData(sce)[, -c(1:5)]
sce <- sce[rowSums(counts(sce) > 1) > 20, ]

# need to specify sample_id
sce$patient_id = sce$sample_id
sce$condition = sce$group_id

qc_dt   = make_qc_dt(sce)
rownames(sce)

# which QC metrics do we want to use? (the most important bit)
qc_names        = c('log_counts', 'log_feats', 'logit_mito')
# which discrete-valued variables do we want to annotate the samples with?
annot_discrete  = c('well_id', 'patient_id', 'condition')
# which continuous-valued variables do we want to annotate the samples with?
annot_cont      = NULL

mmd_list    = calculate_sample_to_sample_MMDs(qc_dt, qc_names, subsample=200, n_times=20, n_cores=16)
mmd_list    = embed_sample_to_sample_MMDs(mmd_list, qc_dt, annot_discrete, annot_cont, n_nhbrs=5)

plot_mmd_heatmap(mmd_list)

library(ggplot2)
ggsave(filename = "MMD_samples.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "~/Desktop/distinct project/distinct Article/v1/images/Kang/",
       width = 6,
       height = 6,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

