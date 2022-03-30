# edit config.yaml first row:
R: "/usr/local/R/R-4.0.0/bin/R"
alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.0 /usr/local/R/R-4.0.0/bin/R'
export R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.0 ; /usr/bin/time -v /usr/local/R/R-4.0.0/bin/R -e '.libPaths()'


cd ~/DS/muscat
mkdir ~/DS/muscat/min_50_cells

R

setwd("~/DS/muscat/sensitivity-analyses/library_size/data/sim_data")

files = list.files()
files
file.exists(files)

clusters = paste0("cluster", 1:3)

library(SingleCellExperiment)

for(file in files){
  sce = readRDS(file)
  cluster_gene =  c()
  for(cluster in clusters){
    sce_one_cl = sce[, sce$cluster_id == cluster]
    
    genes <- rownames(sce_one_cl)[rowSums(counts(sce_one_cl) > 0) >= 50]
    
    cluster_gene = c(cluster_gene, paste(cluster,genes,sep = "-"))
  }
  print(length(cluster_gene))
  
  name = paste0("~/DS/muscat/min_40_cells/cluster_gene-", file)
  saveRDS(cluster_gene,name)
}



setwd("~/DS/muscat/normalizations/normalizations/data/sim_data")

files = list.files()
files
file.exists(files)

clusters = paste0("cluster", 1:3)

library(SingleCellExperiment)

for(file in files){
  sce = readRDS(file)
  cluster_gene =  c()
  for(cluster in clusters){
    sce_one_cl = sce[, sce$cluster_id == cluster]
    
    genes <- rownames(sce_one_cl)[rowSums(counts(sce_one_cl) > 0) >= 50]
    
    cluster_gene = c(cluster_gene, paste(cluster,genes,sep = "-"))
  }
  print(length(cluster_gene))
  
  name = paste0("~/DS/muscat/min_40_cells/cluster_gene-", file)
  saveRDS(cluster_gene,name)
}
