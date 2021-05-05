rm(list = ls())
setwd("~/Desktop/distinct project/SIMULATED data/FULL muscat pipeline")

library(ggplot2)
saving = TRUE
pairing = "paired" 
local = 1 # locally, globally adjusted p-values, 3 overall with filtered

# NEW RESULTS:
res_names = list.files(paste0("results/results_7Samples"), full.names =TRUE)
file.exists(res_names)
methods = c(
  ",distinct.cpm", # distinct.cpm25FALSE, or distinct.cpm25FALSE1000,
  ",distinct.logcounts", # distinct.logcounts25FALSE, or distinct.logcounts25FALSE1000,
  ",distinct.vstresiduals", # distinct.vstresiduals25FALSE, or distinct.vstresiduals25FALSE1000,
  ",edgeR.sum.counts",
  ",edgeR.sum.scalecpm",
  ",limma-voom.sum.counts",
  ",limma-trend.mean.logcounts",
  ",limma-trend.mean.vstresiduals",
  ",MM-dream2",
  ",MM-nbinom",
  ",MM-vst",
  ",scDD.logcounts",
  ",scDD.vstresiduals",
  ",permscdd.logcounts"
  #",permscdd.vstresiduals"
)

library(RColorBrewer);
all_colours = c(
  brewer.pal(4, "Reds")[4:2],    # distinct cpm, logcounts, vstresid
  brewer.pal(3, "Blues")[3:2],    # edgeR 2 methods
  brewer.pal(4, "Greens")[4:2],   # limma-trend 1 method + # limma-voom 2 methods
  brewer.pal(4, "Greys")[4:2],   # MM 3 methods
  brewer.pal(5, "RdPu")[5:2]   # scDD 3 methods
)
# Purples Oranges Greys

all_methods = c(
  ",distinct.cpm", # distinct.cpm25FALSE, or distinct.cpm25FALSE1000,
  ",distinct.logcounts", # distinct.logcounts25FALSE, or distinct.logcounts25FALSE1000,
  ",distinct.vstresiduals", # distinct.vstresiduals25FALSE, or distinct.vstresiduals25FALSE1000,
  ",edgeR.sum.counts",
  ",edgeR.sum.scalecpm",
  ",limma-voom.sum.counts",
  ",limma-trend.mean.logcounts",
  ",limma-trend.mean.vstresiduals",
  ",MM-dream2",
  ",MM-nbinom",
  ",MM-vst",
  ",scDD.logcounts",
  ",scDD.vstresiduals",
  ",permscdd.logcounts",
  ",permscdd.vstresiduals"
)

# methods names:
methods_names = c(
  "distinct.cpm", # distinct.cpm25FALSE, or distinct.cpm25FALSE1000,
  "distinct.log2-cpm", # distinct.logcounts25FALSE, or distinct.logcounts25FALSE1000,
  "distinct.vstresiduals", # distinct.vstresiduals25FALSE, or distinct.vstresiduals25FALSE1000,
  "edgeR.counts",
  "edgeR.cpm",
  "limma-voom.counts",
  "limma-trend.log2-cpm",
  "limma-trend.vstresiduals",
  "MM-dream2",
  "MM-nbinom",
  "MM-vstresiduals",
  "scDD-KS.log2-cpm",
  "scDD-KS.vstresiduals",
  "scDD-perm.log2-cpm",
  "scDD-perm.vstresiduals"
)

# all_methods = sort(all_methods)
match = match( methods,  all_methods )

colours = c(all_colours[match], "white")

CELLS = c(50, 100, 200, 400)
min_sim_mean = c()
min_sim_mean[50] = 0.2
min_sim_mean[100] = 0.2
min_sim_mean[200] = 0.2
min_sim_mean[400] = 0.2

types = c("de", "dp", "dm", "db", "dv")

n_rep = 5

# OVERALL CURVE ACCOUNTING FOR ALL TYPES:
ROC_ALL = FDR_ALL = res = list()
for(n_cells in CELLS){
  
  for(j in seq_along(methods)){
    sel_method = grep(methods[j], res_names, ignore.case = TRUE)
    
    RES = list()
    for(i in seq_along(types)){
      RES[[i]] = list()
      for(rep in 1:n_rep){
        sel_type = grep(paste0(types[i], n_cells, ",", rep), res_names[sel_method])
        
        x = readRDS(res_names[sel_method][sel_type])$tbl
        
        x = x[ order(x$is_de, decreasing = TRUE),]
        head(x); tail(x) # DE at the top.
        
        if(i == grep("db", types) ){ # if type is "db", increase the threshold.
          x$sim_mean.A = x$sim_mean.A/5
          x$sim_mean.B = x$sim_mean.B/5
        }
        
        if(i == grep("dv", types) ){ # if type is "dv", increase the threshold.
          # dv also has different names for groups (a typo in the simulation naming of groups).
          x$sim_mean.A = x$A/5
          
          if( j %in% grep("limma", methods) ){
            x$sim_mean.B = x$B.x/5
          }else{
            x$sim_mean.B = x$B/5
          }
        }
        
        RES[[i]][[rep]] = data.frame( p_val = x$p_val, 
                                      p_adj.loc = x$p_adj.loc,
                                      p_adj.glb = x$p_adj.glb,
                                      sim_mean.A = x$sim_mean.A,
                                      sim_mean.B = x$sim_mean.B,
                                      is_de = x$is_de)
      }
      RES[[i]] = do.call(rbind, RES[[i]])
    }
    res[[j]] = do.call(rbind, RES)
  }
  # get p-vals
  p_val = sapply(res, function(x){
    E = (x$sim_mean.A + x$sim_mean.B)/2
    sel_genes = E > min_sim_mean[n_cells]
    x$p_val[sel_genes]
  })
  p_val = data.frame(p_val)
  colnames(p_val) = methods_names[match]
  
  if(local == 1){
    # get locally adjusted p-vals
    p_adj = sapply(res, function(x){
      E = (x$sim_mean.A + x$sim_mean.B)/2
      sel_genes = E > min_sim_mean[n_cells]
      x$p_adj.loc[sel_genes]
    })
  }
  if(local == 2){
    # get globally adjusted p-vals
    p_adj = sapply(res, function(x){
      E = (x$sim_mean.A + x$sim_mean.B)/2
      sel_genes = E > min_sim_mean[n_cells]
      x$p_adj.glb[sel_genes]
    })
  }
  if(local == 3){
    p_adj = apply(p_val, 2, p.adjust, method = "BH")
    colnames(p_adj) = methods_names[match]
  }
  
  p_adj = data.frame(p_adj)
  colnames(p_adj) = methods_names[match]
  
  E = (res[[1]]$sim_mean.A + res[[1]]$sim_mean.B)/2
  sel_genes = E > min_sim_mean[n_cells]
  
  # set NA's to 1:
  p_val[is.na(p_val)] = 1
  p_adj[is.na(p_adj)] = 1
  
  library(iCOBRA)
  cobra = COBRAData(pval = p_val,
                    padj = p_adj,
                    truth = data.frame(status = res[[1]]$is_de[sel_genes]),
                    object_to_extend = NULL)
  
  cobraperf <- calculate_performance(cobra, binary_truth = "status",
                                     thrs = c(0.01, 0.05, 0.1, 0.2))
  
  cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colours,
                                     facetted = TRUE, incloverall = FALSE,
                                     conditionalfill = FALSE)
  
  ROC_ALL[[n_cells]] = plot_roc(cobraplot)
  
  FDR_ALL[[n_cells]] = plot_fdrtprcurve(cobraplot, plottype = c("points"),
                                        pointsize = 6, linewidth = 2) +
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.2)) +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(3)),
          axis.text.y=element_text(size=rel(3)),
          axis.title.y = element_text(size=rel(3)),
          axis.title.x = element_text(size=rel(3)),
          legend.title=element_text(size=rel(2)),
          legend.text=element_text(size=rel(3)),
          legend.key.width=unit(2, "cm"),
          aspect.ratio = 1, legend.position="bottom",
          legend.box="vertical", legend.margin=margin()) +
    guides(colour = guide_legend(ncol = 3, byrow = FALSE)) + 
    #geom_point(size = 6, alpha = 0.25, 
    #           colour = rep(colours[-length(colours)], 3))
    geom_point(data = cobraplot@fdrtpr, 
               size = 6, alpha = 0.25, 
               colour = rep(colours[-length(colours)], 4))
}


AA = egg::ggarrange( plots = 
                       list(FDR_ALL[[50]] + xlab("") + labs(title = "50 cells") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)) +
                              scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.75)),
                            FDR_ALL[[100]] + xlab("") + ylab("") + labs(title = "100 cells") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)) +
                              scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.75)),
                            FDR_ALL[[200]] + labs(title = "200 cells") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)) +
                              scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.75)),
                            FDR_ALL[[400]] + ylab("") + labs(title = "400 cells") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)) +
                              scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.75))),
                     bottom = ggpubr::get_legend(FDR_ALL[[50]] + 
                                                   guides(colour = guide_legend(ncol = 2, byrow = FALSE))),
                     ncol = 2, nrow = 2)


if(saving){
  ggsave(filename = paste0("FDR-overall-n_cells.pdf"),
         plot = AA,
         device = "pdf",
         path = "~/Desktop/distinct project/distinct Article/v1/images/muscat/n_cells",
         width = 16,
         height = 20,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}
