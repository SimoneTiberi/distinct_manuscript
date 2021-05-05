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

CELLS = c(50, 100, 200, 400)
# use a different threshold for the min counts across cells:
min_sim_mean = c()
min_sim_mean[50] = 0.2
min_sim_mean[100] = 0.2
min_sim_mean[200] = 0.2
min_sim_mean[400] = 0.2

types = c("de", "dp", "dm", "db", "dv")

n_rep = 5

# OVERALL CURVE ACCOUNTING FOR ALL TYPES:
ROC_ALL = FDR_ALL = res = cobraperf_list = list()
for(j in seq_along(methods)){
  sel_method = grep(methods[j], res_names, ignore.case = TRUE)
  
  for(n_cells in CELLS){
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
    res[[n_cells]] = do.call(rbind, RES)
    
    # FILTER GENES THAT RESPECT THE MIN E REQUIREMENT ACROSS ALL:
    E = (res[[n_cells]]$sim_mean.A + res[[n_cells]]$sim_mean.B)/2
    sel_genes = E > min_sim_mean[n_cells]
    
    p_val = data.frame(res[[n_cells]]$p_val[sel_genes])
    colnames(p_val) = factor(n_cells, levels = CELLS)
    
    if(local == 1){
      # get locally adjusted p-vals
      p_adj = data.frame(res[[n_cells]]$p_adj.loc[sel_genes])
    }
    if(local == 2){
      # get globally adjusted p-vals
      p_adj = data.frame(res[[n_cells]]$p_adj.glb[sel_genes])
    }
    if(local == 3){
      p_adj = apply(p_val, 2, p.adjust, method = "BH")
    }
    p_adj = data.frame(p_adj)
    colnames(p_adj) = factor(n_cells, levels = CELLS)
    
    # set NA's to 1:
    p_val[is.na(p_val)] = 1
    p_adj[is.na(p_adj)] = 1
    
    library(iCOBRA)
    cobra = COBRAData(pval = p_val,
                      padj = p_adj,
                      truth = data.frame(status = res[[n_cells]]$is_de[sel_genes]),
                      object_to_extend = NULL)
    
    cobraperf_list[[n_cells]] <- calculate_performance(cobra, binary_truth = "status",
                                                       thrs = c(0.01, 0.05, 0.1, 0.2))
    
  }
  
  cobraperf = cobraperf_list[[n_cells]]
  
  cobraperf@fdrtpr = do.call(rbind,
                             lapply(CELLS, function(n_cells){
                               cobraperf_list[[n_cells]]@fdrtpr
                             }))
  
  cobraperf@fdrtpr = do.call(rbind,
                             lapply(CELLS, function(n_cells){
                               cobraperf_list[[n_cells]]@fdrtpr
                             }))
  
  cobraperf@fdrtprcurve = do.call(rbind,
                                  lapply(CELLS, function(n_cells){
                                    cobraperf_list[[n_cells]]@fdrtprcurve
                                  }))
  
  cobraperf@fdrnbr = do.call(rbind,
                             lapply(CELLS, function(n_cells){
                               cobraperf_list[[n_cells]]@fdrnbr
                             }))
  
  cobraperf@fdrnbrcurve = do.call(rbind,
                                  lapply(CELLS, function(n_cells){
                                    cobraperf_list[[n_cells]]@fdrnbrcurve
                                  }))
  
  cobraperf@tpr = do.call(rbind,
                          lapply(CELLS, function(n_cells){
                            cobraperf_list[[n_cells]]@tpr
                          }))
  
  cobraperf@fpr = do.call(rbind,
                          lapply(CELLS, function(n_cells){
                            cobraperf_list[[n_cells]]@fpr
                          }))
  
  cobraperf@roc = do.call(rbind,
                          lapply(CELLS, function(n_cells){
                            cobraperf_list[[n_cells]]@roc
                          }))
  
  cobraperf@overlap = data.frame()
  
  library(RColorBrewer);
  colours = c(brewer.pal(5, "Blues")[2:5])
  
  cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colours,
                                     facetted = TRUE, incloverall = FALSE,
                                     conditionalfill = FALSE)
  
  # redefine levels of cells so that 50 is first:
  cobraplot@roc$method = factor(cobraplot@roc$method, levels = CELLS)
  
  ROC_ALL[[j]] = plot_roc(cobraplot)
  
  # redefine levels of cells so that 50 is first:
  cobraplot@fdrtpr$method = factor(cobraplot@fdrtpr$method, levels = CELLS)
  cobraplot@fdrtprcurve$method = factor(cobraplot@fdrtprcurve$method, levels = CELLS)
  
  FDR_ALL[[j]] = plot_fdrtprcurve(cobraplot, plottype = c("points"),
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
    guides(colour = guide_legend(ncol = 4, byrow = FALSE)) + 
    geom_point(size = 6, alpha = 0.25, 
               colour = rep(colours, each = 4))
  #geom_point(data = cobraplot@fdrtpr, 
  #           size = 6, alpha = 0.25, 
  #           colour = rep(colours[-length(colours)], 3))
}

# library(ggpubr)
#print( ggpubr::ggarrange(FDR_ALL[[1]], 
#                 FDR_ALL[[2]],
#                 FDR_ALL[[3]], 
#                 FDR_ALL[[4]],
#                 FDR_ALL[[5]], 
#                 FDR_ALL[[6]],
#                 FDR_ALL[[7]], 
#                 FDR_ALL[[8]],
#                 FDR_ALL[[9]], 
#                 FDR_ALL[[10]],
#                 FDR_ALL[[11]],
#                 labels = methods_names,
#                 font.label = list(size = 35),
#                 ncol = 3, nrow = 4,
#                 legend = "bottom", common.legend = TRUE) )

AA = egg::ggarrange( plots = 
                       list(FDR_ALL[[1]] + xlab("") + labs(title = methods_names[1]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=30)),
                            FDR_ALL[[2]] + xlab("") + ylab("") + labs(title = methods_names[2]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=30)),
                            FDR_ALL[[3]] + xlab("") + ylab("") + labs(title = methods_names[3]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=30)),
                            FDR_ALL[[4]] + xlab("") + labs(title = methods_names[4]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=30)),
                            FDR_ALL[[5]] + xlab("") + ylab("") + labs(title = methods_names[5]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=30)),
                            FDR_ALL[[6]] + xlab("") + ylab("") + labs(title = methods_names[6]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=30)),
                            FDR_ALL[[7]] + xlab("") + labs(title = methods_names[7]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=30)),
                            FDR_ALL[[8]] + xlab("") + ylab("") + labs(title = methods_names[8]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=30)),
                            FDR_ALL[[9]] + xlab("") + ylab("") + labs(title = methods_names[9]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=30)),
                            FDR_ALL[[10]] + xlab("") + labs(title = methods_names[10]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=30)),
                            FDR_ALL[[11]] + ylab("") + xlab("") + labs(title = methods_names[11]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=30)),
                            FDR_ALL[[12]] + ylab("") + labs(title = methods_names[12]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=30)) +
                              scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.4, 0.75, 1), limits = c(0,0.5)) ,
                            FDR_ALL[[13]] + labs(title = methods_names[13]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=30)) + 
                              scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.7, 1), limits = c(0,0.75)),
                            FDR_ALL[[14]] + labs(title = methods_names[14]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=30)) + 
                              scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.7, 1), limits = c(0,0.45))),
                     bottom = ggpubr::get_legend(FDR_ALL[[1]]),
                     ncol = 3, nrow = 5)

if(saving){
  ggsave(filename = paste0("FDR-overall-n_cells-by_method.pdf"),
         plot = AA,
         device = "pdf",
         path = "~/Desktop/distinct project/distinct Article/v1/images/muscat/n_cells",
         width = 20,
         height = 37.5,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}
