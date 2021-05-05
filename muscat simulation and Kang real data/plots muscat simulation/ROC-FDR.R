rm(list = ls())
setwd("~/Desktop/distinct project/SIMULATED data/FULL muscat pipeline")

library(ggplot2)
saving = TRUE
pairing = "paired" 
local = 1 # locally, globally adjusted p-values, 3 overall with filtered
min_sim_mean = 0.2
n_cells = 200
n_rep = 5

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
  ",permscdd.logcounts",
  ",permscdd.vstresiduals"
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

gg_roc = gg_fdr = list()
# filter results for n_cells:
res_names = res_names[ grep(n_cells, res_names) ]

types = c("de", "dp", "dm", "db", "dv")

RES_keep = list()
for(i in seq_along(types)){
  RES = res = list()
  
  for(j in seq_along(methods)){
    sel_method = grep(methods[j], res_names, ignore.case = TRUE)
    
    RES[[j]] = list()
    for(rep in 1:n_rep){
      sel_type = grep(paste0(types[i], n_cells, ",", rep), res_names[sel_method], ignore.case = TRUE)
      
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
      
      # for distinct always use GLOBALLY ADJUSTED p-values:
      RES[[j]][[rep]] = data.frame( p_val = x$p_val, 
                                    p_adj.loc = x$p_adj.loc,
                                    p_adj.glb = x$p_adj.glb,
                                    sim_mean.A = x$sim_mean.A,
                                    sim_mean.B = x$sim_mean.B,
                                    is_de = x$is_de)
    }
    res[[j]] = do.call(rbind, RES[[j]])
  }
  # get p-vals
  p_val = sapply(res, function(x){
    E = (x$sim_mean.A + x$sim_mean.B)/2
    sel_genes = E > min_sim_mean
    x$p_val[sel_genes]
  })
  p_val = data.frame(p_val)
  colnames(p_val) = methods_names[match]
  
  if(local == 1){
    # get locally adjusted p-vals
    p_adj = sapply(res, function(x){
      E = (x$sim_mean.A + x$sim_mean.B)/2
      sel_genes = E > min_sim_mean
      x$p_adj.loc[sel_genes]
    })
  }
  if(local == 2){
    # get globally adjusted p-vals
    p_adj = sapply(res, function(x){
      E = (x$sim_mean.A + x$sim_mean.B)/2
      sel_genes = E > min_sim_mean
      x$p_adj.glb[sel_genes]
    })
  }
  if(local == 3){
    p_adj = apply(p_val, 2, p.adjust, method = "BH")
  }
  
  p_adj = data.frame(p_adj)
  colnames(p_adj) = methods_names[match]
  
  E = (res[[1]]$sim_mean.A + res[[1]]$sim_mean.B)/2
  sel_genes = E > min_sim_mean
  
  # set NA's to 1:
  p_val[is.na(p_val)] = 1
  p_adj[is.na(p_adj)] = 1
  
  RES_keep[[i]] = data.frame(truth = res[[1]]$is_de[sel_genes], p_val)
  
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
  
  gg_roc[[i]] = plot_roc(cobraplot)
  
  gg_fdr[[i]] = plot_fdrtprcurve(cobraplot, plottype = c("points"),
                                 pointsize = 6, linewidth = 2) +
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.7)) +
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
    geom_point(data = cobraplot@fdrtpr, 
               size = 6, alpha = 0.25, 
               colour = rep(colours[-length(colours)], 4))
}

library(ggpubr)
legend <- ggpubr::get_legend(gg_fdr[[1]] +
                               guides(colour = guide_legend(ncol = 1, byrow = FALSE)) + 
                               theme(legend.key.size = unit(0.2, "cm")))

ggpubr::as_ggplot(legend)

AA = egg::ggarrange( plots = 
                       list(gg_fdr[[1]] + xlab("") + labs(title = "DE") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                            gg_fdr[[2]] + xlab("") + ylab("") + labs(title = "DP") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                            gg_fdr[[3]] + labs(title = "DM") + xlab("") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                            gg_fdr[[4]] + ylab("") + labs(title = "DB") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                            gg_fdr[[5]] + labs(title = "DV") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                            ggpubr::as_ggplot(legend)),
                     ncol = 2, nrow = 3)

if(saving){
  ggsave(filename = paste0("FDR-by-type-", pairing, ".pdf"),
         plot = AA,
         device = "pdf",
         path = "~/Desktop/distinct project/distinct Article/v1/images/muscat",
         width = 16,
         height = 25,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}
