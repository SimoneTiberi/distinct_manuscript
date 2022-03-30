rm(list = ls())
setwd("~/Desktop/distinct project/SIMULATED data/EXTRA SIMULATIONS/unbalanced/")

library(ggplot2)
pairing = "paired" 
local = 1 # locally, globally adjusted p-values, 3 overall with filtered
min_sim_mean = 0.2
n_cells = 200
n_rep = 5

source("~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/all_methods.R")
methods = all_methods[c(1:15, 19:22)]

# all_methods = sort(all_methods)
match = match( methods,  all_methods )

colours = c(all_colours[match], "white")
# points, borders:
shape_border = c(0, 1, 2, 5, 6, 3, 4, 8)
# points, fill:
shape_fill = c(15, 16, 17, 23, 25,  3, 4, 8)

normalizations = c("basics", "cpm", "linnorm", "logcounts", "vst", ".counts", "-dream2", "-nbinom")
sel_shape = c()
for(i in 1:length(methods)){
  sel_shape[i] = which(sapply(normalizations, function(norm){
    grepl(norm, methods[i], fixed = TRUE)
  }))
}
#shape_border[sel_shape]
#shape_fill[sel_shape]

results = c( "3_3", "3_2", "4_3", "5_3")

min_sim_mean = 0.2

types = c("de", "dp", "dm", "db", "dv")

n_rep = 5

# OVERALL CURVE ACCOUNTING FOR ALL TYPES:
ROC_ALL = FDR_ALL = res = list()
for(rr in 1:length(results)){
  
  if(results[rr] == "3_3"){ # original simulation location"
    # original results
    res_names = list.files("~/Desktop/distinct project/SIMULATED data/FULL muscat pipeline/results/results_7Samples/", full.names =TRUE)
  }else{ # new additional simulations:
    res_names = list.files(paste0("results_", results[rr]), full.names =TRUE)
  }
  
  for(j in seq_along(methods)){
    sel_method = grep(methods[j], res_names, ignore.case = TRUE)
    
    RES = list()
    for(i in seq_along(types)){
      RES[[i]] = list()
      for(rep in 1:n_rep){
        sel_type = grep(paste0(types[i], "200,", rep), res_names[sel_method])
        
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
    colnames(p_adj) = methods_names[match]
  }
  
  p_adj = data.frame(p_adj)
  colnames(p_adj) = methods_names[match]
  
  E = (res[[1]]$sim_mean.A + res[[1]]$sim_mean.B)/2
  sel_genes = E > min_sim_mean
  
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
  
  ROC_ALL[[rr]] = plot_roc(cobraplot)
  
  FDR_ALL[[rr]] = plot_fdrtprcurve(cobraplot, plottype = c("points"),
                                   pointsize = 0, linewidth = 2) +
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,1)) +
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
    #guides(colour = guide_legend(ncol = 3, byrow = FALSE)) + 
    #geom_point(data = cobraplot@fdrtpr, 
    #           size = 6, alpha = 0.25, 
    #           colour = rep(colours[-length(colours)], 4)) +
    guides(colour = guide_legend(ncol = 2, byrow = FALSE,
                                 override.aes = list(shape = shape_fill[sel_shape],
                                                     fill = colours[-length(colours)]) ) ) +
    geom_point(size = 6, aes_string(fill = "method", colour = "method", shape = "method"), 
               shape = rep(shape_border[sel_shape], 4), stroke = 2, alpha = 1) + # stroke = line width
    geom_point(size = 6, aes_string(fill = "method", colour = "method", shape = "method"), 
               shape = rep(shape_fill[sel_shape], 4), stroke = 1, alpha = 0.25)
}

results

AA = egg::ggarrange( plots = 
                       list(FDR_ALL[[1]] + xlab("") + labs(title = "3 vs. 3") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                            FDR_ALL[[2]] + xlab("") + ylab("") + labs(title = "3 vs. 2") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                            FDR_ALL[[3]] + labs(title = "4 vs. 3") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                            FDR_ALL[[4]] + ylab("") + labs(title = "5 vs. 3") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35))),
                     bottom = ggpubr::get_legend(FDR_ALL[[1]]),
                     ncol = 2, nrow = 2)

ggsave(filename = paste0("FDR-overall-unbalanced.pdf"),
       plot = AA,
       device = "pdf",
       path = "~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/plots muscat simulation - sensitivity analysis/Plots/",
       width = 16,
       height = 22,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
