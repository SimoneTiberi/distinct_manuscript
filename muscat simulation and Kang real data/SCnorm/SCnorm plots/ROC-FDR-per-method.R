rm(list = ls())
setwd("~/Desktop/distinct project/SIMULATED data/FULL muscat pipeline")


library(ggplot2)
pairing = "paired" 
local = 1 # locally, globally adjusted p-values, 3 overall with filtered
min_sim_mean = 0.2
n_cells = 200
n_rep = 5

# NEW RESULTS:
res_names = list.files(paste0("results/results_7Samples"), full.names =TRUE)
file.exists(res_names)

source("~/Desktop/distinct project/distinct Article/scripts/v4 - SCnorm only/SCnorm plots/all_methods.R")
METHODS = c("distinct", "edgeR", "limma", "scDD")

min_sim_mean = 0.2
types = c("de", "dp", "dm", "db", "dv")

REP = list(
  c(4),
  c(1,2),
  c(1,3,4),
  c(1:3),
  c(4)
)

# points, borders:
shape_border = c(0, 1, 2, 5, 6, 3, 4, 8, 11)
# points, fill:
shape_fill = c(15, 16, 17, 23, 25,  3, 4, 8, 11)

#shape_border[sel_shape]
#shape_fill[sel_shape]

FDR_ALL = ROC_ALL = list()
for(mm in seq_along(METHODS)){
  methods = all_methods[grep(METHODS[mm], all_methods, ignore.case = TRUE)]
  
  match = match( methods,  all_methods )
  colours = c(all_colours[match], "white")
  
  normalizations = c("basics", "cpm", "linnorm", "logcounts", "vst", ".counts", "-dream2", "-nbinom", "SCnorm")
  sel_shape = c()
  for(i in 1:length(methods)){
    sel_shape[i] = which(sapply(normalizations, function(norm){
      grepl(norm, methods[i], fixed = TRUE)
    }))
  }
  
  # OVERALL CURVE ACCOUNTING FOR ALL TYPES:
  res = list()
  for(j in seq_along(methods)){
    sel_method = grep(methods[j], res_names, ignore.case = TRUE)
    
    RES = list()
    for(i in seq_along(types)){
      RES[[i]] = list()
      for(rep in 1:n_rep){
        if(rep %in% REP[[i]]){
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
  
  ROC_ALL[[mm]] = plot_roc(cobraplot)
  
  FDR_ALL[[mm]] = plot_fdrtprcurve(cobraplot, plottype = c("points"),
                                   pointsize = 0, linewidth = 2) +
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
    geom_point(size = 6, aes_string(fill = "method", colour = "method", shape = "method"), 
               shape = rep(shape_border[sel_shape], 4), stroke = 2, alpha = 1) + # stroke = line width
    geom_point(size = 6, aes_string(fill = "method", colour = "method", shape = "method"), 
               shape = rep(shape_fill[sel_shape], 4), stroke = 1, alpha = 0.25)
}

########################################################################################
# get common legend:
########################################################################################
methods = all_methods[ unlist(sapply(METHODS, grep, all_methods, ignore.case = TRUE)) ] 
match = match( methods,  all_methods )
n_methods = length(methods)

normalizations = c("basics", "cpm", "linnorm", "logcounts", "vst", ".counts", "-dream2", "-nbinom", "SCnorm")
sel_shape = c()
for(i in 1:length(methods)){
  sel_shape[i] = which(sapply(normalizations, function(norm){
    grepl(norm, methods[i], fixed = TRUE)
  }))
}

padj = matrix(runif(100), ncol = n_methods, nrow = 100)
rownames(padj) = paste0("F", 1:100)
colnames(padj) = methods_names[match]
padj = as.data.frame(padj)
head(padj)

library(iCOBRA)
cobra = COBRAData(pval = padj,
                  padj = padj,
                  truth = data.frame(status = round(runif(100)),
                                     row.names = paste0("F", 1:100)),
                  object_to_extend = NULL)

cobraperf <- calculate_performance(cobra, binary_truth = "status",
                                   thrs = c(0.01, 0.05, 0.1, 0.2))

colours = c(all_colours[match], "white")

cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colours,
                                   facetted = TRUE, incloverall = FALSE,
                                   conditionalfill = FALSE)

gg = plot_fdrtprcurve(cobraplot, plottype = c("points"),
                      pointsize = 0, linewidth = 2) +
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
  guides(colour = guide_legend(ncol = 2, byrow = FALSE,
                               override.aes = list(shape = shape_fill[sel_shape],
                                                   fill = colours[-length(colours)]) ) ) +
  geom_point(size = 6, aes_string(fill = "method", colour = "method", shape = "method"), 
             shape = rep(shape_border[sel_shape], 4), stroke = 2, alpha = 1) + # stroke = line width
  geom_point(size = 6, aes_string(fill = "method", colour = "method", shape = "method"), 
             shape = rep(shape_fill[sel_shape], 4), stroke = 1, alpha = 0.25) +
  theme(legend.key.height=unit(2.25,"line"))

ggpubr::get_legend(gg)


########################################################################################
# increase xlim:
########################################################################################
for(i in 1:4){
  FDR_ALL[[i]] = FDR_ALL[[i]] + 
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,1))
}


########################################################################################
# make common plot:
########################################################################################
library(grid) # for unit()
library(ggpubr)


AA = egg::ggarrange( plots = 
                       list(FDR_ALL[[1]] + xlab("") + labs(title = METHODS[1]) + theme(legend.position = "null", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                            FDR_ALL[[2]] + xlab("") + ylab("") + labs(title = METHODS[2]) + theme(legend.position = "null", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                            FDR_ALL[[3]] + labs(title = METHODS[3]) + xlab("") + theme(legend.position = "null", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                            FDR_ALL[[4]] + ylab("") + labs(title = METHODS[4]) + theme(legend.position = "null", plot.title = element_text(hjust = 0.5, face = "bold", size=35))),
                     bottom = ggpubr::get_legend(gg),
                     ncol = 2, nrow = 2)
AA

ggsave(filename = paste0("FDR-by-method-", pairing, "-SCnorm.pdf"),
       plot = AA,
       device = "pdf",
       path = "~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/plots muscat simulation - NEW normalization/Plots/",
       width = 16,
       height = 23,
       units = "in",
       dpi = 300,
       limitsize = TRUE)


