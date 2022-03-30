rm(list = ls())

BATCH = c(FALSE, TRUE)

for(batch in BATCH){
  
  if(batch){
    setwd("~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/splatPOP/BATCH")
  }else{
    setwd("~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/splatPOP/NO-BATCH")
  }
  
  # sim parameters:
  Loc   = c(0.2, 0.5, 1, 1.5, 0, 0, 0, 0)
  Scale = c(0, 0, 0, 0, 0.2, 0.5, 1, 1.5)
  
  n_sim = length(Loc)
  
  library(SingleCellExperiment)
  library(ggplot2)
  
  ROC = FDR = list()
  for(i in 1:n_sim){
    name = paste0("results/sce-Loc",Loc[i],"-Scale",Scale[i],".RDS")
    sce = readRDS(name)
    
    # define DE genes:
    mean(rowData(sce)$GroupDE.Group2 != rowData(sce)$GroupDE.Group1)
    rowData(sce)$DE = rowData(sce)$GroupDE.Group2 != rowData(sce)$GroupDE.Group1
    
    RES = data.frame(gene_id = rownames(sce),
                     status = rowData(sce)$DE)
    
    # distinct:
    name_res = paste0("results/res_distinct-Loc",Loc[i],"-Scale",Scale[i],".RData")
    load(name_res)
    res_distinct = sapply(res_distinct, function(x){
      matchings = match(RES$gene_id, x$gene)
      x$p_val[matchings]
    })
    colnames(res_distinct) = paste0("distinct.", colnames(res_distinct))
    head(res_distinct);
    dim(res_distinct)
    
    # PB:
    name_res = paste0("results/res_PB-Loc",Loc[i],"-Scale",Scale[i],".RData")
    load(name_res)
    res_PB = sapply(res_PB, function(x){
      matchings = match(RES$gene_id, x$gene)
      x$p_val[matchings]
    })
    head(res_PB);
    dim(res_PB)
    
    # remove PB middle part name:
    colnames(res_PB) = sapply(colnames(res_PB), function(x){
      a = strsplit(x, ".", fixed = TRUE)[[1]]
      paste0(a[1], ".", a[3])
    })
    
    # replace "scalecpm" with "cpm":
    colnames(res_PB) = sub("scalecpm", "cpm", colnames(res_PB))
    
    # scDD:
    if(!batch){
      name_res = paste0("results/res_scDD-Loc",Loc[i],"-Scale",Scale[i],".RData")
      load(name_res)
      res_scDD = sapply(res_scDD, function(x){
        matchings = match(RES$gene_id, x$gene)
        x$p_val[matchings]
      })
      colnames(res_scDD) = paste0("scDD-KS.", colnames(res_scDD))
      head(res_scDD);
      dim(res_scDD)
    }
    
    # merge all methods:
    if(batch){
      res = cbind(res_distinct, res_PB)
    }else{
      res = cbind(res_distinct, res_PB, res_scDD)
    }
    res = as.data.frame(res)
    head(res)
    
    # re-oder res columns based on methods names:
    res = res[,order(colnames(res))]
    
    source("~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/all_methods.R")
    methods = colnames(res)
    
    match = match( methods,  methods_names )
    colours = c(all_colours[match], "white")
    
    # points, borders:
    shape_border = c(0, 1, 2, 5, 6, 3, 4, 8)
    # points, fill:
    shape_fill = c(15, 16, 17, 23, 25,  3, 4, 8)
    
    normalizations = c("basics", "cpm", "linnorm", "logcounts", "vst", ".counts", "-dream2", "-nbinom")
    sel_shape = c()
    for(mm in 1:length(methods)){
      sel_shape[mm] = which(sapply(normalizations, function(norm){
        grepl(norm, methods[mm], fixed = TRUE)
      }))
    }
    
    library(iCOBRA)
    cobra = COBRAData(pval = res,
                      truth = RES)
    
    cobra = calculate_adjp(cobra, method = "BH")
    
    cobraperf <- calculate_performance(cobra, binary_truth = "status",
                                       thrs = c(0.01, 0.05, 0.1, 0.2))
    
    cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colours,
                                       facetted = TRUE, incloverall = FALSE,
                                       conditionalfill = FALSE)
    
    ROC[[i]] = plot_roc(cobraplot)
    
    FDR[[i]] = plot_fdrtprcurve(cobraplot, plottype = c("points"),
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
      #guides(colour = guide_legend(ncol = 3, byrow = FALSE)) + 
      #geom_point(data = cobraplot@fdrtpr, 
      #           size = 6, alpha = 0.25, 
      #           colour = rep(colours[-length(colours)], 4)) +
      geom_point(size = 6, aes_string(fill = "method", colour = "method", shape = "method"), 
                 shape = rep(shape_border[sel_shape], 4), stroke = 2, alpha = 1) + # stroke = line width
      geom_point(size = 6, aes_string(fill = "method", colour = "method", shape = "method"), 
                 shape = rep(shape_fill[sel_shape], 4), stroke = 1, alpha = 0.25) + 
      guides(colour = guide_legend(ncol = 2, byrow = FALSE,
                                   override.aes = list(shape = shape_fill[sel_shape],
                                                       fill = colours[-length(colours)]) ) )
  }
  
  title_names = paste("de.facLoc:", Loc[1:4])
  
  AA = egg::ggarrange( plots = 
                         list(FDR[[1]] + xlab("") + labs(title = title_names[1]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                              FDR[[2]] + xlab("") + ylab("") + labs(title = title_names[2]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                              FDR[[3]] + labs(title = title_names[3]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                              FDR[[4]] + ylab("") + labs(title = title_names[4]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35))),
                       bottom = ggpubr::get_legend(FDR[[1]]),
                       ncol = 2, nrow = 2)
  
  if(batch){
    name =  paste0("FDR-splatPOP-BATCH-Differential-Location.pdf")
  }else{
    name =  paste0("FDR-splatPOP-NO-BATCH-Differential-Location.pdf")
  }
  
  ggsave(filename = name,
         plot = AA,
         device = "pdf",
         path = "~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/splatPOP/Figures",
         width = 16,
         height = 20,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
  
  title_names = paste("de.facScale:", Scale[5:8])
  
  AA = egg::ggarrange( plots = 
                         list(FDR[[5]] + xlab("") + labs(title = title_names[1]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                              FDR[[6]] + xlab("") + ylab("") + labs(title = title_names[2]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                              FDR[[7]] + labs(title = title_names[3]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                              FDR[[8]] + ylab("") + labs(title = title_names[4]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35))),
                       bottom = ggpubr::get_legend(FDR[[5]]),
                       ncol = 2, nrow = 2)
  
  if(batch){
    name = paste0("FDR-splatPOP-BATCH-Differential-Scale.pdf")
  }else{
    name = paste0("FDR-splatPOP-NO-BATCH-Differential-Scale.pdf")
  }
  
  ggsave(filename = name,
         plot = AA,
         device = "pdf",
         path = "~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/splatPOP/Figures",
         width = 16,
         height = 20,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
  
  
  title_names = c(paste("de.facLoc:", Loc[1:4]), paste("de.facScale:", Scale[5:8]))[c(1,5,2,6,3,7,4,8)]
  
  AA = egg::ggarrange( plots = 
                         list(FDR[[1]] + xlab("") + labs(title = title_names[1]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                              FDR[[5]] + xlab("") + ylab("") + labs(title = title_names[2]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                              FDR[[2]] + xlab("") + labs(title = title_names[3]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                              FDR[[6]] + xlab("") + ylab("") + labs(title = title_names[4]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                              FDR[[3]] + xlab("") + labs(title = title_names[5]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                              FDR[[7]] + xlab("") + ylab("") + labs(title = title_names[6]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                              FDR[[4]] + labs(title = title_names[7]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                              FDR[[8]] + ylab("") + labs(title = title_names[8]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size=35))),
                       bottom = ggpubr::get_legend(FDR[[1]]),
                       ncol = 2, nrow = 4)
  
  if(batch){
    name =  paste0("FDR-splatPOP-BATCH-ALL.pdf")
  }else{
    name =  paste0("FDR-splatPOP-NO-BATCH-ALL.pdf")
  }
  
  ggsave(filename = name,
         plot = AA,
         device = "pdf",
         path = "~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/splatPOP/Figures",
         width = 15,
         height = 36.25,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}