my_plot_fdrcurve = function (cobraplot,
                             title = "", 
                             stripsize = 15,
                             titlecol = "black",
                             pointsize = 5,
                             xaxisrange = c(0, 1),
                             yaxisrange = c(0, 1),
                             plottype = c("curve", "points"),
                             linewidth = 1,
                             asp = "TPR",
                             shape_fill,
                             shape_border) 
{
  library(ggplot2)
  
  plot_data_lines <- iCOBRA:::fdrtprcurve(cobraplot)
  plot_data_points <- iCOBRA:::fdrtpr(cobraplot)
  xasp <- "FDR"
  yasp <- aspc
  
  pp <- ggplot()
  thresholds <- sort(unique(as.numeric(gsub("thr", "", plot_data_points$thr))))
  plot_data_points$method2.satis <- paste0(plot_data_points$method, 
                                           plot_data_points$satis)
  
  nlevs <- length(unique(plot_data_lines$method))
  nthr <- length(unique(plot_data_points$thr))
  
  pp <- ggplot(plot_data_points, aes_string(x = xasp, y = yasp, 
                                            group = "method")) + 
    geom_vline(xintercept = seq(0, xaxisrange[2], 0.1), colour = "lightgrey", linetype = "dashed") + 
    geom_vline(xintercept = thresholds, linetype = "dashed") + 
    geom_path(size = linewidth, aes_string(colour = "method", linetype = "method")) +
    scale_linetype_manual(values = rep("solid", nlevs), guide = "none") + 
    geom_point(size = pointsize, aes_string(fill = "method", colour = "method", shape = "method"), 
               shape = rep(shape_border, 4), stroke = 2, alpha = 1, linewidth = 1) + # stroke = line width
    geom_point(size = pointsize, aes_string(fill = "method", colour = "method", shape = "method"), 
               shape = rep(shape_fill, 4), stroke = 1, alpha = 0.25) +
    #scale_shape_manual(values = rep(21, nthr), guide = "none") + 
    scale_fill_manual(values = plotcolors(cobraplot), guide = "none", name = "") + 
    scale_color_manual(values = plotcolors(cobraplot), name = "", limits = force) + 
    ylim(ifelse(yasp == "TPR", yaxisrange[1], 0), ifelse(yasp == "TPR", yaxisrange[2], max(c(0, plot_data_lines[[yasp]][plot_data_lines[[xasp]] <= 
                                                                                                                          xaxisrange[2]])))) + 
    scale_x_continuous(breaks = c(thresholds, seq(0, xaxisrange[2], 0.1)), labels = c(thresholds, "", seq(0, xaxisrange[2], 0.1)[-1]), limits = c(xaxisrange[1], 
                                                                                                                                                  xaxisrange[2])) + 
    iCOBRA:::plot_theme(stripsize = stripsize, titlecol = titlecol) + 
    ggtitle(title)
  
  if (isTRUE(facetted(cobraplot))) {
    if (!is.finite(maxsplit(cobraplot))) {
      if (length(plot_data_lines) != 0) 
        msp <- length(unique(plot_data_lines$splitval))
      else msp <- length(unique(plot_data_points$splitval))
    }
    else {
      msp <- maxsplit(cobraplot)
    }
    pp + facet_wrap(~splitval, nrow = ceiling((msp + 1)/3))
  }
  else {
    pp
  }
}