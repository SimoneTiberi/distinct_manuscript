library(ggplot2)

rm(list = ls())
setwd("~/Desktop/distinct project/SIMULATED data/FULL muscat pipeline")

saving = TRUE

pairing = c("paired") # "unpaired"
local = TRUE # locally vs globally adjusted p-values
min_sim_mean = 0.2
n_cells = 200
n_rep = 5

data = "kang" # "magl" or "kang"

# NEW RESULTS with 100 cells:
res_names = list.files(paste0("results/results_7Samples/"), full.names =TRUE)

methods = c(
  "distinct.cpm", # distinct.cpm25FALSE, or distinct.cpm25FALSE1000,
  "distinct.logcounts", # distinct.logcounts25FALSE, or distinct.logcounts25FALSE1000,
  "distinct.vstresiduals", # distinct.vstresiduals25FALSE, or distinct.vstresiduals25FALSE1000,
  "edgeR.sum.counts",
  "edgeR.sum.scalecpm",
  "limma-voom.sum.counts",
  "limma-trend.mean.logcounts",
  "limma-trend.mean.vstresiduals",
  "MM-dream2",
  "MM-nbinom",
  "MM-vst"
)

library(RColorBrewer);
all_colours = c(
  brewer.pal(4, "Reds")[4:2],    # distinct cpm, logcounts, vstresid
  brewer.pal(3, "Blues")[3:2],    # edgeR 2 methods
  brewer.pal(4, "Greens")[4:2],   # limma-trend + # limma-voom 2 methods
  brewer.pal(4, "Greys")[4:2]   # limma-trend + # limma-voom 2 methods
)
# Purples Oranges Greys

all_methods = c(
  "distinct.cpm", # distinct.cpm25FALSE, or distinct.cpm25FALSE1000,
  "distinct.logcounts", # distinct.logcounts25FALSE, or distinct.logcounts25FALSE1000,
  "distinct.vstresiduals", # distinct.vstresiduals25FALSE, or distinct.vstresiduals25FALSE1000,
  "edgeR.sum.counts",
  "edgeR.sum.scalecpm",
  "limma-voom.sum.counts",
  "limma-trend.mean.logcounts",
  "limma-trend.mean.vstresiduals",
  "MM-dream2",
  "MM-nbinom",
  "MM-vst")

# all_methods = sort(all_methods)
match = match( methods,  all_methods )

# method names:
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
  "MM-vstresiduals")

colours = c(all_colours[match])

gg_roc = gg_fdr = list()
res_names = res_names[ grep(data, res_names) ]

# filter results for n_cells:
res_names = res_names[ grep(n_cells, res_names) ]

types = paste0("nill",n_cells, ",", 1:n_rep)

res = list()
for(i in seq_along(types)){
  
  sel_type = grep(types[i], res_names)
  
  RES = list()
  for(j in seq_along(methods)){
    sel_method = grep(methods[j], res_names[sel_type])
    if(length(sel_method) == 0){
      # is method not found:
      RES[[j]] = NULL
    }else{
      x = readRDS(res_names[sel_type][sel_method])$tbl
      
      E = (x$sim_mean.A + x$sim_mean.B)/2
      sel_genes = E > min_sim_mean
      x = x[sel_genes,]
      RES[[j]] = data.frame( p_val = x$p_val, 
                             p_adj.loc = x$p_adj.loc,
                             p_adj.glb = x$p_adj.glb,
                             method = methods_names[j],
                             pairing = pairing[i],
                             replicate = i)
    }
  }
  
  res[[i]] = do.call(rbind, RES)
}
RES = do.call(rbind, res)
head(RES); tail(RES)

.prettify <- function(theme = NULL, ...) {
  if (is.null(theme)) theme <- "classic"
  base <- paste0("theme_", theme)
  base <- getFromNamespace(base, "ggplot2")
  base(base_size = 8) + theme(
    aspect.ratio = 1,
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.2, color = "lightgrey"),
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text = element_text(color = "black"),
    legend.key.size = unit(2, "mm"),
    strip.background = element_rect(fill = NA),
    plot.margin = unit(rep(1, 4), "mm"),
    panel.spacing = unit(0, "mm"),
    legend.margin = margin(0,0,1,0,"mm"),
    ...)}

RES$replicate = factor(RES$replicate)


p <- ggplot(data = RES, aes(x = p_val, y = ..ndensity.., 
                            col = method, fill = method, lty = replicate)) +
  facet_wrap(~ method, ncol = 4) + 
  geom_density(adjust = 0.5, size = 0.3, alpha = 0.1) + 
  scale_alpha_manual(values = c("TRUE" = 0.1, "FALSE" = 0.4)) +
  guides(col = FALSE,
         lty = guide_legend(ncol = 1, order = 2),
         fill = guide_legend(ncol = 2, order = 1,
                             override.aes = list(alpha = 1, col = NA))) +
  scale_x_continuous("p-value", breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
  scale_y_continuous("density", breaks = c(0, 1), expand = c(0, 0.1)) +
  .prettify("bw") + 
  theme(aspect.ratio = 0.6,
        #legend.box.just = "left",
        panel.grid = element_blank(),
        panel.spacing = unit(1, "mm"),
        panel.border = element_rect(color = "grey"),
        strip.text = element_text(size = 4),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(1.5)),
        axis.text.y=element_text(size=rel(1.5)),
        axis.title.y = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(1.5)),
        axis.title.x.top = element_text(size=rel(1.5)),
        legend.title=element_text(size=rel(1.5)),
        legend.text=element_text(size=rel(1.5)),
        legend.key.width=unit(0.4, "cm"),
        legend.position="bottom",
        legend.direction = "horizontal",
        legend.box="horizontal",
        strip.text.x = element_text(size = rel(2)),
        legend.margin=margin()) +
  scale_colour_manual(values = colours) +
  scale_fill_manual(values = colours)
p


if(saving){
  ggsave(filename = paste0("Muscat-null.pdf"),
         plot = last_plot(),
         device = "pdf",
         path = "~/Desktop/distinct project/distinct Article/v1/images/muscat",
         width = 8,
         height = 6,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}
