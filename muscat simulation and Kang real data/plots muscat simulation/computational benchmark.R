rm(list = ls())
setwd("~/Desktop/distinct project/SIMULATED data/FULL muscat pipeline")

saving = TRUE
pairing = "paired" 

# NEW RESULTS:
res_names = list.files(paste0("results/results_7Samples"), full.names =TRUE)
file.exists(res_names)

methods = c(
  ",distinct.cpm", # distinct.cpm25FALSE, or distinct.cpm25FALSE1000,
  ",distinct.logcounts", # distinct.logcounts25FALSE, or distinct.logcounts25FALSE1000,
  ",distinct.vstresiduals", # distinct.vstresiduals25FALSE, or distinct.vstresiduals25FALSE1000,
  ",edgeR.sum.counts",
  ",edgeR.sum.scalecpm",
  ",limma-trend.mean.logcounts",
  ",limma-trend.mean.vstresiduals",
  ",limma-voom.sum.counts",
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
  "limma-trend.log2-cpm",
  "limma-trend.vstresiduals",
  "limma-voom.counts",
  "MM-dream2",
  "MM-nbinom",
  "MM-vstresiduals",
  "scDD-KS.log2-cpm",
  "scDD-KS.vstresiduals",
  "scDD-perm.log2-cpm",
  "scDD-perm.vstresiduals"
)

library(RColorBrewer);
all_colours = c(
  brewer.pal(4, "Reds")[4:2],    # distinct cpm, logcounts, vstresid
  brewer.pal(3, "Blues")[3:2],    # edgeR 2 methods
  brewer.pal(4, "Greens")[4:2],   # limma-trend 1 method + # limma-voom 2 methods
  brewer.pal(4, "Greys")[4:2],   # MM 3 methods
  brewer.pal(5, "RdPu")[5:2]   # scDD 3 methods
)

# filter results for n_cells:
n_cells = 200
res_names = res_names[ grep(n_cells, res_names) ]

n_rep = 5
types = paste0( c("de", "dp", "dm", "db", "dv", "nill") ,n_cells, ",", rep(1:n_rep, each= 5))

TIMES = matrix(NA, nrow = length(methods), ncol = length(types))

apply(TIMES/60, 1, mean)

for(i in seq_along(types)){
  sel_type = grep(types[i], res_names)
  
  for(j in seq_along(methods)){
    sel_method = grep(methods[j], res_names[sel_type], ignore.case = TRUE)
    if(length(sel_method) > 1){
      other_methods = unlist(sapply(methods[-j], grep, res_names[sel_type]))
      sel_method = sel_method[ ! sel_method %in% other_methods]
    }
    TIMES[j,i] = sum(as.numeric(readRDS(res_names[sel_type][sel_method])$rt))
  }
}
rownames(TIMES) = methods_names
colnames(TIMES) = types
TIMES/60

types = c("de", "dp", "dm", "db", "dv", "nill")

TIMES = sapply(types, function(id){
  rowMeans(TIMES[,grep(id, colnames(TIMES))])
})

TIMES_avg = data.frame(TIMES, average = rowMeans(TIMES))/60

TIMES_avg = TIMES_avg[order(TIMES_avg$average, decreasing = TRUE),]

library(xtable)
xtable( TIMES_avg, digits = 1 )

gg_data = data.frame(method = methods_names, minutes = rowMeans(TIMES)/60)
# sort computational cost:
gg_data = gg_data[order(gg_data$minutes, decreasing = TRUE),]

gg_data$method = factor(gg_data$method, levels = gg_data$method)

# colour bars by method:
match = match( gg_data$method,  methods_names )
gg_data$cols = all_colours[match]


library(ggplot2)
ggplot() +
  geom_bar(data = gg_data, aes_string(x = "method", y = "minutes", fill = "method"), stat = "identity") +
  theme_bw() + 
  xlab("") +
  ylab("minutes") + 
  scale_y_sqrt( breaks = c(5, 30, 60, 150, 300, 500, 1000, 2000) ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(1.5)),
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),
        aspect.ratio = 1,
        legend.position = "none") +
  scale_fill_manual(values=gg_data$cols)


if(saving){
  ggsave(filename = paste0("Mean-time-(minutes).pdf"),
         plot = last_plot(),
         device = "pdf",
         path = "~/Desktop/distinct project/distinct Article/v1/images/muscat",
         width = 10,
         height = 10,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}
