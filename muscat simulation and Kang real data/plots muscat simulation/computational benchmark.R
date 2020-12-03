rm(list = ls())
setwd("~/Desktop/distinct project/SIMULATED data/FULL muscat pipeline")

pairing = "paired" 

types = c("de", "dp", "dm", "db", "nill,1", "nill,2", "nill,3")

# NEW RESULTS:
res_names = list.files(paste0("results/results_7Samples"), full.names =TRUE)
file.exists(res_names)

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
  "MM-vstresiduals")


# filter results for n_cells:
n_cells = 200
res_names = res_names[ grep(n_cells, res_names) ]

n_rep = 5
types = paste0( c("de", "dp", "dm", "db", "nill") ,n_cells, ",", rep(1:n_rep, each= 5))

TIMES = matrix(NA, nrow = length(methods), ncol = length(types))

apply(TIMES/60, 1, mean)

for(i in seq_along(types)){
  sel_type = grep(types[i], res_names)
  
  for(j in seq_along(methods)){
    sel_method = grep(methods[j], res_names[sel_type])
    TIMES[j,i] = sum(as.numeric(readRDS(res_names[sel_type][sel_method])$rt))
  }
}
rownames(TIMES) = methods_names
colnames(TIMES) = types
TIMES/60

types = c("de", "dp", "dm", "db", "nill")

TIMES = sapply(types, function(id){
  rowMeans(TIMES[,grep(id, colnames(TIMES))])
})

TIMES_avg = data.frame(TIMES, average = rowMeans(TIMES))/60

TIMES_avg = TIMES_avg[order(TIMES_avg$average, decreasing = TRUE),]

library(xtable)
xtable( TIMES_avg, digits = 1 )
# Supplementary Table 1

gg_data = data.frame(method = methods_names, minutes = rowMeans(TIMES)/60)
# sort computational cost:
gg_data = gg_data[order(gg_data$minutes, decreasing = TRUE),]

gg_data$method = factor(gg_data$method, levels = gg_data$method)

library(ggplot2)
ggplot() +
  geom_bar(data = gg_data, aes_string(x = "method", y = "minutes"), stat = "identity", fill = "red") +
  theme_bw() + 
  xlab("") +
  ylab("minutes") + 
  scale_y_sqrt( breaks = c(0, 1, 5, 30, 60, 150, 300, 600, 800) ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(1.5)),
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),
        aspect.ratio = 0.6 )

ggsave(filename = paste0("Mean-time-(minutes).pdf"),
       plot = last_plot(),
       device = "pdf",
       path = "~/Desktop/distinct project/distinct Article/v1/images/muscat",
       width = 10,
       height = 8,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
