rm(list = ls())
setwd("~/Desktop/distinct project/SIMULATED data/FULL muscat pipeline")

library(ggplot2)
saving = FALSE
pairing = "paired" 

# NEW RESULTS:
res_names = list.files(paste0("results/results_7Samples"), full.names =TRUE)
file.exists(res_names)

source("~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/all_methods.R")
methods = all_methods
methods[1:5] = paste0(methods[1:5], "3")

methods_names = methods_names 

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

ggsave(filename = paste0("Mean-time-(minutes).pdf"),
       plot = last_plot(),
       device = "pdf",
       path = "~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/plots muscat simulation - NEW normalization/Plots/",
       width = 10,
       height = 10,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

