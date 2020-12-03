n_cells = 10^5
n_samples = 1

library(magrittr)

set.seed(61217)
counts = matrix(nrow = 4, ncol = 2*n_samples*n_cells)
rownames(counts) = c("DE", "DP", "DM", "DB")

for(i in seq_len(n_samples)){
  # DE:
  counts[1, seq( (i-1) * n_cells + 1, i*n_cells) ] = rnorm(n_cells, mean = 9)
  counts[1, n_samples * n_cells + seq( (i-1) * n_cells + 1, i*n_cells) ] = rnorm(n_cells, mean = 11)
  
  # DP:
  counts[2, seq( (i-1) * n_cells + 1, i*n_cells) ] = c(rnorm(n_cells * 0.65, mean = 8), rnorm(n_cells * 0.35, mean = 12))
  counts[2, n_samples * n_cells + seq( (i-1) * n_cells + 1, i*n_cells)  ] = c(rnorm(n_cells * 0.35, mean = 8), rnorm(n_cells * 0.65, mean = 12))
  
  # DM:
  counts[3, seq( (i-1) * n_cells + 1, i*n_cells) ] =  rnorm(n_cells, mean = 8)
  counts[3, n_samples * n_cells + seq( (i-1) * n_cells + 1, i*n_cells)  ] = c(rnorm(n_cells * 0.35, mean = 8), rnorm(n_cells * 0.65, mean = 12))
  
  # DB:
  counts[4, seq( (i-1) * n_cells + 1, i*n_cells) ] =rnorm(n_cells, mean = 10)
  counts[4, n_samples * n_cells + seq( (i-1) * n_cells + 1, i*n_cells) ] = c(rnorm(n_cells/2, mean = 8), rnorm(n_cells/2, mean = 12))
}

colData = list(sample_id = rep( seq_len(2*n_samples), each = n_cells ),
               cluster_id = rep( "", 2 * n_samples * n_cells ),
               group_id = rep( c("A", "B") , each = n_samples * n_cells ) )

# create a SingleCellExperiment:
library(SingleCellExperiment)
x = SingleCellExperiment(assays = list(expression = counts),
                         colData = colData)

# try distinct plots:
library(ggplot2); library(distinct)
source("MY_plot_densities.R")
# DE
DE_den = MY_plot_densities(x, cluster = "",
                           gene = "DE",
                           name_assays_expression = "expression", adjust = 2, size = 1.5) +
  labs(title = "") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

source("My_plot_cdf.R")
DE_cdf = MY_plot_cdf(x, cluster = "",
                     gene = "DE",
                     name_assays_expression = "expression", size = 1.5) +
  labs(title = "") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) + theme(legend.position="bottom")


# DP
DP_den = MY_plot_densities(x, cluster = "",
                           gene = "DP",
                           name_assays_expression = "expression", adjust = 2, size = 1.5) +
  labs(title = "") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

DP_cdf = MY_plot_cdf(x, cluster = "",
                     gene = "DP",
                     name_assays_expression = "expression", size = 1.5) +
  labs(title = "") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

# DM
DM_den = MY_plot_densities(x, cluster = "",
                           gene = "DM",
                           name_assays_expression = "expression", adjust = 2, size = 1.5) +
  labs(title = "") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

DM_cdf = MY_plot_cdf(x, cluster = "",
                     gene = "DM",
                     name_assays_expression = "expression", size = 1.5) +
  labs(title = "") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

# DB
DB_den = MY_plot_densities(x, cluster = "",
                           gene = "DB",
                           name_assays_expression = "expression", adjust = 2, size = 1.5) +
  labs(title = "") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

DB_cdf = MY_plot_cdf(x, cluster = "",
                     gene = "DB",
                     name_assays_expression = "expression", size = 1.5) +
  labs(title = "") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())


library(ggpubr)
ggarrange(DE_cdf,
          DP_cdf, 
          DM_cdf, 
          DB_cdf,
          labels = c("DE",
                     "DP",
                     "DM",
                     "DB"),
          ncol = 2, nrow = 2,
          legend = "bottom", common.legend = TRUE)

library(ggpubr)
ggarrange(DE_den,
          DP_den, 
          DM_den, 
          DB_den, 
          labels = c("DE",
                     "DP",
                     "DM",
                     "DB"),
          ncol = 2, nrow = 2,
          legend = "bottom", common.legend = TRUE)

if(TRUE){
  ggsave(filename = "DE_profiles.pdf",
         plot = last_plot(),
         device = "pdf",
         path = "~/Desktop/distinct project/distinct Article/v1/images/DE_profiles",
         width = 6,
         height = 6,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}

library(ggpubr)
ggarrange(DE_den,
          DP_den, 
          DM_den, 
          DB_den, 
          DE_cdf,
          DP_cdf, 
          DM_cdf, 
          DB_cdf,
          labels = c("DE",
                     "DP",
                     "DM",
                     "DB",
                     "",
                     "",
                     "",
                     ""),
          legend.grob = get_legend(DE_cdf),
          ncol = 4, nrow = 2,
          legend = "bottom", common.legend = TRUE)




library(ggpubr)
ggarrange(DE_den, DE_cdf,
          DP_den, DP_cdf, 
          DM_den, DM_cdf, 
          DB_den, DB_cdf,
          labels = c("DE","",
                     "DP","",
                     "DM","",
                     "DB",""),
          legend.grob = get_legend(DE_cdf),
          ncol = 2, nrow = 4,
          legend = "bottom", common.legend = TRUE)

if(TRUE){
  ggsave(filename = "DE_density_cdf.pdf",
         plot = last_plot(),
         device = "pdf",
         path = "~/Desktop/distinct project/distinct Article/v1/images/DE_profiles",
         width = 6,
         height = 8,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}


ggarrange(DE_den + xlab("") + labs(title = "DE") + theme(plot.title = element_text(hjust = 0.5, face = "bold")),
          DP_den + xlab("") + ylab("") + labs(title = "DP") + theme(plot.title = element_text(hjust = 0.5, face = "bold")),
          DM_den + xlab("") + ylab("") + labs(title = "DM") + theme(plot.title = element_text(hjust = 0.5, face = "bold")),
          DB_den + xlab("") + ylab("") + labs(title = "DB") + theme(plot.title = element_text(hjust = 0.5, face = "bold")),
          DE_cdf,
          DP_cdf + ylab(""),
          DM_cdf + ylab(""),
          DB_cdf + ylab(""),
          ncol = 4, nrow = 2,
          legend = "bottom", common.legend = TRUE)

if(TRUE){
  ggsave(filename = "DE_density_cdf_horizontal.pdf",
         plot = last_plot(),
         device = "pdf",
         path = "~/Desktop/distinct project/distinct Article/v1/images/DE_profiles",
         width = 7,
         height = 4,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}


