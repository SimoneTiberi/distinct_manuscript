n_cells = 100
n_samples = 4

library(magrittr)
library(ggplot2)

set.seed(61217)
counts = matrix(nrow = 5, ncol = 2*n_samples*n_cells)
rownames(counts) = c("DE", "DP", "DM", "DB", "DV")

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
  
  # DV:
  counts[5, seq( (i-1) * n_cells + 1, i*n_cells) ] = rnorm(n_cells, mean = 10, sd = 1)
  counts[5, n_samples * n_cells + seq( (i-1) * n_cells + 1, i*n_cells) ] = rnorm(n_cells, mean = 10, sd = 2)
}

colData = list(sample_id = rep( seq_len(2*n_samples), each = n_cells ),
               cluster_id = rep( "", 2 * n_samples * n_cells ),
               group_id = rep( c("A", "B") , each = n_samples * n_cells ) )

# create a SingleCellExperiment:
library(SingleCellExperiment)
x = SingleCellExperiment(assays = list(expression = counts),
                         colData = colData)

# try distinct plots:
library(distinct)
# DE
DE_den = plot_densities(x, cluster = "",
                        gene = "DE",
                        name_assays_expression = "expression") +
  labs(title = "")

DE_cdf = plot_cdfs(x, cluster = "",
                   gene = "DE",
                   name_assays_expression = "expression") +
  labs(title = "")


# DP
DP_den = plot_densities(x, cluster = "",
                        gene = "DP",
                        name_assays_expression = "expression") +
  labs(title = "")

DP_cdf = plot_cdfs(x, cluster = "",
                   gene = "DP",
                   name_assays_expression = "expression") +
  labs(title = "")

# DM
DM_den = plot_densities(x, cluster = "",
                        gene = "DM",
                        name_assays_expression = "expression") +
  labs(title = "")

DM_cdf = plot_cdfs(x, cluster = "",
                   gene = "DM",
                   name_assays_expression = "expression") +
  labs(title = "")

# DB
DB_den = plot_densities(x, cluster = "",
                        gene = "DB",
                        name_assays_expression = "expression") +
  labs(title = "")

DB_cdf = plot_cdfs(x, cluster = "",
                   gene = "DB",
                   name_assays_expression = "expression") +
  labs(title = "")

# DV
DV_den = plot_densities(x, cluster = "",
                        gene = "DV",
                        name_assays_expression = "expression") +
  labs(title = "")

DV_cdf = plot_cdfs(x, cluster = "",
                   gene = "DV",
                   name_assays_expression = "expression") +
  labs(title = "")

library(ggpubr)
ggarrange(DE_den + xlab("") + labs(title = "DE") + theme(plot.title = element_text(hjust = 0.5, face = "bold")),
          DE_cdf+ xlab("") + ylab("ECDF") + labs(title = "DE") + theme(plot.title = element_text(hjust = 0.5, face = "bold")),
          DP_den + xlab("") + labs(title = "DP") + theme(plot.title = element_text(hjust = 0.5, face = "bold")),
          DP_cdf + xlab("") + ylab("ECDF") + labs(title = "DP") + theme(plot.title = element_text(hjust = 0.5, face = "bold")),
          DM_den + xlab("") + labs(title = "DM") + theme(plot.title = element_text(hjust = 0.5, face = "bold")),
          DM_cdf + xlab("") + ylab("ECDF") + labs(title = "DM") + theme(plot.title = element_text(hjust = 0.5, face = "bold")),
          DB_den + labs(title = "DB") + theme(plot.title = element_text(hjust = 0.5, face = "bold")),
          DB_cdf + labs(title = "DB") + ylab("ECDF") + theme(plot.title = element_text(hjust = 0.5, face = "bold")),
          DV_den + labs(title = "DV") + theme(plot.title = element_text(hjust = 0.5, face = "bold")),
          DV_cdf + labs(title = "DV") + ylab("ECDF") + theme(plot.title = element_text(hjust = 0.5, face = "bold")),
          ncol = 2, nrow = 5,
          legend = "right", common.legend = TRUE)

if(TRUE){
  ggsave(filename = "density_cdf_DE_DP_DM_DM.pdf",
         plot = last_plot(),
         device = "pdf",
         path = "~/Desktop/distinct project/distinct Article/v1/images/DE_profiles/",
         width = 8,
         height = 12,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}

