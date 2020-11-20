rm(list = ls())
setwd("~/Desktop/distinct project/EXPERIMENTAL data/")

load("T cells from 12 Colorectal cancer patients - Scientific Data 2019/p_T_cells.RData")
load("Kang (muscData)/p_Kang.RData")

ggpubr::ggarrange(p_T_cells, p_Kang,
                  labels = c("T-cells", "Kang"),
                  ncol = 1, nrow = 2,
                  legend = "bottom", common.legend = TRUE)

# Figure 2b
ggsave(filename = paste0("null-JOINT.pdf"),
       plot = last_plot(),
       device = "pdf",
       path = "~/Desktop/distinct project/distinct Article/v1/images/NULL",
       width = 7.5,
       height = 7,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
