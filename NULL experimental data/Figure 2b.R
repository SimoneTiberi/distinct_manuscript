rm(list = ls())
setwd("~/Desktop/distinct project/EXPERIMENTAL data/")

load("T cells from 12 Colorectal cancer patients - Scientific Data 2019/p_T_cells.RData")
load("Kang (muscData)/p_Kang.RData")


# p_Kang = p_Kang + theme(legend.key.size = unit(0.5, 'cm'))

legend <- ggpubr::get_legend(p_Kang)

AA = egg::ggarrange( plots = 
                       list(p_T_cells + xlab("") + labs(title = "T-cells") + theme(legend.position = "none", aspect.ratio = 2),
                            p_Kang + labs(title = "Kang") + theme(legend.position = "none", aspect.ratio = 2)),
                     bottom = legend,
                     ncol = 1, nrow = 2)
#AA

#ggpubr::ggarrange(p_T_cells, p_Kang,
#                  labels = c("T-cells", "Kang"),
#                  ncol = 1, nrow = 2,
#                  legend = "bottom", common.legend = TRUE)

# try WITHOUT LEGEND TO MAKE MORE ROOM FOR THE IMAGES

library(ggplot2)
# Figure 2b
ggsave(filename = paste0("null-JOINT.pdf"),
       plot = AA,
       device = "pdf",
       path = "~/Desktop/distinct project/distinct Article/scripts/NEW - revision analyses/NULL experimental data/Figures",
       width = 5.75,
       height = 8.5,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
