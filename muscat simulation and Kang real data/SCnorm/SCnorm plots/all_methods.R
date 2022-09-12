library(RColorBrewer);
all_colours = c(
  brewer.pal(7, "Reds")[7:2],   # distinct 5 methods
  brewer.pal(6, "Blues")[6:2],  # edgeR 4 methods
  brewer.pal(8, "Greens")[8:2], # limma-trend 1 method + # limma-voom 5 methods
  brewer.pal(4, "Greys")[4:2],  # MM 3 methods
  brewer.pal(6, "RdPu")[6:2],   # scDD-KS 4 methods
  brewer.pal(5, "BuPu")[5:2]    # scDD-perm 4 methods
)
# Purples Oranges Greys

all_methods = c(
  ",distinct.basics", # distinct.cpm25FALSE, or distinct.cpm25FALSE1000,
  ",distinct.cpm", # distinct.cpm25FALSE, or distinct.cpm25FALSE1000,
  ",distinct.linnorm",
  ",distinct.logcounts", # distinct.logcounts25FALSE, or distinct.logcounts25FALSE1000,
  ",distinct.vstresiduals", # distinct.vstresiduals25FALSE, or distinct.vstresiduals25FALSE1000,
  ",distinct.SCnorm", # distinct.vstresiduals25FALSE, or distinct.vstresiduals25FALSE1000,
  ",edgeR.mean.basics",
  ",edgeR.mean.linnorm",
  ",edgeR.sum.counts",
  ",edgeR.sum.scalecpm",
  ",edgeR.sum.SCnorm",
  ",limma-trend.mean.basics",
  ",limma-trend.mean.linnorm",
  ",limma-trend.mean.logcounts",
  ",limma-trend.mean.vstresiduals",
  ",limma-trend.sum.scalecpm",
  ",limma-voom.sum.counts",
  ",limma-trend.sum.SCnorm",
  ",MM-dream2",
  ",MM-nbinom",
  ",MM-vst",
  ",scDD.basics",
  ",scDD.cpm",
  ",scDD.linnorm",
  ",scDD.vstresiduals",
  ",scDD.SCnorm",
  ",permscdd.basics",
  ",permscdd.cpm",
  ",permscdd.linnorm",
  ",permscdd.vstresiduals"
)

# methods names:
methods_names = c(
  "distinct.basics", # distinct.cpm25FALSE, or distinct.cpm25FALSE1000,
  "distinct.cpm", # distinct.cpm25FALSE, or distinct.cpm25FALSE1000,
  "distinct.linnorm",
  "distinct.logcounts", # distinct.logcounts25FALSE, or distinct.logcounts25FALSE1000,
  "distinct.vstresiduals", # distinct.vstresiduals25FALSE, or distinct.vstresiduals25FALSE1000,
  "distinct.SCnorm", # distinct.vstresiduals25FALSE, or distinct.vstresiduals25FALSE1000,
  "edgeR.basics",
  "edgeR.linnorm",
  "edgeR.counts",
  "edgeR.cpm",
  "edgeR.SCnorm",
  "limma-trend.basics",
  "limma-trend.linnorm",
  "limma-trend.logcounts",
  "limma-trend.vstresiduals",
  "limma-trend.cpm",
  "limma-voom.counts",
  "limma-trend.SCnorm",
  "MM-dream2",
  "MM-nbinom",
  "MM-vstresiduals",
  "scDD-KS.basics",
  "scDD-KS.cpm",
  "scDD-KS.linnorm",
  "scDD-KS.vstresiduals",
  "scDD-KS.SCnorm",
  "scDD-perm.basics",
  "scDD-perm.cpm",
  "scDD-perm.linnorm",
  "scDD-perm.vstresiduals"
)

ordering = order(methods_names)

all_methods = all_methods[ordering]
methods_names = methods_names[ordering]
all_colours = all_colours[ordering]
