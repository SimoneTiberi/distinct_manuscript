MY_plot_cdf = function (x, name_assays_expression = "logcounts", name_cluster = "cluster_id", 
                        name_sample = "sample_id", name_group = "group_id", cluster, 
                        gene, group_level = FALSE, pad = TRUE, size = 0.75) 
{
  stopifnot((is(x, "SummarizedExperiment") | is(x, "SingleCellExperiment")), 
            is.character(name_assays_expression), length(name_assays_expression) == 
              1L, is.character(name_cluster), length(name_cluster) == 
              1L, is.character(name_sample), length(name_sample) == 
              1L, is.character(name_group), length(name_group) == 
              1L, is.character(cluster), length(cluster) == 1L, 
            is.character(gene), length(gene) == 1L, is.logical(group_level), 
            length(group_level) == 1L)
  sel = which(names(assays(x)) == name_assays_expression)
  if (length(sel) == 0) {
    message("'name_assays_expression' not found in names(assays(x))")
    return(NULL)
  }
  if (length(sel) > 1) {
    message("'name_assays_expression' found multiple times in names(assays(x))")
    return(NULL)
  }
  counts = assays(x)[[sel]]
  sel = which(names(colData(x)) == name_cluster)
  if (length(sel) == 0) {
    message("'name_cluster' not found in names(colData(x))")
    return(NULL)
  }
  if (length(sel) > 1) {
    message("'name_cluster' found multiple times in names(colData(x))")
    return(NULL)
  }
  cluster_ids = colData(x)[[sel]]
  if (!group_level) {
    sel = which(names(colData(x)) == name_sample)
    if (length(sel) == 0) {
      message("'name_sample' not found in names(colData(x))")
      return(NULL)
    }
    if (length(sel) > 1) {
      message("'name_sample' found multiple times in names(colData(x))")
      return(NULL)
    }
    sample_ids = factor(colData(x)[[sel]])
  }
  sel = which(names(colData(x)) == name_group)
  if (length(sel) == 0) {
    message("'name_group' not found in names(colData(x))")
    return(NULL)
  }
  if (length(sel) > 1) {
    message("'name_group' found multiple times in names(colData(x))")
    return(NULL)
  }
  group_ids = factor(colData(x)[[sel]])
  sel_gene = rownames(x) == gene
  sel_gene[is.na(sel_gene)] = FALSE
  if (all(!sel_gene)) {
    message("'gene' not found in `rownames(x)`")
    return(NULL)
  }
  sel_cluster = cluster_ids == cluster
  sel_cluster[is.na(sel_cluster)] = FALSE
  if (all(!sel_cluster)) {
    message("'cluster' not found in `colData(x)[[name_cluster]]`")
    return(NULL)
  }
  
  DF = data.frame(x = counts[sel_gene, sel_cluster], 
                  group = group_ids[sel_cluster])

  # Calculate density of data (stat_ecdf plots too many points!)
  dens = split(DF, DF$group) %>% 
    purrr::map_df(function(d) {
      dens = density(d$x, adjust=0.01, from=min(DF$x), # - 0.05*diff(range(DF$x)), 
                     to=max(DF$x))# + 0.05*diff(range(DF$x)))
      data.frame(x=dens$x, y=dens$y, cd=cumsum(dens$y)/sum(dens$y), group=d$group[1])
    })
  
  gg = ggplot(dens, aes(x=cd, group=group, colour=group, fill = group))
  
  # CDF plot:
  gg + geom_hline(yintercept = 0, linetype = "dashed", color = "grey") + 
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey") + 
    #stat_ecdf(pad = pad, size = size, geom = "line") +
    geom_line(data=dens, aes(x, cd, colour=group), size = size) +
    theme_bw() + theme(panel.grid = element_blank()) + 
    labs(title = paste(cluster, "-", gene),
         x = name_assays_expression, 
         y = "CDF")
}
