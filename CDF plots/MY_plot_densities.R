MY_plot_densities = function(x, 
                             name_assays_expression = "logcounts",
                             name_cluster = "cluster_id",
                             name_sample = "sample_id",
                             name_group = "group_id",
                             cluster,
                             gene,
                             group_level = TRUE,
                             adjust = 1,
                             size = 0.75,
                             alpha = 0.25){
  
  stopifnot(
    ( is(x, "SummarizedExperiment") | is(x, "SingleCellExperiment") ),
    is.character(name_assays_expression), length(name_assays_expression) == 1L,
    is.character(name_cluster), length(name_cluster) == 1L,
    is.character(name_sample), length(name_sample) == 1L,
    is.character(name_group), length(name_group) == 1L,
    is.character(cluster), length(cluster) == 1L,
    is.character(gene), length(gene) == 1L,
    is.logical(group_level), length(group_level) == 1L
  )
  
  # count matrix:
  sel = which(names(assays(x)) == name_assays_expression)
  if( length(sel) == 0 ){
    message("'name_assays_expression' not found in names(assays(x))")
    return(NULL)
  }
  if( length(sel) > 1 ){
    message("'name_assays_expression' found multiple times in names(assays(x))")
    return(NULL)
  }
  counts = assays(x)[[sel]]
  
  # cluster ids:
  sel = which(names(colData(x)) == name_cluster)
  if( length(sel) == 0 ){
    message("'name_cluster' not found in names(colData(x))")
    return(NULL)
  }
  if( length(sel) > 1 ){
    message("'name_cluster' found multiple times in names(colData(x))")
    return(NULL)
  }
  cluster_ids = colData(x)[[sel]]
  
  # sample ids:
  if(!group_level){
    sel = which(names(colData(x)) == name_sample)
    if( length(sel) == 0 ){
      message("'name_sample' not found in names(colData(x))")
      return(NULL)
    }
    if( length(sel) > 1 ){
      message("'name_sample' found multiple times in names(colData(x))")
      return(NULL)
    }
    sample_ids = factor(colData(x)[[sel]])
  }
  
  # group ids (1 per cell)
  sel = which(names(colData(x)) == name_group)
  if( length(sel) == 0 ){
    message("'name_group' not found in names(colData(x))")
    return(NULL)
  }
  if( length(sel) > 1 ){
    message("'name_group' found multiple times in names(colData(x))")
    return(NULL)
  }
  group_ids = factor(colData(x)[[sel]])
  
  # sel gene:
  sel_gene = rownames(x) == gene
  sel_gene[is.na(sel_gene)] = FALSE
  
  if( all(!sel_gene) ){ # if all sel_gene FALSE:
    message("'gene' not found in `rownames(x)`")
    return(NULL)
  }
  
  # sel cluster:
  sel_cluster = cluster_ids == cluster
  sel_cluster[is.na(sel_cluster)] = FALSE
  
  if( all(!sel_cluster) ){ # if all sel_gene FALSE:
    message("'cluster' not found in `colData(x)[[name_cluster]]`")
    return(NULL)
  }
  
  if(group_level){
    DF = data.frame(x = counts[sel_gene, sel_cluster], 
                    group = group_ids[sel_cluster])
    gg = ggplot(DF, aes(x=x, group=group, colour=group, fill = group))
  }else{
    DF = data.frame(x = counts[sel_gene, sel_cluster], 
                    group = group_ids[sel_cluster],
                    sample = sample_ids[sel_cluster])
    gg = ggplot(DF, aes(x=x, group=sample, colour=group, fill = group))
  }
  
  # Density plot:
  gg +
    geom_density(aes(x=x, colour=group, fill = group), 
                 adjust = adjust,
                 size = size,
                 # geom="line",
                 position="identity",
                 alpha = alpha) +
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    labs(title = paste(cluster, "-", gene),
         x = name_assays_expression)
}
