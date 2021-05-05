simData_DV_batch <- function(x, nc = 2e3, ns = 3, nk = 3, nb = 2,
                    probs = NULL, p_dd = diag(6)[1, ], paired = FALSE,
                    p_ep = 0.5, p_dp = 0.3, p_dm = 0.5,
                    p_type = 0, lfc = 2, rel_lfc = NULL, 
                    rel_be = NULL, rel_be_c = NULL, lfc_be = 0,
                    phylo_tree = NULL, phylo_pars = c(ifelse(is.null(phylo_tree), 0, 0.1), 3),
                    ng = nrow(x), force = FALSE,
                    batch_design, p_DV = 0.0) {
  
  cats = factor(c("ee", "ep", "de", "dp", "dm", "db", "dv"), 
                levels = c("ee", "ep", "de", "dp", "dm", "db", "dv") )
  names(cats) = cats
  # batch_design is a n_samples * n_groups matrix, specifying the batch each sample-group is to.
  
  # throughout this code...
  # k: cluster ID
  # s: sample ID
  # g: group ID
  # c: DD category
  # 0: reference
  
  # store all input arguments to be returned in final output
  args <- c(as.list(environment()))
  
  # check validity of input arguments
  muscat.batch:::.check_sce(x, req_group = FALSE)
  #args_tmp <- .check_args_simData(as.list(environment()))
  #nk <- args$nk <- args_tmp$nk
  p_dd = c(p_dd, p_DV)
  
  # reference IDs
  nk0 <- length(kids0 <- set_names(levels(x$cluster_id)))
  ns0 <- length(sids0 <- set_names(levels(x$sample_id)))
  nb0 <- length(bids0 <- set_names(levels(x$batch_id)))
  
  # simulation IDs
  nk <- length(kids <- set_names(paste0("cluster", seq_len(nk))))
  sids <- set_names(paste0("sample", seq_len(ns)))
  gids <- set_names(c("A", "B"))
  bids <- set_names(paste0("batch", seq_len(nb)))
  
  # sample reference clusters & samples
  ref_kids <- setNames(sample(kids0, nk, nk > nk0), kids)
  if (paired) { 
    # use same set of reference samples for both groups
    ref_sids <- sample(sids0, ns, ns > ns0)
    ref_sids <- replicate(length(gids), ref_sids)
  } else {
    # draw reference samples at random for each group
    ref_sids <- replicate(length(gids), 
                          sample(sids0, ns, ns > ns0))
  }
  dimnames(ref_sids) <- list(sids, gids)
  ref_bids <- setNames(sample(bids0, nb, nb > nb0), bids)
  
  # match batch and sample (if else to prevent sampling from 1 batch only)
  if (ns > nb) {
    s_by_b <- setNames(c(bids, sample(bids, ns-nb, TRUE)), sids)
  } else {
    s_by_b <- setNames(sample(bids, ns), sids)
  }
  
  
  if (is.null(rel_lfc)) 
    rel_lfc <- rep(1, nk)
  if (is.null(names(rel_lfc))) {
    names(rel_lfc) <- kids
  } else {
    stopifnot(names(rel_lfc) %in% kids0)
  }
  if (is.null(rel_be))
    rel_be <- rep(1, nb) 
  names(rel_be) <- bids
  if (is.null(rel_be_c))
    rel_be_c <- rep(1, nk)
  names(rel_be_c) <- kids
  
  # initialize count matrix
  gs <- paste0("gene", seq_len(ng))
  cs <- paste0("cell", seq_len(nc))
  y <- matrix(0, ng, nc, dimnames = list(gs, cs))
  
  # sample cell metadata
  cd <- muscat.batch:::.sample_cell_md(
    n = nc, probs = probs,
    ids = list(kids, sids, gids))
  rownames(cd) <- cs
  cs_idx <- muscat.batch:::.split_cells(cd, by = colnames(cd))
  n_cs <- modify_depth(cs_idx, -1, length)
  
  # split input cells by cluster-sample
  cs_by_ks <- muscat.batch:::.split_cells(x)
  
  # sample nb. of genes to simulate per category & gene indices
  n_dd <- table(sample(cats, ng, TRUE, p_dd))
  n_dd <- replicate(nk, n_dd)
  colnames(n_dd) <- kids
  gs_idx <- sample_gene_inds_DV(gs, n_dd, cats)
  
  # for ea. cluster, sample set of genes to simulate from
  gs_by_k <- setNames(sample(rownames(x), ng, ng > nrow(x)), gs)
  gs_by_k <- replicate(nk, gs_by_k)
  colnames(gs_by_k) <- kids
  
  # when 'phylo_tree' is specified, induce hierarchical cluster structure
  if (!is.null(phylo_tree)) {                                  
    res <- .impute_shared_type_genes(x, gs_by_k, gs_idx, phylo_tree, phylo_pars)
    gs_by_k <- res$gs_by_k
    class  <- res$class
    specs <- res$specs
    # otherwise, simply impute type-genes w/o phylogeny
  } else if (p_type != 0) {
    res <- .impute_type_genes(x, gs_by_k, gs_idx, p_type)
    stopifnot(!any(res$class == "shared"))
    gs_by_k <- res$gs_by_k
    class <- res$class
    specs <- res$specs
  } else {
    class <- rep("state", ng)
    specs <- rep(NA, ng)
    names(class) <- names(specs) <- gs
  }
  
  # split by cluster & categroy
  gs_by_k <- split(gs_by_k, col(gs_by_k))
  gs_by_k <- setNames(map(gs_by_k, set_names, gs), kids)
  gs_by_kc <- lapply(kids, function(k) 
    lapply(unfactor(cats), function(c) 
      gs_by_k[[k]][gs_idx[[c, k]]])) 
  
  #logFC of batch effect
  # ==========================================================================
  # Test if logFC are estimated from real data during prepSim and provided 
  # within rowData slot. If so pick a reference batch and get logFcs for other
  # batches for each cluster from rowData. If more cluster are initiated than 
  # present in rowData, additional ones are generated from the present ones. 
  # Estimated batch effects can only be used if enough reference batches are
  # present.
  # --------------------------------------------------------------------------
  lfc_all <- names(rowData(x))[grepl("logFC", names(rowData(x)))]
  batch_est <- nb0 >= nb & nb > 0 & length(lfc_all) > 0 & lfc_be == 0
  
  if (batch_est) {
    ref <- bids0[1]              
    lfc_ref <- lfc_all[grepl(paste0(ref,"-"), lfc_all)]
    lfc_b <- lapply(bids, function(b){
      br <- ref_bids[b]
      lfcb <- lfc_ref[grepl(br, lfc_ref)]
      lfcb <- muscat.batch:::.est_lfc_batch(lfcb, kids)
      lfcb_kc <- vapply(kids, function(k)
        lapply(cats, function(c){
          n_gskc <- length(gskc <- gs_by_kc[[k]][[c]])
          if (br %in% ref) {
            lfc_kc <- rep(0, n_gskc) %>% set_names(gskc)
          }else{
            lfc_kc <- lfcb[seq_len(n_gskc), paste0("logFC_", k)] %>% 
              set_names(gskc)
            lfc_kc <- lfc_kc * rel_be[[b]] * rel_be_c[[k]]
          }
          return(lfc_kc)
        }), vector("list", length(cats))) %>% set_rownames(cats)
      lfcb_kc
    })
  }
  
  # ==========================================================================
  # If logFcs are not provided, but a batch effect is defined in the 
  # lfc_be parameter this is used to sample a clusterspecific batch from a 
  # specified distribution.
  # --------------------------------------------------------------------------
  if (nb0 < nb & lfc_be == 0) {
    warning("No batch effect included into simulation.
                To few batches in reference. Please simulate less batches
                or specify 'lfc_be'.")
  }
  
  if (!batch_est) {
    lfc_b <- lapply(bids, function(b){
      lfcb <- muscat.batch:::.sim_lfc_be(lfc_be, gs)
      #celltype_specific logFC
      lfcb_kc <- vapply(kids, function(k)
        lapply(cats, function(c) {
          lfc_kc <- lfcb[names(gs_by_kc[[k]][[c]])]
          lfc_kc <- lfc_kc * rel_be[[b]] * rel_be_c[[k]]
          names(lfc_kc) <- gs_by_kc[[k]][[c]]
          lfc_kc
        }), vector("list", length(cats))) %>%
        set_rownames(cats)
      lfcb_kc
    })
  }
  
  # sample logFCs
  lfc <- vapply(kids, function(k) 
    lapply(unfactor(cats), function(c) { 
      n <- n_dd[c, k]
      if (c == "ee") return(rep(NA, n))
      signs <- sample(c(-1, 1), n, TRUE)
      lfcs <- rgamma(n, 4, 4/lfc) * signs
      names(lfcs) <- gs_by_kc[[k]][[c]]
      lfcs * rel_lfc[k]
    }), vector("list", length(cats)))
  
  # compute NB parameters
  m <- lapply(sids0, function(s) {
    b <- paste0("beta.", s)
    b <- exp(rowData(x)[[b]])
    m <- outer(b, exp(x$offset), "*")
    dimnames(m) <- dimnames(x); m
  })
  d <- rowData(x)$dispersion 
  names(d) <- rownames(x)
  
  # initialize list of depth two to store 
  # simulation means in each cluster & group
  sim_mean <- lapply(kids, function(k) 
    lapply(gids, function(g) 
      setNames(numeric(ng), gs)))
  
  # run simulation -----------------------------------------------------------
  rownames(batch_design) = sids
  
  for (k in kids) {
    for (s in sids) {
      # get reference samples, clusters & cells
      s0 <- ref_sids[s, ]
      k0 <- ref_kids[k]
      cs0 <- cs_by_ks[[k0]][s0]
      
      # get output cell indices
      ci <- cs_idx[[k]][[s]]
      
      for (c in cats[n_dd[, k] != 0]) {
        # sample cells to simulate from
        cs_g1 <- sample(cs0[[1]], n_cs[[k]][[s]][[1]], TRUE)
        cs_g2 <- sample(cs0[[2]], n_cs[[k]][[s]][[2]], TRUE)
        
        # get reference genes & output gene indices
        gs0 <- gs_by_kc[[k]][[c]] 
        gi <- gs_idx[[c, k]]
        
        #Add  batch specific logfoldchange
        b <- batch_design[s, 1] # s_by_b[s]
        lfc_be_s_1 <- lfc_b[[b]][[c, k]][gs0] 
        lfc_be_s_1[is.na(lfc_be_s_1)] <- 0
        
        b <- batch_design[s, 2] # s_by_b[s]
        lfc_be_s_2 <- lfc_b[[b]][[c, k]][gs0] 
        lfc_be_s_2[is.na(lfc_be_s_2)] <- 0
        
        # get NB parameters
        m_g1 <- m[[s0[[1]]]][gs0, cs_g1, drop = FALSE] * 2^lfc_be_s_1
        m_g2 <- m[[s0[[2]]]][gs0, cs_g2, drop = FALSE] * 2^lfc_be_s_2
        d_kc <- d[gs0]
        lfc_kc <- lfc[[c, k]]
        
        re <- sim_DV_batch(c, cs_g1, cs_g2, m_g1, m_g2, d_kc, lfc_kc, p_ep, p_dp, p_dm)
        y[gi, unlist(ci)] <- re$cs
        
        for (g in gids) sim_mean[[k]][[g]][gi] <- ifelse(
          is.null(re$ms[[g]]), NA, list(re$ms[[g]]))[[1]]
      }
    }
  }
  # construct gene metadata table storing ------------------------------------
  # adjust lfc_be to summary table
  red_lfc_be <- lapply(bids, function(b){
    lfc_red <- lfc_b[[b]] %>% unlist()
  }) %>% bind_cols() %>% 
    set_colnames(paste0("lfc_be_", colnames(.))) %>%
    mutate("gene" = rep(gs, nk), "cluster_id" = rep(kids, each = length(gs)))
  
  # gene | cluster_id | category | logFC, gene, disp, mean used for sim.
  gi <- data.frame(
    gene = unlist(gs_idx),
    cluster_id = rep.int(rep(kids, each = length(cats)), c(n_dd)),
    category = rep.int(rep(cats, nk), c(n_dd)),
    logFC = unlist(lfc),
    sim_gene = unlist(gs_by_kc),
    sim_disp = d[unlist(gs_by_kc)]) %>% 
    mutate_at("gene", as.character)
  
  # add lfc_be
  gi <- full_join(gi, red_lfc_be, by = c("gene", "cluster_id"))
  # add true simulation means
  sim_mean <- sim_mean %>%
    map(bind_cols) %>% 
    bind_rows(.id = "cluster_id") %>% 
    mutate_at("cluster_id", factor) %>% 
    mutate(gene = rep(gs, nk))
  gi <- full_join(gi, sim_mean, by = c("gene", "cluster_id")) %>% 
    dplyr::rename("sim_mean.A" = "A", "sim_mean.B" = "B")
  # reorder
  o <- order(as.numeric(gsub("[a-z]", "", gi$gene)))
  gi <- gi[o, ]; rownames(gi) <- NULL
  
  # construct SCE ------------------------------------------------------------
  # cell metadata including group, sample, batch, cluster IDs
  cd <- cd %>% mutate("batch_id" = recode(as.factor(sample_id), !!!s_by_b))
  cd$group_id <- droplevels(cd$group_id)
  cd$sample_id <- factor(paste(cd$sample_id, cd$group_id, sep = "."))
  m <- match(levels(cd$sample_id), cd$sample_id)
  gids <- cd$group_id[m]
  o <- order(gids)
  sids <- levels(cd$sample_id)[o]
  ei <- data.frame(sample_id = sids, group_id = gids[o])
  cd <- cd %>% mutate_at("sample_id", factor, levels = sids)
  # gene metadata storing gene classes & specificities
  rd <- DataFrame(class = factor(class, 
                                 levels = c("state", "shared", "type")))
  rd$specs <- as.list(specs)
  # simulation metadata including used reference samples/cluster, 
  # list of input arguments, and simulated genes' metadata
  md <- list(
    experiment_info = ei,
    n_cells = table(cd$sample_id),
    gene_info = gi,
    ref_sids = ref_sids,
    ref_kids = ref_kids, 
    args = args)
  # return SCE 
  res = SingleCellExperiment(
    assays = list(counts = as.matrix(y)),
    colData = cd, rowData = rd, metadata = md)
  res
}

sample_gene_inds_DV = function (gs, ns, cats){
  cluster_ids <- colnames(ns)
  vapply(cluster_ids, function(k) split(sample(gs), rep.int(cats, ns[, k])), vector("list", length(cats)))
}

sim_DV_batch <- function(
  cat = c("ee", "ep", "de", "dp", "dm", "db", "dv"), 
  cs_g1, cs_g2, m_g1, m_g2, d, lfc, ep, dp, dm) {
  
  cat <- match.arg(cat)
  ng1 <- length(cs_g1)
  ng2 <- length(cs_g2)
  
  re <- switch(cat,
               "ee" = {
                 list(
                   muscat.batch:::.nb(cs_g1, d, m_g1),
                   muscat.batch:::.nb(cs_g2, d, m_g2))
               },
               "ep" = {
                 g1_hi <- sample(ng1, round(ng1 * ep))
                 g2_hi <- sample(ng2, round(ng2 * ep))
                 list(
                   muscat.batch:::.nb(cs_g1[-g1_hi], d, m_g1),
                   muscat.batch:::.nb(cs_g1[ g1_hi], d, m_g1, lfc), # 50% g1 hi
                   muscat.batch:::.nb(cs_g2[-g2_hi], d, m_g2),
                   muscat.batch:::.nb(cs_g2[ g2_hi], d, m_g2, lfc)) # 50% g2 hi
               },
               "de" = {
                 list(
                   muscat.batch:::.nb(cs_g1, d, m_g1, -lfc), # lfc < 0 => all g1 hi
                   muscat.batch:::.nb(cs_g2, d, m_g2,  lfc)) # lfc > 0 => all g2 hi
               },
               "dp" = {
                 props <- sample(c(dp, 1 - dp), 2)
                 g1_hi <- sample(ng1, round(ng1 * props[1]))
                 g2_hi <- sample(ng2, round(ng2 * props[2]))
                 list(                           
                   muscat.batch:::.nb(cs_g1[-g1_hi], d, m_g1), 
                   muscat.batch:::.nb(cs_g1[ g1_hi], d, m_g1,  lfc), # lfc > 0 => dp/(1-dp)% up
                   muscat.batch:::.nb(cs_g2[-g2_hi], d, m_g2), 
                   muscat.batch:::.nb(cs_g2[ g2_hi], d, m_g2, -lfc)) # lfc < 0 => (1-dp)/dp% up
               },
               "dm" = {
                 g1_hi <- sample(ng1, round(ng1 * dm))
                 g2_hi <- sample(ng2, round(ng2 * dm))
                 list(
                   muscat.batch:::.nb(cs_g1[-g1_hi], d, m_g1),
                   muscat.batch:::.nb(cs_g1[ g1_hi], d, m_g1, -lfc), # lfc < 0 => 50% g1 hi
                   muscat.batch:::.nb(cs_g2[-g2_hi], d, m_g2),
                   muscat.batch:::.nb(cs_g2[ g2_hi], d, m_g2,  lfc)) # lfc > 0 => 50% g2 hi
               }, 
               "dv" = {
                 d_1 = d * 2^lfc
                 list(
                   muscat.batch:::.nb(cs_g1, d_1, m_g1),
                   muscat.batch:::.nb(cs_g2, d, m_g2))
               },
               "db" = {
                 if (sample(c(TRUE, FALSE), 1)) {
                   # all g1 mi, 50% g2 hi
                   g2_hi <- sample(ng2, round(ng2 * 0.5))
                   list(
                     muscat.batch:::.nb(cs_g1, d, m_g1, abs(lfc), 0.5),
                     muscat.batch:::.nb(cs_g2[-g2_hi], d, m_g2, -lfc), 
                     muscat.batch:::.nb(cs_g2[ g2_hi], d, m_g2,  lfc)) 
                 } else {
                   # all g2 mi, 50% g1 hi
                   g1_hi <- sample(ng1, round(ng1 * 0.5))
                   list(
                     muscat.batch:::.nb(cs_g2, d, m_g2, abs(lfc), 0.5), 
                     muscat.batch:::.nb(cs_g1[-g1_hi], d, m_g1, -lfc),  
                     muscat.batch:::.nb(cs_g1[ g1_hi], d, m_g1,  lfc))  
                 }
               }
  )
  cs <- map(re, "counts")
  cs <- do.call("cbind", cs)
  ms <- map(re, "means")
  rmv <- vapply(ms, is.null, logical(1))
  ms <- ms[!rmv] %>% 
    map_depth(2, mean) %>% 
    map_depth(1, unlist) %>% 
    bind_cols %>% 
    as.matrix
  ms <- switch(cat, 
               ee = ms,
               de = ms,
               dv = ms,
               db = if (ng2 == 0) {
                 as.matrix(
                   ms[, 1])
               } else {
                 cbind(
                   ms[, 1],
                   rowMeans(ms[, c(2, 3)]))
               }, if (ng2 == 0) {
                 as.matrix(
                   rowMeans(ms[, c(1, 2)]))
               } else {
                 cbind(
                   rowMeans(ms[, c(1, 2)]),
                   rowMeans(ms[, c(3, 4)]))
               })
  ms <- split(ms, col(ms))
  names(ms) <- c("A", "B")[c(ng1, ng2) != 0]
  list(cs = cs, ms = ms)
}

