# this determines which types of methods to include
names(ids) <- ids <- c("pb", "distinct", "scdd") 

# UPDATE PB as:
pb <- dplyr::bind_rows(
  expand.grid(
    stringsAsFactors = FALSE, scale = FALSE,
    assay = c("SCnorm"),
    fun = "mean", method = c("edgeR", "limma-trend")
  )
)
pb$treat[is.na(pb$treat)] <- FALSE
pb$id <- with(pb, sprintf("%s%s.%s.%s%s",
                          method, ifelse(treat, "-treat", ""),
                          fun, ifelse(scale, "scale", ""), assay))

# cdf (distinct) -----------------------------------------------------------------
distinct <- expand.grid(
  stringsAsFactors = FALSE,
  cores = c(3),
  assay = c("SCnorm"))
distinct$id <- with(distinct, paste0("distinct.", assay, cores))

# scDD -------------------------------------------------------------------------
scdd <- data.frame(
    stringsAsFactors = FALSE,
    assay = c("SCnorm"))
scdd$id <- with(scdd, sprintf("scDD.%s", assay))

# write method IDs to .csv -----------------------------------------------------
for (id in ids) {
    ps <- get(id)
    assign(id, split(ps, ps$id))
}

pars <- sapply(ids, get)
typs <- rep.int(ids, sapply(pars, length))
pars <- purrr::flatten(pars)
mids <- names(pars)

mids_df <- data.frame(
    row.names = NULL,
    stringsAsFactors = FALSE,
    id = mids, type = typs)

write.csv(mids_df, config$mids)

# write method parameters to .json ---------------------------------------------
fns <- paste0(mids, ".json")
fns <- paste0(config$meth_pars, fns)
names(fns) <- mids

for (id in mids) {
    new <- as.list(pars[[id]])
    if (file.exists(fns[id])) {
        old <- jsonlite::fromJSON(fns[id])
        if (!identical(old, new))
            write(jsonlite::toJSON(new), fns[id])
    } else {
        write(jsonlite::toJSON(new), fns[id])
    }
}

