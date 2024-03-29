config <- yaml::read_yaml("config.yaml")

de10 <- c(0.9, 0, 0.1, 0, 0, 0)

sim_pars <- list(
  # pure simulations
  de200_significant_20 = list(nr = 5, p_dd = c(0.8, 0, 0.2, 0, 0, 0), seed = 2, lfc = 1, nc = 18 * 200),
  dp200_significant_20 = list(nr = 5, p_dd = c(0.8, 0, 0, 0.2, 0, 0), seed = 3, lfc = 1.5, nc = 18 * 200),
  dm200_significant_20 = list(nr = 5, p_dd = c(0.8, 0, 0, 0, 0.2, 0), seed = 4, lfc = 1.5, nc = 18 * 200),
  db200_significant_20 = list(nr = 5, p_dd = c(0.8, 0, 0, 0, 0, 0.2), seed = 5, lfc = 3, nc = 18 * 200),
  dv200_significant_20 = list(nr = 5, p_dd = c(0.8, 0, 0, 0, 0, 0), seed = 6, lfc = 3, nc = 18 * 200, p_dv = 0.2),
  # pure simulations
  de200_significant_5 = list(nr = 5, p_dd = c(0.95, 0, 0.05, 0, 0, 0), seed = 2, lfc = 1, nc = 18 * 200),
  dp200_significant_5 = list(nr = 5, p_dd = c(0.95, 0, 0, 0.05, 0, 0), seed = 3, lfc = 1.5, nc = 18 * 200),
  dm200_significant_5 = list(nr = 5, p_dd = c(0.95, 0, 0, 0, 0.05, 0), seed = 4, lfc = 1.5, nc = 18 * 200),
  db200_significant_5 = list(nr = 5, p_dd = c(0.95, 0, 0, 0, 0, 0.05), seed = 5, lfc = 3, nc = 18 * 200),
  dv200_significant_5 = list(nr = 5, p_dd = c(0.95, 0, 0, 0, 0, 0), seed = 6, lfc = 3, nc = 18 * 200, p_dv = 0.05),
  # pure simulations
  de200_significant_1 = list(nr = 5, p_dd = c(0.99, 0, 0.01, 0, 0, 0), seed = 2, lfc = 1, nc = 18 * 200),
  dp200_significant_1 = list(nr = 5, p_dd = c(0.99, 0, 0, 0.01, 0, 0), seed = 3, lfc = 1.5, nc = 18 * 200),
  dm200_significant_1 = list(nr = 5, p_dd = c(0.99, 0, 0, 0, 0.01, 0), seed = 4, lfc = 1.5, nc = 18 * 200),
  db200_significant_1 = list(nr = 5, p_dd = c(0.99, 0, 0, 0, 0, 0.01), seed = 5, lfc = 3, nc = 18 * 200),
  dv200_significant_1 = list(nr = 5, p_dd = c(0.99, 0, 0, 0, 0, 0), seed = 6, lfc = 3, nc = 18 * 200, p_dv = 0.01)
)

def_pars <- list(nr = 1, nk = 3, ns = 3, 
                 lfc = 2, ng = 4e3, nc = function(nk, ns) 4*nk*ns*200, 
		p_dd = de10, probs = NULL, seed = 1, p_dv = 0)

sim_pars <- lapply(sim_pars, function(u) {
    v <- !names(def_pars) %in% names(u)
    u[names(def_pars)[v]] <- def_pars[v]
    if (is.function(u$nc))
        u$nc <- u$nc(u$nk, u$ns)
    return(u)
})

sim_ids <- names(sim_pars)
write(jsonlite::toJSON(sim_ids), config$sids)

# write parameters to .json (only if something changed!)

fns <- paste0(sim_ids, ".json")
fns <- paste0(config$sim_pars, fns)
names(fns) <- sim_ids

for (id in sim_ids) {
    new <- sim_pars[[id]]
    if (file.exists(fns[id])) {
        old <- jsonlite::fromJSON(fns[id])
        if (!isTRUE(all.equal(old, new, tolerance = 1e-3)))
            write(jsonlite::toJSON(new, null = "null"), fns[id])
    } else {
        write(jsonlite::toJSON(new, null = "null"), fns[id])
    }
}
