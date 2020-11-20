suppressMessages({
  library(muscat)
  library(muscat.batch)
  library(magrittr)
  })
source(file.path("scripts", "utils.R"))

# load reference SCE
ref <- readRDS(args$input_sce)

# prep. SCE for simulation w/ 'muscat::simData'
sce <- muscat.batch:::prepSim(ref, verbose = FALSE,
    # keep genes w/ count > 1 in >= 10 cells
    min_count = 1, min_cells = 10,   
    # keep cells w/ >= 100 detected genes & cluster w/ > 100 cells
    min_genes = 100, min_size = 100,
    batch = TRUE) 

# write SCE to .rds
saveRDS(sce, args$output_sce)
