This snakemakes was used to perform the `muscat` simulation study with batch effects.

This snakemake was obtained by editing the original snakemake from Crowell et al. (2020), available at: https://github.com/HelenaLC/muscat-comparison

This `snakemake` workflow is organized into

- a `config.yaml` file specify key parameters and directories
- a `scripts` folder housing all utilized scripts (see below)
- a `data` folder containing raw (reference) and simulated data
- a `meta` folder for simulation, runmode, and method parameters
- a `results` folder where all results are generated (as `.rds` files)

The table below summarizes the different R scripts in `scripts`:

script      | description 
:-----------|:-----------------------------------------------
`prep_kang`    | generates a references SCE for simulation by<br>i) keeping samples from one condition only; and,<br>ii) unifying relevant cell metadata names to `"cluster/sample/group_id"`
`prep_sim` | prepares a reference SCE for simulation by<br>i) retaining subpopulation-sample combinations with at least 100 cells; and,<br>ii) estimating cell / gene parameters (offsets / coefficients and dispersions)
`sim_pars`  | for ea. simulation ID, generates a `.json` file in `meta/sim_pars`<br>that specifies simulation parameters (e.g., prob. of DS, nb. of simulation replicates)
`run_pars`  | for ea. reference and simulation ID, generates a `.json` file in `meta/run_pars`<br>that specifies runmode parameters (e.g., nb. of cells/genes to sample, nb. of run replicates) 
`meth_pars` | for ea. method ID, generates a `.json` file in `meta/meth_pars`<br>that specifies method parameters
`sim_data`  | provided with a reference dataset and simulation parameters,<br>simulates data and writes a SCE to `data/sim_data`
`apply_X`   | wrapper to run DS method of type X (`distinct`, `pb`, `mm`)
`run_meth`  | reads in simulated data, method parameters, and performs DS analysis<br>by running the corresponding `apply_X` script
`utils`     | various helpers for data handling, formatting, and plotting
`session_info` | generates a `.txt` file capturing the output of `session_info()`


Two steps are necessary before running the snakemake:

1) run the setup.R file;

``` R
source("setup.R")
```

2) create the `sce0_kang.rds` file needed for the analyses.

``` R
library(muscData)
sce = Kang18_8vs8()

library(scater)
library(sctransform)
library(SingleCellExperiment) 
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)
assays(sce)$cpm <- calculateCPM(sce)

vstresiduals = sctransform::vst( assays(sce)$counts, show_progress = FALSE)$y
assays(sce)$vstresiduals = vstresiduals
sum(is.na(assays(sce)$logvstresiduals))

sce

saveRDS(sce, file ="data/raw_data/sce0_kang.rds")
```

The snakemake can be run in parallel (on `n_cores`) via:
``` bash
snakemake --cores n_cores
```