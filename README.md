This repository contains the code used for the analyses presented in the manuscript "distinct: a novel approach to differential distribution analyses", which presents the differential distribution R/Bioconductor package distinct.

The scripts are organized in five folders:

- `CDF plots` contains the code to plot the differential profiles;

- `diffcyt semi-simulated data` contains the code of the `diffcyt` simulation study;

- `Kang dataset - exploratory plots` contains the code to perform some exploratory plots of the `Kang` dataset;

- `muscat simulation and Kang real data` contains the Snakemakes and the code to perform all `muscat` simulation studies and the `Kang` non-null real data analysis (comparing controls to stimulated samples);

- `NULL experimental data` contains the code to run the null analyses on the real `T-cells` and `Kang` datasets.

The Snakemakes in `muscat simulation and Kang real data` were obtained by editing the Snakemake from: https://github.com/HelenaLC/muscat-comparison .
The folder `muscat simulation and Kang real data` also contains a modified version of muscat R/Bioconductor package, developed by Almut Luetge at the Robinson lab (University of Zurich), that allows simulating batch effects.
