# PCFSiM-Paper

This repository contains all scripts used in the paper "Exploiting pair correlation function to describe biological tissue structure". Here we use functions located in the "List_function_pcf_analysis.R" file script to perform PCF-SiM analysis but we recommend to use the R [PCFSiM package](https://github.com/BOSTLAB/PCFSiM).
Original and processed datasets can be found on a [Zenodo repository](https://zenodo.org/records/17867607). 
The repository is splitted into 4 directories:
- **Processing_annotation** which contains scripts required to annoate the single-cell spatial datasets but also to process the raw IMC .mcd files for the thyroid study.
- **Spatial_analysis_healthy** which contains scripts required for the spatial analysis of the healthy samples, i.e. samples analyzed in Figure 1.
- **Spatial_analysis_clinical** which contains scripts required for the spatial analysis of clinical samples, i.e. samples analyzed in Figure 3, 4, 5 and 6.
- **Other_scripts** which corresponds to all other scripts, notably the ones used for simulations (Figure 1), cross-pcf (Figure S4) and analysis of non-cellular patterns (Figure 2).
