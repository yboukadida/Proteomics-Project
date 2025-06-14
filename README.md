# TMT Analysis

This R project analyzes proteomics data from TMT-labeled experiments. It performs preprocessing, statistical testing, and visualizations to identify differentially expressed peptides/proteins. It also includes MS2 spectrum inspection from raw `.mzML` files.

## Data

This pipeline was developed using a TMT‐labeled dataset from ProteomeXchange (e.g., PXD021704, which generated F063721.dat‐mztab.txt and TMT.mzML). Before running, place your own *.mzML and *.mztab files (or the original files) into the data/ folder.


## Files

- 01_main_pipeline.R – Main pipeline: data loading, cleaning, statistical analysis, and visualizations
- 02_peptide_spectrum_link.R – Plots an MS2 spectrum for a significant peptide
- 03_mzml_exploration.R – General exploration of MS2 scans and TIC chromatogram
- 04_correlation_analysis.R – Optional script for correlation studies
  
## Requirements

Install the required packages in R:
```r
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("rpx", "MSnbase", "mzR", "Biobase"))
install.packages(c("ggplot2", "dplyr", "pheatmap", "FactoMineR", "factoextra"))
```

## How to Run

Run the scripts in order:
```r
source("01_main_pipeline.R")
source("02_peptide_spectrum_link.R")
source("03_mzml_exploration.R")
source("04_correlation_analysis.R")  # Optional
```



## Notes

- This project was developed for academic purposes.
- Figures will be generated by the scripts and saved locally.
- Raw input files (`.mzML`, `.mztab`) must be stored in the `data/` folder (not pushed to GitHub).
