# R package "amplicon" for amplicon data statistics and visualization

Statistics and visualization for amplicon data

![image](http://210.75.224.110/Note/R/amplicon/fig2.png)

## Figure 2. Examples of visualizations. 

A. Boxplot showing richness of alpha diversity. B. Rarefaction curve of richness. C. Principal coordinate analysis (PCoA) of Bray-Curtis distance. D. Heatmap of Bray-Curtis distance. E. Stack plot of taxonomic composition in phylum level. F. Tree map of taxonomic composition. G. Volcano plot showing count per million and fold change of different features between KO and WT groups. H. Manhattan plot showing different features and related taxonomy between KO and WT groups.


Main features:

- Input feature (OTU/ASV/Taxonomic) tables, metadata;
- Supports more than ten analysis methods, provides multiple visualization styles;
    - Alpha diversity: index barplot, boxplot with stat label;
    - Beta diveristy: Supported 47 distance * 8 ordination methods, output scatter plot with stat eclippse and sample label;
    - Taxonomic composition: Stackplot in samples or groups, even support bubble plot.
- Published-ready figures and tables.

![image](http://210.75.224.110/Note/R/amplicon/fig1.png)

## Figure 1. Pipeline of EasyAmplicon for analysis pair-end amplicon data.
The EasyAmplicon pipeline mainly include three steps. Dimensionality reduction: process raw sequencing into feature table; Analysis: provide functional prediction, alpha, and beta diversity


## Install

Please open RStudio， and run the following code for install.

    library(devtools)
    install_github("microbiota/amplicon")