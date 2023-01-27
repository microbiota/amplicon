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
    install_github("liaochenlanruo/amplicon")
    
Version: 1.14.2
Date: 2021-12-12
Author: Yong-Xin Liu, Tao Wen, Tong Chen
Maintainer: Yong-Xin Liu <yxliu@genetics.ac.cn>; Tao Wen <2018203048@njau.edu.cn>; Tong Chen <chent@nrc.ac.cn>
Description: A basic statistics and plotting pipeline for amplicon data.
    1. Alpha boxplot with ANOVA + Tukey test / LSD.test;
    2. Alpha rarefraction curve in samples and groups with standard error;
    3. PCoA plot with confidence ellipse and Adonis P-value;
    4. Constrained PCoA with variant, P-value and confidence ellipse;
    5. Taxonomy stackplot in group and sample.
    1.0.1 add beta_cpcoa_dis.R support constrained distance matrix
    1.1.0 add more than 10 functions
    1.1.1 2020-06-15 Update citation and beta_pcao example
    1.1.2 2020-06-20 Add alpha_sample_rare plot sample rarefaction boxplot
    1.1.3 2020-08-03 Debug code in R 4.0.2, such tax_stackplot
    1.1.4 2020-09-03 Update tax_wordcloud, add tax_stack_clust for tree+stackplot
    1.1.5 2020-09-05 Fix bug "下标出界" in beta_pcoa_stat.R, add intersect metadata and distance matrix 
    1.11.0 2021-04-03 Fix bug metadata only 2 columns cause "Error in metadata[, groupID] : 量度数目不对" 
    1.11.1 2021-04-04 Using multcomp 解决字母顺序相反的问题(alpha_boxplot);tax_stackplot中修改Unassigned至Other组;compare增加不标准化的normlize参数
    1.14.1 2021-10-16 陈同新增堆叠柱状图分面功能
    1.14.2 2021-12-12 文涛更新tax_treemap.R
