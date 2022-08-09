**MVA-DNF**
========

## Analysis Pipeline and Code Development: 

This MVA-DNF analysis pipeline and code was concieved, designed and developed by **Deena M.A. Gendoo** for the following publication: 

**_Computational pharmacogenomics screen identifies synergistic statin-compound combinations as anti-breast cancer therapies_**

Further Bioinformatics Analysis on Drug combinations is also developed by **Wail B-Alawi** and can be found in the subfolder Wail_ComboAnalysis

**Questions or Comments:** 
Please email d.gendoo@bham.ac.uk or deena.gendoo1984@gmail.com

## Publication: 

**Please cite:** 
Jenna van Leeuwen, Wail Ba-Alawi, Emily Branchard, Joseph Longo, Jennifer Silvester, David W. Cescon, Benjamin Haibe-Kains, Linda Z. Penn, Deena M.A. Gendoo. bioRxiv 2020.09.07.286922; doi: https://doi.org/10.1101/2020.09.07.286922 


## Introduction to the Analysis

This repository hosts code for the MVA-DNF algorithm, using a pathway-centric approach to identify top drug agents, like dipyridamole, that potentiate statin-induced tumor cell death by targeting the mevalonate pathway. 

An integrative pharmacogenomics pipeline has been developed to identify agents that were similar to dipyridamole at the level of drug structure, in vitro sensitivity and molecular perturbation, while enriching for compounds expected to target the mevalonate pathway. This resulted in the MVA-DNF (mevalonate drug network fusion). 

Further code assesses the synergistic ability of the top DP-like compounds with Fluvastatin, and their effect on the mevalonate pathway. 

## The Analysis 

We describe how to reproduce the statistical analysis as reported in the manuscript. To do this, please proceed to:

1. Set up the software environment
2. Run the R scripts

### Set up the software environment

We developed and tested our analysis pipeline using R running on Mac OS X platforms.

To mimic our software environment the following R packages should be installed. All these packages are available on CRAN or Bioconductor.


```
R version 3.6.0 (2019-04-26)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] RColorBrewer_1.1-2   gplots_3.0.1.1       ggplot2_3.3.2        GSVA_1.32.0          GSA_1.03.1          
 [6] piano_2.0.2          snowfall_1.84-6.1    snow_0.4-3           proxy_0.4-24         reshape2_1.4.3      
[11] survcomp_1.34.0      prodlim_2019.10.13   survival_3.1-7       ROCR_1.0-11          SNFtool_2.3.0       
[16] org.Hs.eg.db_3.8.2   annotate_1.62.0      XML_3.98-1.20        AnnotationDbi_1.46.1 IRanges_2.18.3      
[21] S4Vectors_0.22.1     Biobase_2.44.0       BiocGenerics_0.30.0  fingerprint_3.5.7    rcdk_3.5.0          
[26] rcdklibs_2.3         rJava_0.9-13         apcluster_1.4.8      PharmacoGx_1.14.2    fmsb_0.6.3          
[31] pheatmap_1.0.12     

loaded via a namespace (and not attached):
  [1] fgsea_1.10.1         colorspace_1.4-1     ellipsis_0.3.0       lsa_0.73.1           rstudioapi_0.10     
  [6] SnowballC_0.6.0      DT_0.10              bit64_0.9-7          splines_3.6.0        geneplotter_1.62.0  
 [11] shinythemes_1.1.2    SuppDists_1.1-9.4    heatmap.plus_1.3     itertools_0.1-3      jsonlite_1.6        
 [16] magicaxis_2.0.7      alluvial_0.1-2       cluster_2.1.0        png_0.1-7            graph_1.62.0        
 [21] shinydashboard_0.7.1 shiny_1.4.0          ExPosition_2.8.23    mapproj_1.2.6        compiler_3.6.0      
 [26] prettyGraphs_2.1.6   assertthat_0.2.1     Matrix_1.2-17        fastmap_1.0.1        limma_3.40.6        
 [31] later_1.0.0          visNetwork_2.0.8     htmltools_0.4.0      tools_3.6.0          igraph_1.2.4.1      
 [36] gtable_0.3.0         glue_1.3.1           RANN_2.6.1           dplyr_0.8.3          maps_3.3.0          
 [41] fastmatch_1.1-0      Rcpp_1.0.3           slam_0.1-46          vctrs_0.3.1          gdata_2.18.0        
 [46] iterators_1.0.12     stringr_1.4.0        mime_0.7             lifecycle_0.2.0      gtools_3.8.1        
 [51] MASS_7.3-51.4        scales_1.1.1         promises_1.1.0       relations_0.6-9      sets_1.0-18         
 [56] memoise_1.1.0        gridExtra_2.3        downloader_0.4       rmeta_3.0            stringi_1.4.3       
 [61] RSQLite_2.1.2        NISTunits_1.0.1      plotrix_3.7-6        caTools_1.17.1.2     BiocParallel_1.18.1 
 [66] lava_1.6.6           rlang_0.4.6          pkgconfig_2.0.3      bitops_1.0-6         pracma_2.2.5        
 [71] lattice_0.20-38      purrr_0.3.3          survivalROC_1.0.3    htmlwidgets_1.5.1    bit_1.1-14          
 [76] tidyselect_1.1.0     GSEABase_1.46.0      plyr_1.8.6.9000      magrittr_1.5         R6_2.4.0            
 [81] bootstrap_2019.6     DBI_1.0.0            sm_2.2-5.6           withr_2.1.2          pillar_1.4.4        
 [86] RCurl_1.95-4.12      tibble_3.0.1         crayon_1.3.4         KernSmooth_2.23-16   grid_3.6.0          
 [91] data.table_1.12.6    marray_1.62.0        blob_1.2.0           digest_0.6.22        xtable_1.8-4        
 [96] httpuv_1.5.2         munsell_0.5.0        celestial_1.4.6      tcltk_3.6.0          shinyjs_1.0           

```

### Running the R Scripts

Once the packages are installed, please download this github repository to run the code. 

The **main folder** contains scripts to run the MVA-DNF algorithm and associated output of the manuscript for:  
*Figure 1*   
*Supplementary Figure S1*   
*Supplementary Table1*  

1. The scripts  **DeenaGendoo_Generate_MVA_DNF.R** and **DeenaGendoo_PermutationTestAndFiltering.R** are used to generate the MVA-DNF matrix, and then identify top drug agents to Dipryidamole, using permutation testing
2. The script **DeenaGendoo_Heatmap_DrugPertSigs.R** generate heatmaps to show up/down regulated genes due to drug treatment (drug perturbation signatures) for Dipryidamole (DP) and DP-like drugs
3. The script **DeenaGendoo_CompareLayerContributions_Dec2021.R** is used to compare for any two drugs, whether the strength of the drug-drug relationships is a reflection of perturbation, sensitivity, or structural similarity. Specific files (Cytoscape files for Figure 1B) are also found in the Data folder. 


The **Wail_ComboAnalysis** subfolder is used to generate associated output of the manuscript for:  
*Figure 4*  
*Supplementary Figure S6*  
*Supplementary Figure S7*  
*Supplementary Figure S8*    

1. The script **combo_analysis.R** contains the analysis of the synergy between the identified DP-like drugs and Fluvastatin.  
