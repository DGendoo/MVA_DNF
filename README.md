**MVA-DNF**
========

## Analysis Pipeline and Code Development: 

This MVA-DNF analysis pipeline and code was concieved, designed and developed by **Deena M.A. Gendoo** for the following publication: 

**_Computational pharmacogenomics screen identifies synergistic statin-compound combinations as anti-breast cancer therapies_**

Further Bioinformatics Analysis on Drug combinations is also developed by **Wail B-Alawi** and can be found in the subfolder Wail_ComboAnalysis

**Publication:** 
Jenna van Leeuwen, Wail Ba-Alawi, Emily Branchard, Joseph Longo, Jennifer Silvester, David W. Cescon, Benjamin Haibe-Kains, Linda Z. Penn, Deena M.A. Gendoo. bioRxiv 2020.09.07.286922; doi: https://doi.org/10.1101/2020.09.07.286922 

**Questions or Comments:** 
Please email d.gendoo@bham.ac.uk or deena.gendoo1984@gmail.com

## Introduction to the Analysis

This repository hosts code to analyze 



## The Analysis 

We describe how to reproduce the statistical analysis as reported in the manuscript. To do this, please proceed to:

1. Set up the software environment
2. Run the R scripts

### Set up the software environment

We developed and tested our analysis pipeline using R running on Mac OS X platforms.

To mimic our software environment the following R packages should be installed. All these packages are available on CRAN or Bioconductor.


```
R version 3.3.1 (2016-06-21)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.10.3 (Yosemite)

locale:
[1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8

attached base packages:
 [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] GenVisR_1.0.4              plotrix_3.6-4              vcfR_1.4.0                
 [4] scales_0.4.1               VennDiagram_1.6.17         futile.logger_1.4.3       
 [7] VariantAnnotation_1.18.7   Rsamtools_1.24.0           Biostrings_2.40.2         
[10] XVector_0.12.1             SummarizedExperiment_1.2.3 Biobase_2.32.0            
[13] GenomicRanges_1.24.3       GenomeInfoDb_1.8.7         IRanges_2.6.1             
[16] S4Vectors_0.10.3           copynumber_1.12.0          BiocGenerics_0.18.0       
[19] gplots_3.0.1               RCircos_1.2.0              pheatmap_1.0.8            
[22] reshape2_1.4.2            

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.9             ape_4.0                 lattice_0.20-34         gtools_3.5.0           
 [5] rprojroot_1.2           assertthat_0.1          digest_0.6.11           plyr_1.8.4             
 [9] backports_1.0.5         futile.options_1.0.0    evaluate_0.10           RSQLite_1.1-2          
[13] ggplot2_2.2.1           zlibbioc_1.18.0         GenomicFeatures_1.24.5  lazyeval_0.2.0         
[17] gdata_2.17.0            vegan_2.4-2             Matrix_1.2-7.1          rmarkdown_1.6          
[21] pinfsc50_1.1.0          BiocParallel_1.6.6      stringr_1.2.0           RCurl_1.95-4.8         
[25] biomaRt_2.28.0          munsell_0.4.3           rtracklayer_1.32.2      mgcv_1.8-16            
[29] htmltools_0.3.5         gridExtra_2.2.1         tibble_1.2              XML_3.98-1.5           
[33] permute_0.9-4           viridisLite_0.1.3       GenomicAlignments_1.8.4 MASS_7.3-45            
[37] bitops_1.0-6            nlme_3.1-128            gtable_0.2.0            DBI_0.5-1              
[41] magrittr_1.5            KernSmooth_2.23-15      stringi_1.1.2           viridis_0.3.4          
[45] lambda.r_1.1.9          RColorBrewer_1.1-2      tools_3.3.1             FField_0.1.0           
[49] BSgenome_1.40.1         AnnotationDbi_1.34.4    colorspace_1.3-2        cluster_2.0.5          
[53] caTools_1.17.1          memoise_1.0.0           knitr_1.15.1                         

```

### Running the R Scripts

Once the packages are installed, please download this github repository. 

This repository contains three folders: 
