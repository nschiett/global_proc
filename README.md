
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Drivers and vulnerability of global coral reef fish functions

<!-- badges: start -->

<!-- badges: end -->

The goal of this R project is to reproduce all analyses, figures, and
tables of a global analyses on reef fish function. The project uses a
drake workflow, which binds all elements of the analyses together with a
plan.

## Content

The directory contains:

  - [:file\_folder: R](/R): Folder containing all functions, packages
    and drake plan.  
    \- **plan.R**: File containing the drake plan.  
    \- **packages.R**: File to load all needed packages.  
    \- **functions\_wrangling.R**: Contains functions to prepare the
    data for analyses.  
    \- **functions\_analysis.R**: Contains functions to reproduce
    analysis.  
    \- **functions\_plots.R**: Contains functions to create figures
    presenting the results.  
  - [:file\_folder: data](/data): Folder containing all data.
  - [:file\_folder: text](/text): All files to create output documents
    for publication.
  - make.R File: Script to run the whole project.

## Working environment

``` r
sessionInfo()
#> R version 3.6.3 (2020-02-29)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 18.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
#> 
#> locale:
#>  [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8    
#>  [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=en_US.utf8    
#>  [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] officer_0.3.14    flextable_0.5.11  broom_0.5.2      
#>  [4] readr_1.3.1       ggrepel_0.8.1     forcats_0.4.0    
#>  [7] ggtree_2.0.1      tidytree_0.3.0    ggnewscale_0.4.0 
#> [10] fishualize_0.2.0  patchwork_1.0.0   brms_2.13.5      
#> [13] Rcpp_1.0.5        e1071_1.7-2       drake_7.12.5     
#> [16] tidyr_1.0.0       dplyr_0.8.3       Rphylopars_0.2.11
#> [19] ape_5.3           fishflux_0.0.1.1  ggplot2_3.2.1    
#> [22] tidybayes_1.1.0   purrr_0.3.3      
#> 
#> loaded via a namespace (and not attached):
#>   [1] uuid_0.1-4                backports_1.1.5          
#>   [3] fastmatch_1.1-0           systemfonts_0.1.1        
#>   [5] plyr_1.8.4                igraph_1.2.4.1           
#>   [7] lazyeval_0.2.2            svUnit_0.7-12            
#>   [9] storr_1.2.1               crosstalk_1.0.0          
#>  [11] listenv_0.7.0             rstantools_2.1.1         
#>  [13] inline_0.3.15             digest_0.6.22            
#>  [15] htmltools_0.4.0           rsconnect_0.8.15         
#>  [17] magrittr_1.5              phytools_0.6-99          
#>  [19] memoise_1.1.0             base64url_1.4            
#>  [21] globals_0.12.4            RcppParallel_5.0.2       
#>  [23] matrixStats_0.55.0        xts_0.11-2               
#>  [25] prettyunits_1.0.2         colorspace_1.4-1         
#>  [27] xfun_0.10                 callr_3.4.4              
#>  [29] crayon_1.3.4              jsonlite_1.6             
#>  [31] zoo_1.8-6                 phangorn_2.5.5           
#>  [33] glue_1.3.1                gtable_0.3.0             
#>  [35] geiger_2.0.6.2            V8_3.2.0                 
#>  [37] pkgbuild_1.0.6            rstan_2.21.2             
#>  [39] future.apply_1.3.0        maps_3.3.0               
#>  [41] abind_1.4-5               scales_1.1.0             
#>  [43] mvtnorm_1.0-11            miniUI_0.1.1.1           
#>  [45] plotrix_3.7-6             xtable_1.8-4             
#>  [47] progress_1.2.2            ggstance_0.3.3           
#>  [49] subplex_1.5-4             txtq_0.2.3               
#>  [51] deSolve_1.24              stats4_3.6.3             
#>  [53] StanHeaders_2.21.0-6      DT_0.9                   
#>  [55] animation_2.6             htmlwidgets_1.5.1        
#>  [57] httr_1.4.1                threejs_0.3.1            
#>  [59] arrayhelpers_1.0-20160527 pkgconfig_2.0.3          
#>  [61] loo_2.3.1                 tidyselect_1.1.0         
#>  [63] rlang_0.4.7               reshape2_1.4.3           
#>  [65] later_1.0.0               munsell_0.5.0            
#>  [67] tools_3.6.3               cli_1.1.0                
#>  [69] generics_0.0.2            ggridges_0.5.1           
#>  [71] evaluate_0.14             stringr_1.4.0            
#>  [73] fastmap_1.0.1             yaml_2.2.0               
#>  [75] processx_3.4.1            knitr_1.25               
#>  [77] zip_2.1.1                 gh_1.0.1                 
#>  [79] future_1.14.0             nlme_3.1-149             
#>  [81] mime_0.7                  xml2_1.2.2               
#>  [83] doBy_4.6-2                rfishbase_3.0.4          
#>  [85] compiler_3.6.3            bayesplot_1.7.0          
#>  [87] shinythemes_1.1.2         filelock_1.0.2           
#>  [89] curl_4.2                  png_0.1-7                
#>  [91] treeio_1.10.0             clusterGeneration_1.3.4  
#>  [93] tibble_2.1.3              stringi_1.4.3            
#>  [95] ps_1.3.0                  Brobdingnag_1.2-6        
#>  [97] gdtools_0.2.1             lattice_0.20-41          
#>  [99] Matrix_1.2-18             markdown_1.1             
#> [101] shinyjs_1.0               vctrs_0.3.4              
#> [103] pillar_1.4.2              lifecycle_0.1.0          
#> [105] BiocManager_1.30.8        combinat_0.0-8           
#> [107] bridgesampling_0.7-2      data.table_1.12.8        
#> [109] httpuv_1.5.2              R6_2.4.0                 
#> [111] promises_1.1.0            gridExtra_2.3            
#> [113] mvnmle_0.1-11.1           codetools_0.2-16         
#> [115] colourpicker_1.0          MASS_7.3-53              
#> [117] gtools_3.8.1              assertthat_0.2.1         
#> [119] withr_2.1.2               phylolm_2.6              
#> [121] shinystan_2.5.0           mnormt_1.5-5             
#> [123] expm_0.999-4              parallel_3.6.3           
#> [125] hms_0.5.1                 quadprog_1.5-7           
#> [127] grid_3.6.3                coda_0.19-3              
#> [129] class_7.3-17              rvcheck_0.1.5            
#> [131] rmarkdown_2.1             numDeriv_2016.8-1.1      
#> [133] scatterplot3d_0.3-41      shiny_1.4.0              
#> [135] base64enc_0.1-3           dygraphs_1.1.1.6
```
