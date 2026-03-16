```
> sessionInfo()

```

```
R version 4.2.0 (2022-04-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Rocky Linux 9.6 (Blue Onyx)

Matrix products: default

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] readxl_1.4.3                writexl_1.5.0              
 [3] openxlsx_4.2.7.1            magrittr_2.0.3             
 [5] data.table_1.18.0           purrr_1.0.4                
 [7] forcats_1.0.0               stringr_1.5.1              
 [9] tidyr_1.3.1                 dplyr_1.1.4                
[11] patchwork_1.3.2             ggrepel_0.9.6              
[13] scales_1.4.0                systemfonts_1.1.0          
[15] svglite_2.1.3               circlize_0.4.16            
[17] ComplexHeatmap_2.14.0       ggplot2_4.0.1              
[19] emmeans_1.10.4              compositions_2.0-8         
[21] msigdbr_25.1.1              clusterProfiler_4.12.6     
[23] harmony_1.2.4               Rcpp_1.1.1                 
[25] scds_1.14.0                 scDblFinder_1.12.0         
[27] SingleCellExperiment_1.20.1 SummarizedExperiment_1.28.0
[29] Biobase_2.58.0              GenomicRanges_1.50.2       
[31] GenomeInfoDb_1.34.9         IRanges_2.32.0             
[33] S4Vectors_0.36.2            BiocGenerics_0.44.0        
[35] MatrixGenerics_1.16.0       matrixStats_1.4.1          
[37] SoupX_1.6.2                 Seurat_5.1.0               
[39] SeuratObject_5.0.2          sp_2.1-4                   

loaded via a namespace (and not attached):
  [1] rsvd_1.0.5                ica_1.0-3                
  [3] Rsamtools_2.14.0          foreach_1.5.2            
  [5] lmtest_0.9-40             crayon_1.5.3             
  [7] MASS_7.3-56               nlme_3.1-168             
  [9] GOSemSim_2.30.2           rlang_1.1.6              
 [11] XVector_0.38.0            HDO.db_0.99.1            
 [13] ROCR_1.0-11               irlba_2.3.5.1            
 [15] limma_3.54.2              scater_1.32.1            
 [17] xgboost_1.7.8.1           BiocParallel_1.32.6      
 [19] rjson_0.2.23              bit64_4.5.2              
 [21] glue_1.8.0                sctransform_0.4.1        
 [23] parallel_4.2.0            vipor_0.4.7              
 [25] spatstat.sparse_3.1-0     AnnotationDbi_1.60.2     
 [27] SeuratDisk_0.0.0.9021     dotCall64_1.1-1          
 [29] DOSE_3.24.2               spatstat.geom_3.3-3      
 [31] tidyselect_1.2.1          fitdistrplus_1.2-1       
 [33] XML_3.99-0.17             zoo_1.8-12               
 [35] spatstat.univar_3.0-1     GenomicAlignments_1.34.1 
 [37] xtable_1.8-4              RcppHNSW_0.6.0           
 [39] scuttle_1.8.4             cli_3.6.5                
 [41] zlibbioc_1.44.0           rstudioapi_0.16.0        
 [43] miniUI_0.1.1.1            fastmatch_1.1-4          
 [45] treeio_1.22.0             fastDummies_1.7.4        
 [47] shiny_1.9.1               BiocSingular_1.14.0      
 [49] xfun_0.47                 clue_0.3-65              
 [51] gson_0.1.0                cluster_2.1.3            
 [53] tidygraph_1.3.1           KEGGREST_1.38.0          
 [55] tibble_3.2.1              ape_5.8                  
 [57] listenv_0.9.1             Biostrings_2.66.0        
 [59] png_0.1-8                 future_1.40.0            
 [61] withr_3.0.1               bitops_1.0-8             
 [63] ggforce_0.4.2             cellranger_1.1.0         
 [65] plyr_1.8.9                dqrng_0.4.1              
 [67] pROC_1.18.5               coda_0.19-4.1            
 [69] pillar_1.9.0              GlobalOptions_0.1.2      
 [71] cachem_1.1.0              multcomp_1.4-26          
 [73] fs_1.6.4                  hdf5r_1.3.12             
 [75] GetoptLong_1.0.5          DelayedMatrixStats_1.20.0
 [77] vctrs_0.6.5               generics_0.1.3           
 [79] tools_4.2.0               beeswarm_0.4.0           
 [81] tweenr_2.0.3              fgsea_1.24.0             
 [83] DelayedArray_0.24.0       fastmap_1.2.0            
 [85] compiler_4.2.0            abind_1.4-8              
 [87] httpuv_1.6.15             rtracklayer_1.58.0       
 [89] plotly_4.10.4             GenomeInfoDbData_1.2.9   
 [91] gridExtra_2.3             edgeR_3.40.2             
 [93] lattice_0.20-45           deldir_2.0-4             
 [95] utf8_1.2.4                later_1.3.2              
 [97] jsonlite_1.8.9            ScaledMatrix_1.6.0       
 [99] tidytree_0.4.6            pbapply_1.7-2            
[101] sparseMatrixStats_1.10.0  estimability_1.5.1       
[103] lazyeval_0.2.2            promises_1.3.0           
[105] doParallel_1.0.17         R.utils_2.12.3           
[107] goftest_1.2-3             spatstat.utils_3.1-0     
[109] reticulate_1.39.0         sandwich_3.1-1           
[111] cowplot_1.1.3             statmod_1.5.0            
[113] Rtsne_0.17                dichromat_2.0-0.1        
[115] uwot_0.2.2                igraph_2.0.3             
[117] survival_3.7-0            yaml_2.3.10              
[119] htmltools_0.5.8.1         memoise_2.0.1            
[121] BiocIO_1.8.0              locfit_1.5-9.10          
[123] graphlayouts_1.2.0        viridisLite_0.4.2        
[125] digest_0.6.37             assertthat_0.2.1         
[127] mime_0.12                 rappdirs_0.3.3           
[129] httr2_1.0.4               spam_2.10-0              
[131] RSQLite_2.3.7             yulab.utils_0.2.0        
[133] future.apply_1.11.2       blob_1.2.4               
[135] R.oo_1.26.0               splines_4.2.0            
[137] RCurl_1.98-1.16           colorspace_2.1-1         
[139] ggbeeswarm_0.7.2          shape_1.4.6.1            
[141] aplot_0.2.3               RANN_2.6.2               
[143] mvtnorm_1.3-1             enrichplot_1.24.4        
[145] fansi_1.0.6               parallelly_1.44.0        
[147] R6_2.5.1                  ggridges_0.5.6           
[149] lifecycle_1.0.4           zip_2.3.1                
[151] bluster_1.8.0             curl_7.0.0               
[153] leiden_0.4.3.1            robustbase_0.99-4        
[155] Matrix_1.6-5              qvalue_2.30.0            
[157] RcppAnnoy_0.0.22          TH.data_1.1-2            
[159] RColorBrewer_1.1-3        iterators_1.0.14         
[161] spatstat.explore_3.3-2    S7_0.2.1                 
[163] htmlwidgets_1.6.4         beachmat_2.14.2          
[165] polyclip_1.10-7           shadowtext_0.1.4         
[167] gridGraphics_0.5-1        globals_0.17.0           
[169] spatstat.random_3.3-2     tensorA_0.36.2.1         
[171] progressr_0.14.0          codetools_0.2-18         
[173] GO.db_3.16.0              metapod_1.6.0            
[175] RSpectra_0.16-2           R.methodsS3_1.8.2        
[177] gtable_0.3.6              DBI_1.2.3                
[179] bayesm_3.1-6              ggfun_0.1.6              
[181] tensor_1.5                httr_1.4.7               
[183] KernSmooth_2.23-20        stringi_1.8.4            
[185] reshape2_1.4.4            farver_2.1.2             
[187] viridis_0.6.5             ggtree_3.12.0            
[189] BiocNeighbors_1.16.0      restfulr_0.0.15          
[191] ggplotify_0.1.2           scattermore_1.2          
[193] scran_1.32.0              DEoptimR_1.1-3           
[195] bit_4.5.0                 scatterpie_0.2.4         
[197] spatstat.data_3.1-2       ggraph_2.2.1             
[199] pkgconfig_2.0.3           babelgene_22.9           
[201] knitr_1.48
```
