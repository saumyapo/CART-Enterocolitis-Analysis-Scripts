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
 [1] CellChat_2.1.2              synchronicity_1.3.10       
 [3] bigmemory_4.6.4             igraph_2.0.3               
 [5] readxl_1.4.3                writexl_1.5.0              
 [7] openxlsx_4.2.7.1            magrittr_2.0.3             
 [9] data.table_1.18.0           purrr_1.0.4                
[11] forcats_1.0.0               stringr_1.5.1              
[13] tidyr_1.3.1                 dplyr_1.1.4                
[15] patchwork_1.3.2             ggrepel_0.9.6              
[17] scales_1.4.0                systemfonts_1.1.0          
[19] svglite_2.1.3               circlize_0.4.16            
[21] ComplexHeatmap_2.14.0       ggplot2_4.0.1              
[23] emmeans_1.10.4              compositions_2.0-8         
[25] msigdbr_25.1.1              clusterProfiler_4.12.6     
[27] harmony_1.2.4               Rcpp_1.1.1                 
[29] scDblFinder_1.12.0          SingleCellExperiment_1.20.1
[31] SummarizedExperiment_1.28.0 Biobase_2.58.0             
[33] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
[35] IRanges_2.32.0              S4Vectors_0.36.2           
[37] BiocGenerics_0.44.0         MatrixGenerics_1.16.0      
[39] matrixStats_1.4.1           SoupX_1.6.2                
[41] Seurat_5.1.0                SeuratObject_5.0.2         
[43] sp_2.1-4                   

loaded via a namespace (and not attached):
  [1] statnet.common_4.9.0      rsvd_1.0.5               
  [3] ica_1.0-3                 Rsamtools_2.14.0         
  [5] foreach_1.5.2             lmtest_0.9-40            
  [7] crayon_1.5.3              MASS_7.3-56              
  [9] backports_1.5.0           nlme_3.1-168             
 [11] GOSemSim_2.30.2           rlang_1.1.6              
 [13] XVector_0.38.0            HDO.db_0.99.1            
 [15] ROCR_1.0-11               irlba_2.3.5.1            
 [17] limma_3.54.2              scater_1.32.1            
 [19] xgboost_1.7.8.1           BiocParallel_1.32.6      
 [21] rjson_0.2.23              bit64_4.5.2              
 [23] glue_1.8.0                rngtools_1.5.2           
 [25] sctransform_0.4.1         parallel_4.2.0           
 [27] vipor_0.4.7               spatstat.sparse_3.1-0    
 [29] AnnotationDbi_1.60.2      SeuratDisk_0.0.0.9021    
 [31] dotCall64_1.1-1           DOSE_3.24.2              
 [33] spatstat.geom_3.3-3       tidyselect_1.2.1         
 [35] fitdistrplus_1.2-1        XML_3.99-0.17            
 [37] zoo_1.8-12                ggpubr_0.6.0             
 [39] spatstat.univar_3.0-1     GenomicAlignments_1.34.1 
 [41] ggnetwork_0.5.13          xtable_1.8-4             
 [43] RcppHNSW_0.6.0            scuttle_1.8.4            
 [45] cli_3.6.5                 zlibbioc_1.44.0          
 [47] rstudioapi_0.16.0         miniUI_0.1.1.1           
 [49] bslib_0.8.0               fastmatch_1.1-4          
 [51] treeio_1.22.0             fastDummies_1.7.4        
 [53] shiny_1.9.1               BiocSingular_1.14.0      
 [55] xfun_0.47                 clue_0.3-65              
 [57] gson_0.1.0                cluster_2.1.3            
 [59] tidygraph_1.3.1           KEGGREST_1.38.0          
 [61] tibble_3.2.1              ape_5.8                  
 [63] listenv_0.9.1             Biostrings_2.66.0        
 [65] png_0.1-8                 future_1.40.0            
 [67] withr_3.0.1               bitops_1.0-8             
 [69] ggforce_0.4.2             cellranger_1.1.0         
 [71] plyr_1.8.9                dqrng_0.4.1              
 [73] coda_0.19-4.1             pillar_1.9.0             
 [75] GlobalOptions_0.1.2       cachem_1.1.0             
 [77] multcomp_1.4-26           fs_1.6.4                 
 [79] hdf5r_1.3.12              GetoptLong_1.0.5         
 [81] DelayedMatrixStats_1.20.0 vctrs_0.6.5              
 [83] generics_0.1.3            NMF_0.28                 
 [85] tools_4.2.0               beeswarm_0.4.0           
 [87] tweenr_2.0.3              fgsea_1.24.0             
 [89] DelayedArray_0.24.0       fastmap_1.2.0            
 [91] compiler_4.2.0            abind_1.4-8              
 [93] httpuv_1.6.15             rtracklayer_1.58.0       
 [95] plotly_4.10.4             GenomeInfoDbData_1.2.9   
 [97] gridExtra_2.3             edgeR_3.40.2             
 [99] lattice_0.20-45           deldir_2.0-4             
[101] utf8_1.2.4                later_1.3.2              
[103] jsonlite_1.8.9            ScaledMatrix_1.6.0       
[105] carData_3.0-5             tidytree_0.4.6           
[107] pbapply_1.7-2             sparseMatrixStats_1.10.0 
[109] estimability_1.5.1        lazyeval_0.2.2           
[111] promises_1.3.0            car_3.1-2                
[113] doParallel_1.0.17         R.utils_2.12.3           
[115] goftest_1.2-3             sna_2.8                  
[117] spatstat.utils_3.1-0      reticulate_1.39.0        
[119] sandwich_3.1-1            cowplot_1.1.3            
[121] statmod_1.5.0             Rtsne_0.17               
[123] dichromat_2.0-0.1         uwot_0.2.2               
[125] survival_3.7-0            yaml_2.3.10              
[127] htmltools_0.5.8.1         memoise_2.0.1            
[129] BiocIO_1.8.0              locfit_1.5-9.10          
[131] graphlayouts_1.2.0        viridisLite_0.4.2        
[133] digest_0.6.37             assertthat_0.2.1         
[135] mime_0.12                 rappdirs_0.3.3           
[137] httr2_1.0.4               registry_0.5-1           
[139] bigmemory.sri_0.1.8       spam_2.10-0              
[141] RSQLite_2.3.7             yulab.utils_0.2.0        
[143] future.apply_1.11.2       blob_1.2.4               
[145] R.oo_1.26.0               splines_4.2.0            
[147] RCurl_1.98-1.16           broom_1.0.6              
[149] colorspace_2.1-1          BiocManager_1.30.25      
[151] ggbeeswarm_0.7.2          shape_1.4.6.1            
[153] aplot_0.2.3               sass_0.4.9               
[155] RANN_2.6.2                mvtnorm_1.3-1            
[157] enrichplot_1.24.4         fansi_1.0.6              
[159] parallelly_1.44.0         R6_2.5.1                 
[161] ggridges_0.5.6            lifecycle_1.0.4          
[163] zip_2.3.1                 bluster_1.8.0            
[165] ggsignif_0.6.4            curl_7.0.0               
[167] jquerylib_0.1.4           leiden_0.4.3.1           
[169] robustbase_0.99-4         Matrix_1.6-5             
[171] qvalue_2.30.0             RcppAnnoy_0.0.22         
[173] TH.data_1.1-2             RColorBrewer_1.1-3       
[175] iterators_1.0.14          spatstat.explore_3.3-2   
[177] S7_0.2.1                  htmlwidgets_1.6.4        
[179] beachmat_2.14.2           polyclip_1.10-7          
[181] network_1.18.2            shadowtext_0.1.4         
[183] gridGraphics_0.5-1        globals_0.17.0           
[185] spatstat.random_3.3-2     tensorA_0.36.2.1         
[187] progressr_0.14.0          codetools_0.2-18         
[189] FNN_1.1.4.1               GO.db_3.16.0             
[191] metapod_1.6.0             gridBase_0.4-7           
[193] RSpectra_0.16-2           R.methodsS3_1.8.2        
[195] gtable_0.3.6              DBI_1.2.3                
[197] bayesm_3.1-6              ggalluvial_0.12.5        
[199] ggfun_0.1.6               tensor_1.5               
[201] httr_1.4.7                KernSmooth_2.23-20       
[203] stringi_1.8.4             uuid_1.2-1               
[205] reshape2_1.4.4            farver_2.1.2             
[207] viridis_0.6.5             ggtree_3.12.0            
[209] BiocNeighbors_1.16.0      restfulr_0.0.15          
[211] ggplotify_0.1.2           scattermore_1.2          
[213] scran_1.32.0              DEoptimR_1.1-3           
[215] bit_4.5.0                 scatterpie_0.2.4         
[217] spatstat.data_3.1-2       ggraph_2.2.1             
[219] pkgconfig_2.0.3           babelgene_22.9           
[221] rstatix_0.7.2             knitr_1.48
```
