# Persistence of mucosal CAR-T cells and inflammatory remodeling in enterocolitis associated with BCMA CAR-T cell therapy

* Repository contains code used to analyze single-cell RNA-seq data published in [Nature Medicine](https://www.nature.com/articles/s41591-026-04632-y). Alternatively, PMID for publication is: PMID 42778764; [PubMed URL.](https://pubmed.ncbi.nlm.nih.gov/42778764/)<br>
* Raw single-cell data available under GEO accession number [GSE322796.](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE322796)<br>
* Per-sample CellChat data (prior to modeling) is available on FigShare: [LR Interaction Tables.](https://doi.org/10.6084/m9.figshare.32600346)

## Software & Packages
- R was used for the core single-cell workflow (QC through CAR+ mapping and CellChat).
- Python was only used for downstream modeling of the per-sample CellChat outputs (per-interaction regression, rank aggregation, and EDF7 and EDF8 figures).

<table>
  <tr>
    <td valign="top">
<b>Single-Cell Analysis</b><br>
<a href="https://www.r-project.org/"><img src="https://img.shields.io/badge/R-v4.2.0-4e9af1?style=flat-square&logo=r&logoColor=white"/></a> 
<a href="https://satijalab.org/seurat/"><img src="https://img.shields.io/badge/Seurat-v5.1.0-0288d1?style=flat-square"/></a> 
<a href="https://github.com/constantAmateur/SoupX"><img src="https://img.shields.io/badge/SoupX-v1.6.2-1a73e8?style=flat-square"/></a>
<a href="https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html"><img src="https://img.shields.io/badge/scDblFinder-v1.12.0-2962ff?style=flat-square"/></a>
<a href="https://portals.broadinstitute.org/harmony/"><img src="https://img.shields.io/badge/harmony-v1.2.4-276DC3?style=flat-square"/></a><br><br>
<b>Pathway Enrichment</b><br>
<a href="https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html"><img src="https://img.shields.io/badge/clusterProfiler-v4.12.6-26a69a?style=flat-square"/></a>
<a href="https://igordot.github.io/msigdbr/"><img src="https://img.shields.io/badge/msigdbr-v25.1.1-00897b?style=flat-square"/></a><br><br>
<b>Visualisation</b><br>
<a href="https://ggplot2.tidyverse.org/"><img src="https://img.shields.io/badge/ggplot2-v4.0.1-e91e63?style=flat-square"/></a>
<a href="https://jokergoo.github.io/ComplexHeatmap-reference/book/"><img src="https://img.shields.io/badge/ComplexHeatmap-v2.14.0-c2185b?style=flat-square"/></a>
<a href="https://svglite.r-lib.org/"><img src="https://img.shields.io/badge/svglite-v2.1.3-ad1457?style=flat-square"/></a>
    </td>
    <td valign="top">
<b>Cell Frequency</b><br>
<a href="https://www.stat.boogaart.de/compositions/"><img src="https://img.shields.io/badge/compositions-v2.0.8-f57c00?style=flat-square"/></a>
<a href="https://rvlenth.github.io/emmeans/"><img src="https://img.shields.io/badge/emmeans-v1.10.4-ef6c00?style=flat-square"/></a><br><br>
<b>Cell-Cell Communication</b><br>
<i>R</i><br>
<a href="https://github.com/jinworks/CellChat"><img src="https://img.shields.io/badge/CellChat-v2.1.2-b39ddb?style=flat-square"/></a><br>
<i>Python</i><br>
<a href="https://www.python.org/"><img src="https://img.shields.io/badge/Python-v3.8.16-9575cd?style=flat-square&logo=python&logoColor=white"/></a>
<a href="https://www.statsmodels.org/"><img src="https://img.shields.io/badge/statsmodels-v0.14.0-7e57c2?style=flat-square"/></a>
<a href="https://scipy.org/"><img src="https://img.shields.io/badge/SciPy-v1.10.1-5e35b1?style=flat-square"/></a>
<a href="https://matplotlib.org/"><img src="https://img.shields.io/badge/matplotlib-v3.7.1-4527a0?style=flat-square"/></a>
    </td>
  </tr>
</table>


## Repository Structure

```
CART-Enterocolitis-Analysis-Scripts/
├── R Scripts/
│   ├── 01_QC.Rmd                               # Ambient RNA removal, doublet detection, QC filtering
│   ├── 02_Integration.Rmd                      # Batch correction
│   ├── 03_Clustering.Rmd                       # Global clustering, annotation, condition UMAPs
│   ├── 04_Subclustering_Stroma.Rmd             # Sub-clustering workflow (stromal shown; same workflow applied to T cells, myeloid, B cells, plasma cells, epithelial, and endothelial subsets)
│   ├── 05_DEGs_and_GSEA_Stroma_Condition.Rmd   # Condition-wise DEGs + GSEA
│   ├── 06_CLR_CellFreq_Stroma.Rmd              # CLR cell-frequency analysis + forest plots (stromal shown; same workflow applied globally  and to all other subsets)
│   ├── 07_Figs_Stroma.Rmd                      # Figures for stromal subset (bubble plots, stacked bars, condition plots, etc.)
│   └── 08_CAR_Mapping_and_Figs.Rmd             # Certomics CAR+ cell mapping onto T cell subset, heatmaps, pie plots, etc.
│   └── 09_CellChat.Rmd                         # Sample-wise CellChat that was aggregated and analyzed further using regression model
│   └── sessionInfo.md                          # Full R session information and package versions                     
├── Python Scripts/
│   ├── colitis_analysis_script.ipynb           # End-to-end CellChat downstream analysis notebook
│   └── clplots/                                # Plotting helpers (volcano, circle, RRA/regression boxplots, etc.)
│   └── clutils/                                # Data preparation, ranking, summarization, and regression utilities
│   └── package_versions.txt                    # Full Jupyter session information and package versions
└── README.md
```
