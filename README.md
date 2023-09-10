# makino

This repository contains R-codes used in Makino et al's paper (in preparation). All analyses were performed using R ver. 4.2.3.

## Prerequisite 
To execute the codes, [edgeR](https://bioconductor.org/packages/edgeR/), [DESeq2](https://bioconductor.org/packages/DESeq2/), [TCC](https://bioconductor.org/packages/TCC/), [ROC](https://bioconductor.org/packages/ROC/), [compcodeR](https://bioconductor.org/packages/compcodeR/), [PROPER](https://bioconductor.org/packages/PROPER/), [recount](https://bioconductor.org/packages/recount/), [stringr](https://bioconductor.org/packages/stringr/), and [MBCluster.Seq](https://CRAN.R-project.org/package=MBCluster.Seq) are needed:  
```r
if (!requireNamespace("BiocManager", quietly=T))
    install.packages("BiocManager")
BiocManager::install("edgeR", update=F)
BiocManager::install("DESeq2", update=F)
BiocManager::install("TCC", update=F)
BiocManager::install("ROC", update=F)
BiocManager::install("compcodeR", update=F)
BiocManager::install("PROPER", update=F)
BiocManager::install("recount", update=F)
BiocManager::install("stringr", update=F)
install.packages("https://cran.r-project.org/src/contrib/Archive/MBCluster.Seq/MBCluster.Seq_1.0.tar.gz")
```

###  functions.R  ###
This file contains functions to execute the four methods ([edgeR](https://bioconductor.org/packages/edgeR/), [DESeq2](https://bioconductor.org/packages/DESeq2/), [TCC](https://bioconductor.org/packages/TCC/)-DE, and [MBCdeg](https://pubmed.ncbi.nlm.nih.gov/34670485/)) in a unified manner.

###  simu1_tcc.R  ###
This file contains scripts to perform the [TCC](https://bioconductor.org/packages/TCC/) simulation analysis. It includes (1) generation of TCC simulated data, (2) execution of differential expression analysis, and (3) calculation of evaluation metrics (i.e., AUC, Accuracy, and so on).

###  simu2_com.R  ###
This file contains scripts to perform the [compcodeR](https://bioconductor.org/packages/compcodeR/) simulation analysis. It includes (1) generation of compcodeR simulated data, (2) execution of differential expression analysis, and (3) calculation of evaluation metrics (i.e., AUC, Accuracy, and so on).

###  simu3_pro.R  ###
This file contains scripts to perform the [PROPER](https://bioconductor.org/packages/PROPER/) simulation analysis. It includes (1) generation of PROPER simulated data, (2) execution of differential expression analysis, and (3) calculation of evaluation metrics (i.e., AUC, Accuracy, and so on).

###  real1_Blekhman_2010.R  ###
This file contains scripts to analyze the "[Blekhman_2010](https://pubmed.ncbi.nlm.nih.gov/20009012/)" data.

###  real2_Pickrell_2010.R  ###
This file contains scripts to analyze the "[Pickrell_2010](https://pubmed.ncbi.nlm.nih.gov/20220758/)" data.

###  real3_Riaz_2017.R  ###
This file contains scripts to analyze the "[Riaz_2017](https://pubmed.ncbi.nlm.nih.gov/29033130/)" data.

###  real4_Li_2014.R  ###
This file contains scripts to analyze the "[Li_2014](https://pubmed.ncbi.nlm.nih.gov/24441097/)" data.

