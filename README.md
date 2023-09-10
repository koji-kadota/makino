# makino

This repository contains R-codes used in Makino et al's paper (inpreparation).

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
This file contains functions to execute the four methods ([edgeR](https://bioconductor.org/packages/edgeR/), [DESeq2](https://bioconductor.org/packages/DESeq2/), [TCC](https://bioconductor.org/packages/TCC/), and [MBCdeg](https://pubmed.ncbi.nlm.nih.gov/34670485/)) in a unified manner.
