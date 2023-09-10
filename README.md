# makino

This repository contains R-codes used in Makino et al's paper (inpreparation).

## Prerequisite 
To execute the codes, [edgeR](https://bioconductor.org/packages/edgeR/), [DESeq2](https://bioconductor.org/packages/DESeq2/), [TCC](https://bioconductor.org/packages/TCC/), [ROC](https://bioconductor.org/packages/ROC/), [compcodeR](https://bioconductor.org/packages/compcodeR/), [PROPER](https://bioconductor.org/packages/PROPER/), [MBCluster.Seq](https://CRAN.R-project.org/package=MBCluster.Seq), [recount](https://bioconductor.org/packages/recount/), and [stringr](https://bioconductor.org/packages/stringr/) are needed:  
```r
if (!requireNamespace("BiocManager", quietly=T))
    install.packages("BiocManager")
BiocManager::install("compcodeR", update=F)
BiocManager::install("PROPER", update=F)
BiocManager::install("recount", update=F)
BiocManager::install("ggsci", update=F)
BiocManager::install("TCC", update=F)
BiocManager::install("tidyverse", update=F)
install.packages("https://cran.r-project.org/src/contrib/Archive/MBCluster.Seq/MBCluster.Seq_1.0.tar.gz")
BiocManager::install("DESeq2", update=F)
```
