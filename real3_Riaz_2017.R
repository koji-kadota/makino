setwd("C:/Users/kadota/OneDrive - The University of Tokyo/ドキュメント/2023/研究/牧野論文/new/Real3_Riaz_2017")
#########################################
### This file contains full R scripts for analyzing the "Riaz_2017" data (i.e., Additional23b).
### The count data can be obtained from the GEO website:
### https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061
### Specifically, as follows:
### https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz
#########################################
###  Parameters  ###
in_f <- "GSE91061_BMS038109Sample.hg19KnownGene.raw.csv" #input filename
N_trial = 5                #number of trials
n1 = 51                    #number of replicates for group 1
n2 = 58                    #number of replicates for group 2
filtering = TRUE           #filtering of zero-count genes: TRUE or FALSE
FDR = 0.1                  #FDR threshold (fixed)

###  Load packages  ###
library(edgeR)
library(DESeq2)
library(TCC)
library(MBCluster.Seq)
library(recount)
library(stringr)

###  Read functions  ###
source("functions.R")

###  Preparation of placeholder  ###
res.list <- list()
resall.PDEG <- NULL
resall.P1 <- NULL

###  Read and preprocess the input file  ###
data.all <- read.csv(in_f, row.names = 1)
Col_pre <- str_detect(colnames(data.all), pattern="Pre") 
Col_on <- str_detect(colnames(data.all), pattern="On")　 
data_pre <- data.all[,Col_pre] 
data_on <- data.all[,Col_on]   
data <- cbind(data_pre,data_on)

data.cl <- c(rep(1, n1), rep(2, n2))
if(filtering == TRUE){
  data <- data[rowSums(data) != 0, ]
}
print(dim(data))

methods <- c("edgeR", "DESeq2", "TCC",
             "MBCdeg1_K3", "MBCdeg2_K3", "MBCdeg3_K3",
             "MBCdeg1_K4", "MBCdeg2_K4", "MBCdeg3_K4")
res_PDEG <- matrix(0, nrow = N_trial, ncol = length(methods))
res_P1 <- matrix(0, nrow = N_trial, ncol = length(methods))
colnames(res_PDEG) <- colnames(res_P1) <- methods

###  Loop ###
for (i in 1:N_trial){
print(i)
print(date())
set.seed(i)

###  1. edgeR  ###
res <- my.edger(data, data.cl, FDR = FDR)
res_PDEG[i, 1] <- res$PDEG
res_P1[i, 1] <- res$P1
res.list$edger[[i]] <- res

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1
res.list$deseq2[[i]] <- res

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors
res.list$tcc[[i]] <- res

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1
res.list$mbcdeg1k3[[i]] <- res

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1
res.list$mbcdeg2k3[[i]] <- res

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1
res.list$mbcdeg3k3[[i]] <- res

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1
res.list$mbcdeg1k4[[i]] <- res

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1
res.list$mbcdeg2k4[[i]] <- res

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1
res.list$mbcdeg3k4[[i]] <- res

}

tmp_rownames <- paste("trial", 1:N_trial, sep="")
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)


###  Output1 (PDEG)  ###
tmp <- cbind(rownames(resall.PDEG), resall.PDEG)
out_f <- "real3_Additional23b_PDEG.txt"
write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)

###  Output2 (P1)  ###
tmp <- cbind(rownames(resall.P1), resall.P1)
out_f <- "real3_Additional23b_P1.txt"
write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)

###  Output3 (all)  ###
out_f <- "real3_Additional23b.obj"
saveRDS(res.list, out_f)

