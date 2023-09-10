#########################################
### This file contains full R scripts for analyzing the "Pickrell_2010" data (i.e., Additional23a).
### The count data can be obtained from the recount2 website:
###  https://jhubiostatistics.shinyapps.io/recount/
### setwd("C:/Users/kadota/OneDrive - The University of Tokyo/ドキュメント/2023/研究/牧野論文/program")
#########################################
###  Parameters  ###
accession = "SRP001540"    #SRA ID
N_trial = 5                #number of trials
n1 = 40                    #number of replicates for group 1
n2 = 29                    #number of replicates for group 2
filtering = TRUE           #filtering of zero-count genes: TRUE or FALSE
FDR = 0.1                  #FDR threshold (fixed)

###  Load packages  ###
library(edgeR)
library(DESeq2)
library(TCC)
library(MBCluster.Seq)
library(recount)

###  Read functions  ###
source("functions.R")

###  Preparation of placeholder  ###
res.list <- list()
resall.PDEG <- NULL
resall.P1 <- NULL

###  Read and preprocess the input file  ###
download_study(accession, type="rse-gene", download=T) 
in_f <- paste(accession, "/rse_gene.Rdata", sep="")
load(in_f)
rse <- rse_gene
rse <- scale_counts(rse)
x <- assays(rse)$counts
x <- as.data.frame(x)
data <- cbind(x$SRR031822, x$SRR031953 + x$SRR031873,             #Female1-2 
              x$SRR031952 + x$SRR031871, x$SRR031868,             #Female3-4
              x$SRR031819, x$SRR031897 + x$SRR031857,             #Female5-6
              x$SRR031823, x$SRR031959, x$SRR031955,              #Female7-9
              x$SRR031954, x$SRR031956, x$SRR031838,              #Female10-12
              x$SRR031918, x$SRR031817, x$SRR031949 + x$SRR031852,#Female13-15
              x$SRR031841, x$SRR031865, x$SRR031896,              #Female16-18
              x$SRR031853, x$SRR031820, x$SRR031874,              #Female19-21
              x$SRR031895, x$SRR031870, x$SRR031839,              #Female22-24
              x$SRR031958, x$SRR031867, x$SRR031848,              #Female25-27
              x$SRR031847, x$SRR031818, x$SRR031919,              #Female28-30
              x$SRR031866, x$SRR031849, x$SRR031877,              #Female31-33
              x$SRR031814, x$SRR031914, x$SRR031812,              #Female34-36
              x$SRR031842, x$SRR031843, x$SRR031860, x$SRR031837, #Female37-40
              x$SRR031917, x$SRR031821 + x$SRR031898,             #Male1-2
              x$SRR031950 + x$SRR031850, x$SRR031876 + x$SRR031862,#Male3-4
              x$SRR031875, x$SRR031915, x$SRR031878 + x$SRR031863,#Male5-7
              x$SRR031869, x$SRR031864, x$SRR031845,              #Male8-10
              x$SRR031951 + x$SRR031851, x$SRR031846,             #Male11-12
              x$SRR031916, x$SRR031844, x$SRR031813,              #Male13-15
              x$SRR031894, x$SRR031854, x$SRR031858,              #Male16-18
              x$SRR031859, x$SRR031872, x$SRR031816,              #Male19-21
              x$SRR031815, x$SRR031920 + x$SRR031899,             #Male22-23
              x$SRR031957 + x$SRR031855, x$SRR031840,             #Male24-25
              x$SRR031948, x$SRR031893, x$SRR031811, x$SRR031861) #Male26-29

data.cl <- c(rep(1, n1), rep(2, n2))
if(filtering == TRUE){
  data <- data[rowSums(data) != 0, ]
}
print(nrow(data))

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
out_f <- "real2_Additional23a_PDEG.txt"
write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)

###  Output2 (P1)  ###
tmp <- cbind(rownames(resall.P1), resall.P1)
out_f <- "real2_Additional23a_P1.txt"
write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)

###  Output3 (all)  ###
out_f <- "real3_Additional23a.obj"
saveRDS(res.list, out_f)

