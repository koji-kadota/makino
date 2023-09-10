setwd("C:/Users/kadota/OneDrive - The University of Tokyo/ドキュメント/2023/研究/牧野論文/new/Real1_Blekhman_2010")
#########################################
### This file contains full information for analyzing the "Blekhman_2010" data.
###  the raw data file can be obtained from the following URL:
###  https://genome.cshlp.org/content/suppl/2009/12/16/gr.099226.109.DC1/suppTable1.xls
#########################################
###  Parameters  ###
N_trial = 5                #number of trials
filtering = TRUE           #filtering of zero-count genes: TRUE or FALSE
FDR = 0.1                  #FDR threshold

###  Load packages  ###
library(edgeR)
library(DESeq2)
library(TCC)
library(MBCluster.Seq)

###  Read functions  ###
source("functions.R")

###  Preparation of placeholder  ###
resall.PDEG <- NULL
resall.P1 <- NULL

###  Read and preprocess the input file  ###
in_f <- "https://genome.cshlp.org/content/suppl/2009/12/16/gr.099226.109.DC1/suppTable1.xls"
dataOri <- read.table(in_f, header=TRUE)
data18 <- cbind(
  dataOri$R1L4.HSF1 + dataOri$R4L2.HSF1, dataOri$R2L7.HSF2 + dataOri$R3L2.HSF2,
  dataOri$R8L1.HSF3 + dataOri$R8L2.HSF3, dataOri$R1L1.HSM1 + dataOri$R5L2.HSM1,
  dataOri$R2L3.HSM2 + dataOri$R4L8.HSM2, dataOri$R3L6.HSM3 + dataOri$R4L1.HSM3,
  dataOri$R1L2.PTF1 + dataOri$R4L4.PTF1, dataOri$R2L4.PTF2 + dataOri$R6L6.PTF2,
  dataOri$R3L7.PTF3 + dataOri$R5L3.PTF3, dataOri$R1L6.PTM1 + dataOri$R3L3.PTM1,
  dataOri$R2L8.PTM2 + dataOri$R4L6.PTM2, dataOri$R6L2.PTM3 + dataOri$R6L4.PTM3,
  dataOri$R1L7.RMF1 + dataOri$R5L1.RMF1, dataOri$R2L2.RMF2 + dataOri$R5L8.RMF2,
  dataOri$R3L4.RMF3 + dataOri$R4L7.RMF3, dataOri$R1L3.RMM1 + dataOri$R3L8.RMM1,
  dataOri$R2L6.RMM2 + dataOri$R5L4.RMM2, dataOri$R3L1.RMM3 + dataOri$R4L3.RMM3)
rownames(data18) <- dataOri[, 1]
colnames(data18) <- c(
  "HSF1", "HSF2", "HSF3", "HSM1", "HSM2", "HSM3", 
  "PTF1", "PTF2", "PTF3", "PTM1", "PTM2", "PTM3", 
  "RMF1", "RMF2", "RMF3", "RMM1", "RMM2", "RMM3")

HS <- data18[, 1:6]
PT <- data18[, 7:12]
RM <- data18[, 13:18]
HSF <- data18[, 1:3]
HSM <- data18[, 4:6]
PTF <- data18[, 7:9]
PTM <- data18[, 10:12]
RMF <- data18[, 13:15]
RMM <- data18[, 16:18]

####################################
###  Table 1(a): 1. HSF vs. PTF  ###
####################################
G1 <- HSF
G2 <- PTF
tmp_rownames <- paste("HSF vs. PTF (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)


####################################
###  Table 1(a): 2. HSF vs. PTM  ###
####################################
G1 <- HSF
G2 <- PTM
tmp_rownames <- paste("HSF vs. PTM (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)

####################################
###  Table 1(a): 3. HSM vs. PTF  ###
####################################
G1 <- HSM
G2 <- PTF
tmp_rownames <- paste("HSM vs. PTF (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)

####################################
###  Table 1(a): 4. HSM vs. PTM  ###
####################################
G1 <- HSM
G2 <- PTM
tmp_rownames <- paste("HSM vs. PTM (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)


####################################
###  Table 1(a): 5. HSF vs. RMF  ###
####################################
G1 <- HSF
G2 <- RMF
tmp_rownames <- paste("HSF vs. RMF (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)


####################################
###  Table 1(a): 6. HSF vs. RMM  ###
####################################
G1 <- HSF
G2 <- RMM
tmp_rownames <- paste("HSF vs. RMM (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)


####################################
###  Table 1(a): 7. HSM vs. RMF  ###
####################################
G1 <- HSM
G2 <- RMF
tmp_rownames <- paste("HSM vs. RMF (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)


####################################
###  Table 1(a): 8. HSM vs. RMM  ###
####################################
G1 <- HSM
G2 <- RMM
tmp_rownames <- paste("HSM vs. RMM (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)


####################################
###  Table 1(a): 9. PTF vs. RMF  ###
####################################
G1 <- PTF
G2 <- RMF
tmp_rownames <- paste("PTF vs. RMF (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)



#####################################
###  Table 1(a): 10. PTF vs. RMM  ###
#####################################
G1 <- PTF
G2 <- RMM
tmp_rownames <- paste("PTF vs. RMM (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)



#####################################
###  Table 1(a): 11. PTM vs. RMF  ###
#####################################
G1 <- PTM
G2 <- RMF
tmp_rownames <- paste("PTM vs. RMF (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)



#####################################
###  Table 1(a): 12. PTM vs. RMM  ###
#####################################
G1 <- PTM
G2 <- RMM
tmp_rownames <- paste("PTM vs. RMM (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)


#####################################
###  Table 1(b): 13. HSF vs. HSM  ###
#####################################
G1 <- HSF
G2 <- HSM
tmp_rownames <- paste("HSF vs. HSM (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
if(filtering == TRUE){
  data <- data[rowSums(data) != 0, ]
}
print(nrow(data))

methods <- c("edgeR", "DESeq2", "TCC",
             "MBCdeg1_K3", "MBCdeg2_K3", "MBCdeg2_K3",
             "MBCdeg1_K4", "MBCdeg2_K4", "MBCdeg2_K4")
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)


#####################################
###  Table 1(b): 14. PTF vs. PTM  ###
#####################################
G1 <- PTF
G2 <- PTM
tmp_rownames <- paste("PTF vs. PTM (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
if(filtering == TRUE){
  data <- data[rowSums(data) != 0, ]
}
print(nrow(data))

methods <- c("edgeR", "DESeq2", "TCC",
             "MBCdeg1_K3", "MBCdeg2_K3", "MBCdeg2_K3",
             "MBCdeg1_K4", "MBCdeg2_K4", "MBCdeg2_K4")
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)


#####################################
###  Table 1(b): 15. RMF vs. RMM  ###
#####################################
G1 <- RMF
G2 <- RMM
tmp_rownames <- paste("RMF vs. RMM (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
if(filtering == TRUE){
  data <- data[rowSums(data) != 0, ]
}
print(nrow(data))

methods <- c("edgeR", "DESeq2", "TCC",
             "MBCdeg1_K3", "MBCdeg2_K3", "MBCdeg2_K3",
             "MBCdeg1_K4", "MBCdeg2_K4", "MBCdeg2_K4")
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)


###################################
###  Table 1(c): 16. HS vs. PT  ###
###################################
G1 <- HS
G2 <- PT
tmp_rownames <- paste("HS vs. PT (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
if(filtering == TRUE){
  data <- data[rowSums(data) != 0, ]
}
print(nrow(data))

methods <- c("edgeR", "DESeq2", "TCC",
             "MBCdeg1_K3", "MBCdeg2_K3", "MBCdeg2_K3",
             "MBCdeg1_K4", "MBCdeg2_K4", "MBCdeg2_K4")
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)

###################################
###  Table 1(c): 17. HS vs. RM  ###
###################################
G1 <- HS
G2 <- RM
tmp_rownames <- paste("HS vs. RM (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
if(filtering == TRUE){
  data <- data[rowSums(data) != 0, ]
}
print(nrow(data))

methods <- c("edgeR", "DESeq2", "TCC",
             "MBCdeg1_K3", "MBCdeg2_K3", "MBCdeg2_K3",
             "MBCdeg1_K4", "MBCdeg2_K4", "MBCdeg2_K4")
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)

###################################
###  Table 1(c): 18. PT vs. RM  ###
###################################
G1 <- PT
G2 <- RM
tmp_rownames <- paste("PT vs. RM (trial", 1:N_trial, ")", sep="")

data <- cbind(G1, G2)
data.cl <- c(rep(1, ncol(G1)), rep(2, ncol(G2)))
if(filtering == TRUE){
  data <- data[rowSums(data) != 0, ]
}
print(nrow(data))

methods <- c("edgeR", "DESeq2", "TCC",
             "MBCdeg1_K3", "MBCdeg2_K3", "MBCdeg2_K3",
             "MBCdeg1_K4", "MBCdeg2_K4", "MBCdeg2_K4")
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

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_PDEG[i, 2] <- res$PDEG
res_P1[i, 2] <- res$P1

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_PDEG[i, 3] <- res$PDEG
res_P1[i, 3] <- res$P1
norm.factors <- res$norm.factors

K <- 3
###  4. MBCdeg1 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 4] <- res$PDEG
res_P1[i, 4] <- res$P1

###  5. MBCdeg2 (K = 3)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 5] <- res$PDEG
res_P1[i, 5] <- res$P1

###  6. MBCdeg3 (K = 3)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 6] <- res$PDEG
res_P1[i, 6] <- res$P1

K <- 4
###  7. MBCdeg1 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_PDEG[i, 7] <- res$PDEG
res_P1[i, 7] <- res$P1

###  8. MBCdeg2 (K = 4)  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_PDEG[i, 8] <- res$PDEG
res_P1[i, 8] <- res$P1

###  9. MBCdeg3 (K = 4)  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_PDEG[i, 9] <- res$PDEG
res_P1[i, 9] <- res$P1

}
rownames(res_PDEG) <- tmp_rownames
rownames(res_P1) <- tmp_rownames
resall.PDEG <- rbind(resall.PDEG, res_PDEG)
resall.P1 <- rbind(resall.P1, res_P1)


########################
###  Output1 (PDEG)  ###
########################
tmp <- cbind(rownames(resall.PDEG), resall.PDEG)
out_f <- "Additional22_PDEG.txt"
write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)

########################
###  Output2 (P1)    ###
########################
tmp <- cbind(rownames(resall.P1), resall.P1)
out_f <- "Additional22_P1.txt"
write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)
