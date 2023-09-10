#########################################
### This file can be used to obtain the results of simulated data (TCC).
#########################################
###  Parameters (variable)  ###
N_trial = 20                #number of trials
n1 = 3                      #number of replicates for group 1: n1 = {3, 6, 9, 12}
PDEG = 0.05                 #percentage of DEG: PDEG = {0.05, 0.25, 0.45, 0.55, 0.65, 0.75}
P1 = 0.5                    #percentage of up-regulated DEGs in group 1: P1 = {0.5, 0.7, 0.9, 1.0]
FC = 4.0                    #degree of fold change: FC = {1.5, 2.2, 4.0}
K = 3                       #preselected number of clusters: K = {2, 3, 4, 5, 10, 15, 20}
filtering = FALSE           #filtering of zero-count genes: TRUE or FALSE

###  Parameters (fixed)  ###
name = "tcc"                #package name for generating simulated data: name = {"tcc", "com", "pro"} (fixed)
n2 = n1                     #number of replicates for group 2 (fixed)
P2 = (1 - P1)               #percentage of up-regulated DEGs in group 2 (fixed)
G = 10000                                  #number of genes (fixed)
FDR = 0.1                                  #FDR threshold (fixed)
fdr_thres = seq(0, 1, 0.001)               #bins of FDR thresholds (fixed)
methods = c("edgeR", "DESeq2", "TCC",      #order of the methods (fixed)
             paste("MBCdeg1_K", K, sep=""),#order of the methods (fixed)
             paste("MBCdeg2_K", K, sep=""),#order of the methods (fixed)
             paste("MBCdeg3_K", K, sep=""))#order of the methods (fixed)

###  Load packages  ###
library(edgeR)
library(DESeq2)
library(TCC)
library(MBCluster.Seq)
library(ROC)

###  Read functions  ###
source("functions.R")

###  Preparation of placeholder  ###
res.list <- list()
res_auc <- matrix(0, nrow = N_trial, ncol = length(methods))
res_acc <- matrix(0, nrow = N_trial, ncol = length(methods))
res_maxq <- matrix(0, nrow = N_trial, ncol = length(methods))
res_pdeg <- matrix(0, nrow = N_trial, ncol = length(methods))
res_p1 <- matrix(0, nrow = N_trial, ncol = length(methods))
res_fdr <- NULL
fdr_hoge1 <- fdr_hoge2 <- fdr_hoge3 <- fdr_hoge4 <- fdr_hoge5 <- fdr_hoge6 <- NULL

###################
###  Main loop  ###
###################
for(i in 1:N_trial){
print(i)

###  Generation of simulated data  ###
set.seed(i)
if(FC == 2.2){
  fc.matrix <- makeFCMatrix(Ngene=G,
                 PDEG = PDEG,
                 DEG.assign = c(P1, P2),
                 replicates = c(n1, n2))
  tccData <- simulateReadCounts(Ngene = G,
             PDEG = PDEG,
             DEG.assign = c(P1, P2),
             replicates = c(n1, n2),
             fc.matrix = fc.matrix)
}else{
  tccData <- simulateReadCounts(Ngene = G,
             PDEG = PDEG,
             DEG.assign = c(P1, P2),
             replicates = c(n1, n2),
             DEG.foldchange = c(FC, FC))
}

truth <- as.numeric(tccData$simulation$trueDEG != 0)
data <- tccData$count
data.cl <- tccData$group$group
if(filtering == TRUE){
  truth <- truth[rowSums(data) != 0]
  data <- data[rowSums(data) != 0, ]
}
data <- na.omit(data)
print(dim(data))
counts <- data
group <- data.cl

###  1. edgeR  ###
res <- my.edger(data, data.cl, FDR = FDR)
res_auc[i, 1] <- AUC(rocdemo.sca(truth = truth, data = -res$ranking))
res_acc[i, 1] <- accuracy(truth, res$estimatedDEG)
res_maxq[i, 1] <- max(res$q.value)
res_pdeg[i, 1] <- res$PDEG
res_p1[i, 1] <- res$P1
fdr_hoge1 <- cbind(fdr_hoge1, actualfdr(truth, res$q.value))
res.list$edgeR[[i]] <- res

###  2. DESeq2  ###
res <- my.deseq2(data, data.cl, FDR = FDR)
res_auc[i, 2] <- AUC(rocdemo.sca(truth = truth, data = -res$ranking))
res_acc[i, 2] <- accuracy(truth, res$estimatedDEG)
res_maxq[i, 2] <- max(res$q.value)
res_pdeg[i, 2] <- res$PDEG
res_p1[i, 2] <- res$P1
fdr_hoge2 <- cbind(fdr_hoge2, actualfdr(truth, res$q.value))
res.list$DESeq2[[i]] <- res

###  3. TCC  ###
res <- my.tcc(data, data.cl, FDR = FDR)
res_auc[i, 3] <- AUC(rocdemo.sca(truth = truth, data = -res$ranking))
res_acc[i, 3] <- accuracy(truth, res$estimatedDEG)
res_maxq[i, 3] <- max(res$q.value)
res_pdeg[i, 3] <- res$PDEG
res_p1[i, 3] <- res$P1
fdr_hoge3 <- cbind(fdr_hoge3, actualfdr(truth, res$q.value))
res.list$TCC[[i]] <- res
norm.factors <- res$norm.factors

###  4. MBCdeg1  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = NULL, geneid = rownames(data))
res_auc[i, 4] <- AUC(rocdemo.sca(truth = truth, data = -res$ranking))
res_acc[i, 4] <- accuracy(truth, res$estimatedDEG)
res_maxq[i, 4] <- max(res$q.value)
res_pdeg[i, 4] <- res$PDEG
res_p1[i, 4] <- res$P1
fdr_hoge4 <- cbind(fdr_hoge4, actualfdr(truth, res$q.value))
res.list$MBCdeg1[[i]] <- res

###  5. MBCdeg2  ###
res <- MBCdeg(data, data.cl, K = K, normalizer = log(norm.factors), geneid = rownames(data))
res_auc[i, 5] <- AUC(rocdemo.sca(truth = truth, data = -res$ranking))
res_acc[i, 5] <- accuracy(truth, res$estimatedDEG)
res_maxq[i, 5] <- max(res$q.value)
res_pdeg[i, 5] <- res$PDEG
res_p1[i, 5] <- res$P1
fdr_hoge5 <- cbind(fdr_hoge5, actualfdr(truth, res$q.value))
res.list$MBCdeg2[[i]] <- res

###  6. MBCdeg3  ###
nf <- 1000000/colSums(data)
res <- MBCdeg(data, data.cl, K = K, normalizer = log(nf), geneid = rownames(data))
res_auc[i, 6] <- AUC(rocdemo.sca(truth = truth, data = -res$ranking))
res_acc[i, 6] <- accuracy(truth, res$estimatedDEG)
res_maxq[i, 6] <- max(res$q.value)
res_pdeg[i, 6] <- res$PDEG
res_p1[i, 6] <- res$P1
fdr_hoge6 <- cbind(fdr_hoge6, actualfdr(truth, res$q.value))
res.list$MBCdeg3[[i]] <- res

}

###  Calculating the average actual FDR values across trials  ###
res_fdr <- cbind(res_fdr, rowMeans(fdr_hoge1))
res_fdr <- cbind(res_fdr, rowMeans(fdr_hoge2))
res_fdr <- cbind(res_fdr, rowMeans(fdr_hoge3))
res_fdr <- cbind(res_fdr, rowMeans(fdr_hoge4))
res_fdr <- cbind(res_fdr, rowMeans(fdr_hoge5))
res_fdr <- cbind(res_fdr, rowMeans(fdr_hoge6))
colnames(res_fdr) <- methods


###  Output1 (AUC)  ###
colnames(res_auc) <- methods
out_f <- paste("res_auc_", name, "_", PDEG, "_", sprintf("%.1f", P1), "_n",
               sprintf("%02d", n1), "_K", K, ".txt", sep="")
write.table(res_auc, out_f, sep="\t", append=F, quote=F, row.names=F)

###  Output2 (Accuracy)  ###
colnames(res_acc) <- methods
out_f <- paste("res_acc_", name, "_", PDEG, "_", sprintf("%.1f", P1), "_n",
               sprintf("%02d", n1), "_K", K, ".txt", sep="")
write.table(res_acc, out_f, sep="\t", append=F, quote=F, row.names=F)

###  Output3 (FDR control)  ###
tmp <- cbind(fdr_thres, res_fdr)
out_f <- paste("res_fdr_", name, "_", PDEG, "_", sprintf("%.1f", P1), "_n",
               sprintf("%02d", n1), "_K", K, ".txt", sep="")
write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)

###  Output4 (Full results)  ###
out_f <- paste("res_all_", name, "_", PDEG, "_", sprintf("%.1f", P1), "_n",
               sprintf("%02d", n1), "_K", K, ".obj", sep="")
saveRDS(res.list, out_f)

###  Output5 (max(q) values)  ###
colnames(res_maxq) <- methods
out_f <- paste("res_maxq_", name, "_", PDEG, "_", sprintf("%.1f", P1), "_n",
               sprintf("%02d", n1), "_K", K, ".txt", sep="")
write.table(res_maxq, out_f, sep="\t", append=F, quote=F, row.names=F)

###  Output6 (PDEG)  ###
colnames(res_pdeg) <- methods
out_f <- paste("res_pdeg_", name, "_", PDEG, "_", sprintf("%.1f", P1), "_n",
               sprintf("%02d", n1), "_K", K, ".txt", sep="")
write.table(res_pdeg, out_f, sep="\t", append=F, quote=F, row.names=F)

###  Output7 (P1)  ###
colnames(res_p1) <- methods
out_f <- paste("res_p1_", name, "_", PDEG, "_", sprintf("%.1f", P1), "_n",
               sprintf("%02d", n1), "_K", K, ".txt", sep="")
write.table(res_p1, out_f, sep="\t", append=F, quote=F, row.names=F)
