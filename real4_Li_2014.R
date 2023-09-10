setwd("C:/Users/kadota/OneDrive - The University of Tokyo/ドキュメント/2023/研究/牧野論文/program")
#########################################
### This file contains full R scripts for analyzing the "Li_2014" data (i.e., Additional23c).
### The count data can be obtained from the recount2 website:
###  https://jhubiostatistics.shinyapps.io/recount/
#########################################
###  Parameters  ###
accession <- "SRP035988"    #SRA ID
N_trial = 5                #number of trials
n1 = 92                    #number of replicates for group 1
n2 = 81                    #number of replicates for group 2
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
data <- cbind(x$SRR1146077, x$SRR1146078, x$SRR1146079,#lesional01-03 
              x$SRR1146080, x$SRR1146081, x$SRR1146082,#lesional04-06
              x$SRR1146083 + x$SRR1146084, x$SRR1146087, x$SRR1146089,#lesional07-09
              x$SRR1146090, x$SRR1146091, x$SRR1146092,#lesional10-12
              x$SRR1146093, x$SRR1146094, x$SRR1146095,#lesional13-15
              x$SRR1146097, x$SRR1146098 + x$SRR1146099, x$SRR1146100,#lesional16-18
              x$SRR1146102, x$SRR1146103, x$SRR1146104,#lesional19-21
              x$SRR1146105, x$SRR1146110, x$SRR1146112,#lesional22-24
              x$SRR1146119, x$SRR1146123, x$SRR1146124,#lesional25-27
              x$SRR1146128, x$SRR1146129, x$SRR1146132,#lesional28-30
              x$SRR1146133, x$SRR1146134, x$SRR1146135,#lesional31-33
              x$SRR1146136, x$SRR1146140, x$SRR1146141,#lesional34-36
              x$SRR1146147, x$SRR1146148, x$SRR1146149,#lesional37-39
              x$SRR1146150, x$SRR1146151, x$SRR1146152,#lesional40-42
              x$SRR1146154, x$SRR1146155, x$SRR1146156,#lesional43-45
              x$SRR1146157, x$SRR1146158, x$SRR1146159,#lesional46-48
              x$SRR1146160, x$SRR1146161, x$SRR1146162,#lesional49-51
              x$SRR1146163, x$SRR1146164, x$SRR1146165,#lesional52-54
              x$SRR1146168, x$SRR1146199, x$SRR1146202,#lesional55-57
              x$SRR1146203, x$SRR1146204, x$SRR1146205,#lesional58-60
              x$SRR1146206, x$SRR1146208, x$SRR1146209,#lesional61-63
              x$SRR1146210, x$SRR1146211, x$SRR1146212,#lesional64-66
              x$SRR1146214, x$SRR1146215, x$SRR1146216 + x$SRR1146217,#lesional67-69
              x$SRR1146218, x$SRR1146219, x$SRR1146220,#lesional70-72
              x$SRR1146221, x$SRR1146222, x$SRR1146223,#lesional73-75
              x$SRR1146224, x$SRR1146225, x$SRR1146226,#lesional76-78
              x$SRR1146227, x$SRR1146228, x$SRR1146229,#lesional79-81
              x$SRR1146233, x$SRR1146234, x$SRR1146235,#lesional82-84
              x$SRR1146237, x$SRR1146239, x$SRR1146240,#lesional85-87
              x$SRR1146241, x$SRR1146242, x$SRR1146252,#lesional88-90
              x$SRR1146253, x$SRR1146254,              #lesional91-92
              x$SRR1146076, x$SRR1146085 + x$SRR1146086, x$SRR1146088,#normal01-03
              x$SRR1146096, x$SRR1146101, x$SRR1146106,#normal04-06
              x$SRR1146107, x$SRR1146108, x$SRR1146109,#normal07-09
              x$SRR1146111, x$SRR1146113, x$SRR1146114,#normal10-12
              x$SRR1146115, x$SRR1146116, x$SRR1146117,#normal13-15
              x$SRR1146118, x$SRR1146120, x$SRR1146121,#normal16-18
              x$SRR1146122, x$SRR1146125, x$SRR1146126,#normal19-21
              x$SRR1146127, x$SRR1146130, x$SRR1146131,#normal22-24
              x$SRR1146137, x$SRR1146138, x$SRR1146139,#normal25-27
              x$SRR1146142, x$SRR1146143, x$SRR1146144,#normal28-30
              x$SRR1146145, x$SRR1146146, x$SRR1146153,#normal31-33
              x$SRR1146166, x$SRR1146167, x$SRR1146169,#normal34-36
              x$SRR1146170, x$SRR1146171, x$SRR1146172,#normal37-39
              x$SRR1146173, x$SRR1146174, x$SRR1146175,#normal40-42
              x$SRR1146176, x$SRR1146177, x$SRR1146178,#normal43-45
              x$SRR1146179, x$SRR1146180, x$SRR1146181,#normal46-48
              x$SRR1146182, x$SRR1146183, x$SRR1146184,#normal49-51
              x$SRR1146185, x$SRR1146186, x$SRR1146187,#normal52-54
              x$SRR1146188, x$SRR1146189, x$SRR1146190,#normal55-57
              x$SRR1146191, x$SRR1146192, x$SRR1146193,#normal58-60
              x$SRR1146194, x$SRR1146195, x$SRR1146196,#normal61-63
              x$SRR1146197, x$SRR1146198, x$SRR1146200,#normal64-66
              x$SRR1146201, x$SRR1146207, x$SRR1146213,#normal67-69
              x$SRR1146230, x$SRR1146231 + x$SRR1146232, x$SRR1146236,#normal70-72
              x$SRR1146238, x$SRR1146243, x$SRR1146244,#normal73-75
              x$SRR1146245, x$SRR1146246, x$SRR1146247,#normal76-78
              x$SRR1146248, x$SRR1146249, x$SRR1146250)#normal79-81

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
out_f <- "real4_Additional23c_PDEG.txt"
write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)

###  Output2 (P1)  ###
tmp <- cbind(rownames(resall.P1), resall.P1)
out_f <- "real4_Additional23c_P1.txt"
write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)

###  Output3 (all)  ###
out_f <- "real4_Additional23c.obj"
saveRDS(res.list, out_f)

