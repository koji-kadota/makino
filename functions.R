###  MBCdeg  ###
MBCdeg <- function(counts, treatment, K, normalizer = NULL, geneid = NULL) {
  
  message("1. Preprocessing data...")
  mbc <- RNASeq.Data(counts, Normalizer = normalizer,
                     Treatment = treatment, GeneID = geneid)
  message("2. Initializing centers...")
  c0 <- KmeansPlus.RNASeq(data = mbc, nK = K,
                          model = "nbinom", print.steps = F)
  message("3. Clustering data...")
  capture.output({
    cls <- Cluster.RNASeq(data = mbc, model = "nbinom",
                          centers = c0$centers, method = "EM")
  })
  
  ###  Identify the k value that corresponds to the non-DEG cluster  ###
  cluster <- cls$cluster
  centers <- cls$centers
  L2norm <- sqrt(rowSums(abs(centers)^2))
  #k <- which.min(L2norm)
  obj <- is.element(1:nrow(centers), names(table(cluster)))
  k <- names(table(cluster))[which.min(L2norm[obj])]
  k <- as.integer(k)

  ###  Obtaining other results such as q-values  ###
  pp <- cls$probability
  pp_nonDEG <- cls$probability[, k]
  estimatedDEG <- rep(0, length(cluster))
  estimatedDEG[cluster != k] <- 1
  obj <- (estimatedDEG == 1)

  serial <- 1:length(cluster)
  res <- data.frame(pp_nonDEG = pp_nonDEG,
                    serial = serial)
  res.sort <- res[order(pp_nonDEG), ]
  cumsum <- cumsum(res.sort$pp_nonDEG)
  qval <- cumsum/serial
  res.sort <- cbind(res.sort, qval)
  q.value <- res.sort[order(res.sort$serial), ]$qval
  PDEG <- 100*sum(obj)/length(obj)
  pattern <- rep("DEG2", nrow(centers))
  pattern[centers[, 1] > 0] <- "DEG1"
  pattern[k] <- "nonDEG"
  for(j in 1:length(pattern)){
    cluster[cluster == j] <- pattern[j]
  }
  P1 <- 100*sum(cluster == "DEG1")/sum(obj)

  rtn <- list()
  #rtn$centers <- cls$centers
  rtn$centers <- centers      # changed
  rtn$pattern <- pattern      # added
  rtn$PP <- cls$probability
  rtn$PP_nonDEG <- cls$probability[, k]
  rtn$ranking <- rank(cls$probability[, k])
  rtn$cluster <- cls$cluster
  rtn$cluster_pat <- cluster  # added
  rtn$estimatedDEG <- estimatedDEG
  rtn$q.value <- q.value
  rtn$PDEG <- PDEG
  rtn$P1 <- P1
  
  class(rtn) <- "mbc"
  return(rtn)
}


###  edgeR  ###
my.edger <- function(counts, group, FDR){
  dge <- edgeR::DGEList(counts = counts, group = group)
  dge <- edgeR::calcNormFactors(dge)
  design <- model.matrix(~group)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design)
  qlf <- edgeR::glmQLFTest(fit, coef = 2,)
  table <- topTags(qlf, n = nrow(counts), sort.by = "none", p.value = 1)@.Data[[1]]

  p.value <- table$PValue
  p.value[is.na(p.value)] <- 1
  q.value <- table$FDR
  q.value[is.na(q.value)] <- 1
  ranking <- rank(p.value)
  estimatedDEG <- rep(0, length(q.value))
  estimatedDEG[q.value < FDR] <- 1
  logratio <- table$logFC
  norm.factors <- dge$samples$norm.factors
  ef.libsizes <- colSums(counts)*norm.factors
  size.factors <- ef.libsizes/mean(ef.libsizes)
  PDEG <- 100*sum(estimatedDEG)/length(estimatedDEG)
  obj <- (q.value < FDR)
  P1 <- 100*sum(logratio[obj] < 0)/sum(obj)

  rtn <- list()
  rtn$ranking <- ranking
  rtn$p.value <- p.value
  rtn$q.value <- q.value
  rtn$estimatedDEG <- estimatedDEG
  rtn$logratio <- logratio
  rtn$norm.factors <- norm.factors
  rtn$size.factors <- size.factors
  rtn$PDEG <- PDEG
  rtn$P1 <- P1

  return(rtn)
}


###  DESeq2  ###
my.deseq2 <- function(counts, group, FDR){
  colData <- data.frame(condition=as.factor(group))
  dds <- DESeqDataSetFromMatrix(countData=counts,
                                colData=colData, design=~condition)
  d <- DESeq(dds)
  tmp <- results(d)

  p.value <- tmp$pvalue
  p.value[is.na(p.value)] <- 1
  q.value <- tmp$padj
  q.value[is.na(q.value)] <- 1
  ranking <- rank(p.value)
  estimatedDEG <- rep(0, length(q.value))
  estimatedDEG[q.value < FDR] <- 1
  logratio <- tmp$log2FoldChange
  size.factors <- sizeFactors(d)
  PDEG <- 100*sum(estimatedDEG)/length(estimatedDEG)
  obj <- (q.value < FDR)
  P1 <- 100*sum(logratio[obj] < 0)/sum(obj)

  rtn <- list()
  rtn$ranking <- ranking
  rtn$p.value <- p.value
  rtn$q.value <- q.value
  rtn$estimatedDEG <- estimatedDEG
  rtn$logratio <- logratio
  rtn$size.factors <- size.factors
  rtn$PDEG <- PDEG
  rtn$P1 <- P1

  return(rtn)
}


###  TCC  ###
my.tcc <- function(counts, group, FDR){
  tcc <- new('TCC', counts, group)
  tcc <- calcNormFactors(tcc, norm.method = "tmm",
                         test.method = "edger", iteration = 3)
  tcc <- estimateDE(tcc, test.method = "edger")
  #tcc <- calcNormFactors(tcc)
  #tcc <- estimateDE(tcc)
  df <- getResult(tcc, sort=FALSE)
  
  p.value <- df$p.value
  p.value[is.na(p.value)] <- 1
  q.value <- df$q.value
  q.value[is.na(q.value)] <- 1
  ranking <- rank(p.value)
  estimatedDEG <- rep(0, length(q.value))
  estimatedDEG[q.value < FDR] <- 1
  logratio <- df$m.value
  norm.factors <- tcc$norm.factors
  ef.libsizes <- colSums(counts)*norm.factors
  size.factors <- ef.libsizes/mean(ef.libsizes)
  PDEG <- 100*sum(estimatedDEG)/length(estimatedDEG)
  obj <- (q.value < FDR)
  P1 <- 100*sum(logratio[obj] < 0)/sum(obj)

  rtn <- list()
  rtn$ranking <- ranking
  rtn$p.value <- p.value
  rtn$q.value <- q.value
  rtn$estimatedDEG <- estimatedDEG
  rtn$logratio <- logratio
  rtn$norm.factors <- norm.factors
  rtn$size.factors <- size.factors
  rtn$PDEG <- PDEG
  rtn$P1 <- P1

  return(rtn)
}

###  Accuracy  ###
accuracy <- function(truth, data){
  res <- sum(truth == data)/length(data)
  return(res)
}

###  FDR control  ###
actualfdr <- function(truth, q.value){
  nominalFDR <- seq(0, 1, 0.001)
  actualFDR <- NULL
  for(t in nominalFDR){
    DEGnum.t <- sum(q.value < t)
    falseDEG.t <- sum(truth[(q.value < t)] == 0)
    if(DEGnum.t == 0){
      actualFDR.t <- 0
    }else{
      actualFDR.t <- falseDEG.t/DEGnum.t
    }
    actualFDR <- rbind(actualFDR, actualFDR.t)
  }
  return(actualFDR)
}

##### Generic functions #####

#' print mbc object
#'
#' @param obj A mbc object
print.mbc <- function(obj) {
  cat("This is a mbc class object.\n")
  cat("It contains some values below.\n")
  cat("  centers: centers of each cluster.\n")
  cat("  PP: posterior probability.\n")
  cat("  PP_nonDEG: posterior probability of nonDEG cluster.\n")
  cat("  ranking: ranking of genes by PP_nonDEG.\n")
  cat("  cluster: cluster to which each gene is assigned.\n")
}
