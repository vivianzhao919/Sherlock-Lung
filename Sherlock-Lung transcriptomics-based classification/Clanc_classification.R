
rm(list = ls())


library(tidyverse)
library("DescTools")
library("sake")
source("~/scripts/Clanc_classification_CV.r")
source("~/scripts/clanc.R")

####################################################################################################
# Estimate error rate of subtype prediction in the Sherlock/TCGA data set

topNgenes <- 5000
### input data: gene expression data of protein coding genes excluding mitochondria genes
#infilestr <- "/data/Sherlock_Lung/RNASeq/htseq/Sherlock_phase1+2+TCGA-LUAD_neversmoker_htseq_genes.Combat-corrected.edgeR_normalized_CPM_GeneName_coding_xMT_update.txt"
infilestr <- "/data/Sherlock_Lung/RNASeq/htseq/Sherlock1+2+TCGA-LUAD_RNAseq_htseq.genes.results.count.Combat-corrected_GeneName.edgeR_normalized_CPM_coding_xMT.txt"
tmp <- read.table(infilestr,sep="\t",header=T,row.names=1,check.names=FALSE)
htseqfile <- apply(tmp,2,as.numeric)
rownames(htseqfile) <- rownames(tmp)
colnames(htseqfile) <- colnames(tmp)

tempMedian <- apply(log2(htseqfile + 1), 1, median)
V <- as.matrix(htseqfile[tempMedian >= 1, ])
tempVar <- apply(V, 1, var)
V <- V[order(tempVar, decreasing = TRUE)[1:topNgenes], ]
logV <- log2(V+1)


infilestr <- "/data/Sherlock_Lung/RNASeq/Clinical/Sherlock_RNA_clinical_data.txt"
metadata <- read_tsv(infilestr)
metadata <- metadata%>% column_to_rownames(var="RNAseq_SampleID")

V <- V[,rownames(metadata)]

its=1234

setwd("/data/Sherlock_Lung/RNASeq/clanc/")
cv_build <- list()
options <- list()
test.err <- train.err<- {}
for(its in 1:100){
  print(its)
  set.seed(its)
  train.flag <- sample(1:nrow(metadata),round(0.7*nrow(metadata)))
  test.flag <- -train.flag
  
  train.data <- V[,train.flag]
  test.data <- V[,-train.flag]
  train.metadata <- metadata[train.flag,]
  test.metadata <- metadata[-train.flag,]
  
  file.prefix=paste("Sherlock1+2+TCGA-LUAD_CPM_coding_xMT_Top5000_rank3_20features_CV_seed",its,sep="")
  
  
  clanc.tmp <- run_ClaNC(train.data, train.metadata$CPM3, active.features=30, est.num=100, select.features=20,  file.prefix=file.prefix,skip.est=F)
  cv_build[[its]] <- clanc.tmp
  options[[its]]=list(active.features=30, est.num=100, select.features=20)
  
  test.out <- predictClanc(test.data, geneNames=rownames(test.data), fit=clanc.tmp)
  test.truth <- test.metadata$CPM3
  test.mat <- table(test.out,test.truth)
  test.err <- c(test.err,1-sum(diag(test.mat))/sum(test.mat))
  train.err <- c(train.err, clanc.tmp$overallError)
}

outfilestr<-"/data/Sherlock_Lung/RNASeq/clanc/Sherlock1+2+TCGA-LUAD_neversmoker_CPM_coding_xMT_Top5000_rank3_CV_20features.RData"

save(list=c("V","cv_build","options","train.err","test.err"),file=outfilestr)


## predict the expression subtypes in the complete Sherlock/TCGA cohort (training set)
its=1234
set.seed(its)
file.prefix=paste("Sherlock1+2+TCGA-LUAD_neversmoker_CPM_coding_xMT_Top5000_rank3_20features_CV_seed",its,sep="")

clanc_build <- run_ClaNC(V, metadata$CPM3, active.features=30, est.num=100, select.features=20,  file.prefix=file.prefix,skip.est=F)
Options=list(active.features=30, est.num=100, select.features=20)

clanc_out <- predictClanc(V, geneNames=rownames(V), fit=clanc_build)
names(clanc_out) <- colnames(V)
outfilestr<-"/data/Sherlock_Lung/RNASeq/clanc/Sherlock1+2+TCGA-LUAD_CPM_coding_xMT_Top5000_rank3_20features_final.RData"
save(list=c("V","clanc_build","Options","clanc_out"),file=outfilestr)

###
## Calculate distance to expression centroids in the complete Sherlock/TCGA cohort (training set)


n = ncol(V)
fit = clanc_build
cntrds = fit$cntrds

active.idx = match(fit$geneNames, rownames(V))
sd = fit$pooledSD
prior = fit$prior

if(any(is.na(active.idx)))
  stop("Gene names do not match those in classifier.")


distClanc_v2 <- function(data, cntrds, sd, prior) {
  vv = sd ^ 2
  pi.k = prior
  
  if(length(vv) > 1)
    dd =  drop(t(data ^ 2) %*% (1 / vv)) + drop(t(cntrds ^ 2) %*% (1 / vv)) - 2 * drop(t(data * cntrds) %*% (1 / vv)) - 2 * log(pi.k)
  else
    dd =  drop(data ^ 2 / vv)  + drop(cntrds ^ 2 / vv) - 2 * drop(data * cntrds / vv) - 2 * log(pi.k)
  
  return(dd)
}


cls = rep(NA, n)
dist.mat <- {}
for(i in 1:n) {
  dd = distClanc_v2(data = V[active.idx, i], cntrds = cntrds, sd = sd, prior = prior)
  dist.mat <- rbind(dist.mat,dd)
  cls[i] = match(min(dd), dd)
}

fit$className <- c("steady","proliferative","chaotic")
colnames(dist.mat) <- fit$className
out <- data.frame(cls,dist.mat,check.names=FALSE) %>% mutate(Samples=colnames(V),.before=1)
outfilestr <-  "/data/Sherlock_Lung/RNASeq/clanc/Sherlock1+2+TCGA-LUAD_CPM_coding_xMT_Top5000_rank3_20features_prediction_dist.txt"
write_tsv(out,outfilestr)



#######
## Apply subtype predictor to the GIS data set (testing set)




infilestr <- "/data/Sherlock_Lung/RNASeq/htseq/Sherlock1+2+2R+TCGA-LUAD+GIS_RNAseq_htseq.genes.results.count.Combat-corrected.edgeR_normalized_CPM_GeneName_coding_xMT.txt"

tmp <- read.table(infilestr,sep="\t",header=T,row.names=1,check.names=FALSE)
newFile <- apply(tmp,2,as.numeric)
rownames(newFile) <- rownames(tmp)
colnames(newFile) <- colnames(tmp)

infilestr <-"/data/Sherlock_Lung/RNASeq/Clinical/Sherlock_RNA_TN_GIS_info.txt"
metadata <- read_tsv(infilestr)
newV <- newFile[,metadata$SAMPLE_ID]

infilestr <- "/data/Sherlock_Lung/RNASeq/GeneExpAnalysis/NMF/clanc/Sherlock1+2+TCGA-LUAD_CPM_coding_xMT_Top5000_rank3_20features_final.RData"
load(infilestr)
clanc_out <- predictClanc(newV, geneNames=rownames(newV), fit=clanc_build) 
names(clanc_out) <- colnames(newV)
outfilestr<- "/data/Sherlock_Lung/RNASeq/GeneExpAnalysis/NMF/clanc/GIS_CPM_coding_xMT_Top5000_rank3_20features_final.RData"
save(list=c("V","clanc_build","Options","clanc_out"),file=outfilestr)

###
## Calculate distance to expression centroids in the GIS data set (testing set)


n = ncol(newV)
fit = clanc_build
cntrds = fit$cntrds

active.idx = match(fit$geneNames, rownames(newV))
sd = fit$pooledSD
prior = fit$prior

if(any(is.na(active.idx)))
  stop("Gene names do not match those in classifier.")



distClanc_v2 <- function(data, cntrds, sd, prior) {
  vv = sd ^ 2
  pi.k = prior
  
  if(length(vv) > 1)
    dd =  drop(t(data ^ 2) %*% (1 / vv)) + drop(t(cntrds ^ 2) %*% (1 / vv)) - 2 * drop(t(data * cntrds) %*% (1 / vv)) - 2 * log(pi.k)
  else
    dd =  drop(data ^ 2 / vv)  + drop(cntrds ^ 2 / vv) - 2 * drop(data * cntrds / vv) - 2 * log(pi.k)
  
  return(dd)
}


cls = rep(NA, n)
dist.mat <- {}
for(i in 1:n) {
  dd = distClanc_v2(data = newV[active.idx, i], cntrds = cntrds, sd = sd, prior = prior)
  dist.mat <- rbind(dist.mat,dd)
  cls[i] = match(min(dd), dd)
}

fit$className <- c("steady","proliferative","chaotic")
colnames(dist.mat) <- fit$className
out <- data.frame(cls,dist.mat,check.names=FALSE) %>% mutate(Samples=colnames(newV),.before=1)
outfilestr <-  "/data/Sherlock_Lung/RNASeq/GeneExpAnalysis/NMF/clanc/GIS_CPM_coding_xMT_Top5000_rank3_20features_prediction_dist.txt"
write_tsv(out,outfilestr)

