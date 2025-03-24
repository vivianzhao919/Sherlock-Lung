
#####################################
### Batch correction for new data

library("sva") 
library("edgeR")
library("UpSetR")
library(tidyverse)



infilestr <- "/data/Sherlock_Lung/RNASeq/htseq/Sherlock1+2+TCGA-LUAD_RNAseq_htseq.genes.results.count.Combat-corrected_GeneName.txt"
tmp <- read.table(infilestr,sep="\t",header=T,row.names=1,check.names=F)
infile.sherlock <- apply(tmp,2,as.numeric)
rownames(infile.sherlock) <- rownames(tmp)
colnames(infile.sherlock) <- colnames(tmp)



infilestr <-  "/data/Sherlock_Lung/RNASeq/htseq/GIS_htseq.genes.results.count_GeneName.txt"
tmp <- read.table(infilestr,sep="\t",header=T,row.names=1,check.names=F)
infile.GIS <- apply(tmp,2,as.numeric)
rownames(infile.GIS) <- rownames(tmp)
colnames(infile.GIS) <- colnames(tmp)


genelist <- intersect(rownames(infile.sherlock),rownames(infile.GIS))

infile.sherlock <- infile.sherlock[genelist,]
infile.GIS <- infile.GIS[genelist,]

infile.merge <- cbind(infile.sherlock,infile.GIS)




infilestr <-"/data/Sherlock_Lung/RNASeq/Clinical/Sherlock_RNA_TN_GIS_info.txt"
metadata <- read.table(infilestr,sep="\t",header=T,row.names=1,check.names=F,stringsAsFactors = FALSE)
metadata <- metadata[colnames(infile.merge),]

setClass <- c(rep("Sherlock",ncol(infile.sherlock)),rep("GIS",ncol(infile.GIS)))
setClass <- factor(setClass,levels=c("Sherlock","GIS"))

TNtype <- factor(as.character(metadata[,"Type"]))

uncorrected_data <- infile.merge
pca_uncorrected_obj = prcomp(uncorrected_data)
pca_uncorrected = as.data.frame(pca_uncorrected_obj[2]$rotation)
pca_uncorrected[,"Batch"] <- setClass
pca_uncorrected[,"Type"] <-  TNtype

#######
## corrected by CombBat_Seq
groups = sapply(as.character(metadata[,"Type"]), switch, "Normal" = 1, "Tumor" = 2, USE.NAMES = F)
batches = sapply(as.character(setClass), switch, "Sherlock" = 1, "GIS"=2, USE.NAMES = F)
corrected_data = ComBat_seq(counts = as.matrix(uncorrected_data), batch = batches, group = groups)

pca_corrected_obj = prcomp(corrected_data)
pca_corrected = as.data.frame(pca_corrected_obj[2]$rotation)
pca_corrected[,"Batch"] <- setClass
pca_corrected[,"Type"] <-  TNtype

newData <- corrected_data

outfilestr <-  "/data/Sherlock_Lung/RNASeq/htseq/Sherlock1+2+2R+TCGA-LUAD+GIS_RNAseq_htseq.genes.results.count.GeneName.Combat-corrected.txt"
newData <- data.frame(newData,check.names=FALSE) %>% mutate(UID=rownames(corrected_data),.before=1)
write_tsv(newData,outfilestr)



######### normalize training and merged set with edgeR



library(edgeR)
library(tidyverse)
infilestr <- "/data/Sherlock_Lung/RNASeq/htseq/Sherlock1+2+TCGA-LUAD_RNAseq_htseq.genes.results.count.Combat-corrected_GeneName.txt"
tmp <- read.table(infilestr,sep="\t",header=T,row.names=1,check.names=F)
counts <- apply(tmp,2,as.numeric)
rownames(counts) <- rownames(tmp)
colnames(counts) <- colnames(tmp)

keep <- rowSums(counts) >= 10
counts <- counts[keep,]
d0 <- DGEList(counts)
sname <- colnames(counts)
d0 <- calcNormFactors(d0,method="TMM")


out.cpm <- cpm(d0, normalized.lib.sizes=TRUE, log=FALSE)

outfilestr <-  "/data/Sherlock_Lung/RNASeq/htseq/Sherlock1+2+TCGA-LUAD_RNAseq_htseq.genes.results.count.Combat-corrected_GeneName.edgeR_normalized_CPM.txt"
out <- data.frame(out.cpm,check.names=FALSE) %>% mutate(Sample=rownames(out.cpm),.before=1)
write_tsv(out, outfilestr)

######


infilestr <- "/data/Sherlock_Lung/RNASeq/htseq/Sherlock1+2+2R+TCGA-LUAD+GIS_RNAseq_htseq.genes.results.count.GeneName.Combat-corrected.txt"
tmp <- read.table(infilestr,sep="\t",header=T,row.names=1,check.names=F)
counts <- apply(tmp,2,as.numeric)
rownames(counts) <- rownames(tmp)
colnames(counts) <- colnames(tmp)

keep <- rowSums(counts) >= 10
counts <- counts[keep,]
d0 <- DGEList(counts)
sname <- colnames(counts)
d0 <- calcNormFactors(d0,method="TMM")


out.cpm <- cpm(d0, normalized.lib.sizes=TRUE, log=FALSE)
outfilestr <-  "/data/Sherlock_Lung/RNASeq/htseq/Sherlock1+2+2R+TCGA-LUAD+GIS_RNAseq_htseq.genes.results.count.GeneName.Combat-corrected.edgeR_normalized_CPM.txt"
out.cpm <- data.frame(out.cpm,check.names=FALSE) %>% mutate(UID=rownames(out.cpm),.before=1)
write_tsv(out.cpm,outfilestr)



