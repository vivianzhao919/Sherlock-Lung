### Cell type deconvolution using CAMTHC

library(CAMTHC)
library(GSA)
library(tidyverse)

infilestr <-  "/data/Sherlock_Lung/RNASeq/htseq/Sherlock_phase1+2+TCGA-LUAD_htseq_genes.Combat-corrected.TPM_GeneName.txt"
tmp <- read.table(infilestr,sep="\t",header=T,row.names=1,check.names=FALSE)
htseqfile <- apply(tmp,2,as.numeric)
rownames(htseqfile) <- rownames(tmp)
colnames(htseqfile) <- colnames(tmp)


indir <- "/data/Sherlock_Lung/RNASeq/Genelist/"

prefix.list <- c("Cell-2020-Maynard-General","Cell-2020-Maynard-Detailed")

for(i in 1:length(prefix.list)){
  prefix <- prefix.list[i]
  infilestr <- paste(indir,prefix,"_cell_type_markers.gmt",sep="")
  gmtfile<- GSA.read.gmt(infilestr)
  
  genelist <- {}
  MGlist <- list()
  for(k in  1:length(gmtfile$geneset.names)){
    genelist <- c(genelist,gmtfile$genesets[[k]][!gmtfile$genesets[[k]]==""])
    MGlist[[gmtfile$geneset.names[k]]] <- gmtfile$genesets[[k]][!gmtfile$genesets[[k]]==""]
  }
  tmp <- htseqfile[rownames(htseqfile)%in%genelist,]
  Aest <- AfromMarkers(tmp, MGlist)
  rownames(Aest) <- colnames(tmp)
  colnames(Aest) <- names(MGlist)
  
  outfilestr <- paste("/data/Sherlock_Lung/RNASeq/deconvolution/Sherlock1+2+TCGA-LUAD_htseq_Combat_TPM_",prefix,"_CAMTHC.txt",sep="")
  tmp <- t(Aest)
  out <- data.frame(tmp,check.names=FALSE) %>% mutate(signatures=rownames(tmp),.before=1)
  write_tsv(out,outfilestr)
  
}

############
### fibroblast cell type deconvolution using Bisque

library(tidyverse)
library(loomR)

library(Seurat)
library(SeuratDisk)
setwd("/data/Sherlock_Lung/scRNAseq/NatMed_Lambrechts/")
infilestr <- "/data/Sherlock_Lung/scRNAseq/NatMed_Lambrechts/Thienpont_Tumors_52k_v4_R_fixed.loom"

refData <-  connect(infilestr, mode = "r+")

library(Biobase)
library(BisqueRNA)
# gene list
n_genes <- refData[['row_attrs']][['Gene']][['dims']]
gns <- refData[['row_attrs']][['Gene']][1:n_genes]

# cell ID list
n_cells <- refData[['col_attrs']][['CellID']][['dims']]
cellids <- refData[['col_attrs']][['CellID']][1:n_cells]
ptsids <- refData[['col_attrs']][['PatientNumber']][1:n_cells]

# get raw counts matrix
#sc.counts.matrix <- as.data.frame(t(refData[['matrix']][,] ))

sc.counts.matrix <- t(refData[['matrix']][,] )

colnames(sc.counts.matrix) <- cellids
rownames(sc.counts.matrix) <- gns

#metadata <- data.frame(
#  cellID = cellids,
#  ClusterID = refData[['col_attrs']][['ClusterID']][1:n_cells],
#  annotation = refData[['col_attrs']][['ClusterName']][1:n_cells]
#)
#rownames(metadata) <- cellids

sample.ids <- colnames(sc.counts.matrix)
cell.type.labels <- refData[['col_attrs']][['ClusterName']][1:n_cells]
# individual.ids and cell.types should be in the same order as in sample.ids
sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=sample.ids,
                       SubjectName=ptsids,
                       cellType=cell.type.labels)
sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))
sc.pdata <- new("AnnotatedDataFrame",
                data=sc.pheno,
                varMetadata=sc.meta)


infilestr <- "/data/Sherlock_Lung/RNASeq/htseq/Sherlock1+2+TCGA-LUAD_RNAseq_htseq.genes.results.count.Combat-corrected_GeneName.txt"
tmp <- read.table(infilestr,sep="\t",header=T,row.names=1,check.names=FALSE)
bulk.matrix <- apply(tmp,2,as.numeric)
rownames(bulk.matrix) <- rownames(tmp)
colnames(bulk.matrix) <- colnames(tmp)

genelist <- intersect(rownames(sc.counts.matrix),rownames(bulk.matrix))
sc.counts.matrix <- sc.counts.matrix[genelist,]
bulk.matrix <- bulk.matrix[genelist,]

sc.eset <- Biobase::ExpressionSet(assayData=sc.counts.matrix,
                                  phenoData=sc.pdata)

bulk.eset <- Biobase::ExpressionSet(assayData = bulk.matrix)

res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, cell.types="cellType",subject.names="SubjectName",markers=NULL, use.overlap=FALSE)
save(res, file="/data/Sherlock_Lung/RNASeq/deconvolution/Sherlock1+r+TCGA-LUAD_htseq_Combat_Bisque_NatMed_Lambrechts.RData")

out <- data.frame( res$bulk.props,check.names=FALSE) %>% mutate(CellType=rownames( res$bulk.props),.before=1)
outfilestr <- "/data/Sherlock_Lung/RNASeq/deconvolution/Sherlock1+r+TCGA-LUAD_htseq_Combat_Bisque_NatMed_Lambrechts_detailed.txt"
write_tsv(out,outfilestr)
