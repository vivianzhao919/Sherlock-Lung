#################
### Multivariate prediction in the Sherlock/TCGA cohort (Training set)
###################
library(glmnet)
library(survival)
library(tidyverse)
setwd("C:\\zhaow\\Sherlock\\Script\\Submission\\")
infilestr <-  "Sherlock1+2+TCGA-LUAD_CPM_coding_xMT_Top5000_rank3_20features_prediction_dist.txt"
infile <- read_tsv(infilestr)
infile <- infile%>% select(!starts_with("Predicted_subtype")) %>% rename(Cluster1="steady",Cluster2="proliferative",Cluster3="chaotic")
#infile <- tmp%>% column_to_rownames(var="Samples")



infilestr <- "Sherlock_RNA_clinical_data.txt"
annfile <- read_tsv(infilestr)
annfile <- annfile  %>% filter(!is.na(Differentiation)) %>% filter(STAGE_group!=99999 & STAGE_group!=0) %>% filter(VITAL_STATUS_10yrs!=99999 & SURVIVAL_DAYS_10yrs!=0) %>% filter(AGE_AT_DIAGNOSIS!=99999)
annfile <- annfile%>% left_join(infile,by=join_by(SAMPLE_ID==Samples))


annfile <- annfile%>% mutate(STAGE=dplyr::recode(as.character(annfile$STAGE_group),"1"="I","2"="II","3"="III","4"="IV"))


dataPheno <- annfile%>% mutate(time=SURVIVAL_DAYS_10yrs/7) %>% mutate(status=recode(VITAL_STATUS_10yrs,"Deceased"=1,"Alive"=0)) %>% mutate(SEX=recode(SEX,"Male"=1,"Female"=2))
dataPheno$SURVIVAL<- Surv(dataPheno$time, dataPheno$status)

dataPheno$STAGE <- factor(dataPheno$STAGE,levels=c("I","II","III","IV"))  
dataPheno$SEX <-  factor(dataPheno$SEX,levels=c(1,2))
dataPheno$Class <- factor(dataPheno$PREDICTED_SUBTYPE)
dataPheno$Differentiation <- factor(dataPheno$Differentiation,levels=c("poorly","moderate","well"))

dataPheno$Cluster1_std <- scale(dataPheno$Cluster1)
dataPheno$Cluster2_std <- scale(dataPheno$Cluster2)
dataPheno$Cluster3_std <- scale(dataPheno$Cluster3)


summary(tempCoxph1 <- coxph(SURVIVAL ~ Differentiation, data = dataPheno))
summary(tempCoxph2 <- coxph(SURVIVAL ~ Differentiation + Cluster1_std+Cluster2_std+Cluster3_std, data = dataPheno))

pred1 <- predict(tempCoxph1)
pred2 <- predict(tempCoxph2)

anova(tempCoxph1,tempCoxph2)

out <- cbind(pred1, pred2)
colnames(out) <- c("Grade","Grade_Subtype")
outfile <- cbind(time=dataPheno$time,status=dataPheno$status,out)
outfile <- data.frame(outfile,check.names=FALSE) %>% mutate(Sample=dataPheno$SAMPLE_ID,.before=1)

outfilestr <- "Sherlock+TCGA-LUAD_pathology_riskscore_prediction.txt"
write_tsv(outfile,outfilestr)


#################
### Multivariate prediction in the Sherlock/TCGA cohort (Training set)
###################
library(glmnet)
library(survival)
library(tidyverse)


setwd("C:\\zhaow\\Sherlock\\Script\\Submission\\")
infilestr <-  "Sherlock1+2+TCGA-LUAD_CPM_coding_xMT_Top5000_rank3_20features_prediction_dist.txt"
infile <- read_tsv(infilestr)
infile <- infile%>% select(!starts_with("Predicted_subtype")) %>% rename(Cluster1="steady",Cluster2="proliferative",Cluster3="chaotic")

infilestr <-  "Sherlock_WGS_driver_mut_fusion_bin.txt"
infile.mut <- read_tsv(infilestr)



infilestr <- "Sherlock_RNA_clinical_data.txt"
annfile <- read_tsv(infilestr)
annfile <- annfile  %>%  filter(STAGE_group!=99999 & STAGE_group!=0) %>% filter(VITAL_STATUS_10yrs!=99999 & SURVIVAL_DAYS_10yrs!=0) %>% filter(AGE_AT_DIAGNOSIS!=99999)
annfile <- annfile %>% left_join(infile,by=join_by(SAMPLE_ID==Samples)) %>% inner_join(infile.mut,by=join_by(SAMPLE_ID))


annfile <- annfile%>% mutate(STAGE=dplyr::recode(as.character(annfile$STAGE_group),"1"="I","2"="II","3"="III","4"="IV"))


dataPheno <- annfile%>% mutate(time=SURVIVAL_DAYS_10yrs/7) %>% mutate(status=recode(VITAL_STATUS_10yrs,"Deceased"=1,"Alive"=0)) %>% mutate(SEX=recode(SEX,"Male"=1,"Female"=2))
dataPheno$SURVIVAL<- Surv(dataPheno$time, dataPheno$status)

dataPheno$STAGE <- factor(dataPheno$STAGE,levels=c("I","II","III","IV"))  
dataPheno$SEX <-  factor(dataPheno$SEX,levels=c(1,2))
dataPheno$Class <- factor(dataPheno$PREDICTED_SUBTYPE)
dataPheno$Differentiation <- factor(dataPheno$Differentiation,levels=c("poorly","moderate","well"))

dataPheno$Cluster1_std <- scale(dataPheno$Cluster1)
dataPheno$Cluster2_std <- scale(dataPheno$Cluster2)
dataPheno$Cluster3_std <- scale(dataPheno$Cluster3)



dataPheno$EGFR <- factor(dataPheno$EGFR,levels=c(0,1))
dataPheno$TP53 <- factor(dataPheno$TP53,levels=c(0,1))
dataPheno$KRAS <- factor(dataPheno$KRAS,levels=c(0,1))
dataPheno$EML4__ALK <- factor(dataPheno$EML4__ALK,levels=c(0,1))



summary(tempCoxph1 <- coxph(SURVIVAL ~  STAGE, data = dataPheno))
summary(tempCoxph2 <- coxph(SURVIVAL ~ EGFR + TP53 + KRAS + EML4__ALK, data = dataPheno))
summary(tempCoxph3 <- coxph(SURVIVAL ~ Cluster1_std+Cluster2_std+Cluster3_std, data = dataPheno))
summary(tempCoxph4 <- coxph(SURVIVAL ~ STAGE + EGFR + TP53 + KRAS + EML4__ALK, data = dataPheno))
summary(tempCoxph5 <- coxph(SURVIVAL ~ STAGE +Cluster1_std + Cluster2_std + Cluster3_std , data = dataPheno))
summary(tempCoxph6 <- coxph(SURVIVAL ~ EGFR + TP53 + KRAS + EML4__ALK + Cluster1_std + Cluster2_std + Cluster3_std, data = dataPheno))
summary(tempCoxph7 <- coxph(SURVIVAL ~ STAGE + EGFR + TP53 + KRAS + EML4__ALK + Cluster1_std + Cluster2_std + Cluster3_std, data = dataPheno))

pred1 <- predict(tempCoxph1)
pred2 <- predict(tempCoxph2)
pred3 <- predict(tempCoxph3)
pred4 <- predict(tempCoxph4)
pred5 <- predict(tempCoxph5)
pred6 <- predict(tempCoxph6)
pred7 <- predict(tempCoxph7)


anova(tempCoxph1,tempCoxph5)
anova(tempCoxph2,tempCoxph6)


out <- cbind(pred1,pred2,pred3, pred4, pred5, pred6, pred7)
colnames(out) <- c("Stage","Mutation","Subtype", "Stage_Mut","Stage_Subtype","Mut_Subtype","Stage_Mut_Subtype")
outfile <- cbind(time=dataPheno$time,status=dataPheno$status,out)
outfile <- data.frame(outfile,check.names=FALSE) %>% mutate(Sample=dataPheno$SAMPLE_ID,.before=1)

outfilestr <- "Sherlock+TCGA-LUAD_Genomic_Clinical_riskscore_prediction.txt"
write_tsv(outfile,outfilestr)



#############
infilestr <-"Sherlock+TCGA-LUAD_Genomic_Clinical_riskscore_prediction.txt"
infile <- read_tsv(infilestr)
indata <- infile[,4:ncol(infile)]
infile$SURVIVAL <- Surv(time=infile$time,event=infile$status)

cIndex <- {}
for(i in 1:ncol(indata)){
  
  cIndex <- c(cIndex,  1-concordance(infile$SURVIVAL ~ pull(indata,i))$concordance)
}
outfile <- data.frame(Model=colnames(indata),cIndex=cIndex)
outfilestr <-"Sherlock+TCGA-LUAD_Genomic_Clinical_riskscore_prediction_cIndex.txt"
write_tsv(outfile,outfilestr)







####################
### Multivariate prediction in the GIS cohort (Testing set)
library(glmnet)
library(survival)
library(tidyverse)


infilestr <- "GIS_CPM_coding_xMT_Top5000_rank3_20features_prediction_dist.txt"
infile <- read_tsv(infilestr)
infile <- infile%>% select(!starts_with("Predicted_subtype")) %>% rename(Cluster1="steady",Cluster2="proliferative",Cluster3="chaotic")


infilestr <-  "GIS_driver_mut_fusion_bin.txt"
infile.mut <- read_tsv(infilestr)

infilestr <- "GIS_clinical_data.txt"
annfile <- read_tsv(infilestr)
annfile <- annfile  %>%  filter(STAGE_group!=99999 & STAGE_group!=0) %>% filter(VITAL_STATUS_10yrs!=99999 & SURVIVAL_DAYS_10yrs!=0) %>% filter(AGE_AT_DIAGNOSIS!=99999)
annfile <- annfile %>% left_join(infile,by=join_by(SAMPLE_ID==Samples)) %>% inner_join(infile.mut,by=join_by(SAMPLE_ID))


annfile <- annfile%>% mutate(STAGE=dplyr::recode(as.character(annfile$STAGE_group),"1"="I","2"="II","3"="III","4"="IV"))


dataPheno <- annfile%>% mutate(time=SURVIVAL_DAYS_10yrs/7) %>% mutate(status=recode(VITAL_STATUS_10yrs,"Deceased"=1,"Alive"=0)) %>% mutate(SEX=recode(SEX,"Male"=1,"Female"=2))
dataPheno$SURVIVAL<- Surv(dataPheno$time, dataPheno$status)

dataPheno$STAGE <- factor(dataPheno$STAGE,levels=c("I","II","III","IV"))  
dataPheno$SEX <-  factor(dataPheno$SEX,levels=c(1,2))
dataPheno$Class <- factor(dataPheno$PREDICTED_SUBTYPE)

dataPheno$Cluster1_std <- scale(dataPheno$Cluster1)
dataPheno$Cluster2_std <- scale(dataPheno$Cluster2)
dataPheno$Cluster3_std <- scale(dataPheno$Cluster3)



dataPheno$EGFR <- factor(dataPheno$EGFR,levels=c(0,1))
dataPheno$TP53 <- factor(dataPheno$TP53,levels=c(0,1))
dataPheno$KRAS <- factor(dataPheno$KRAS,levels=c(0,1))
dataPheno$EML4__ALK <- factor(dataPheno$EML4__ALK,levels=c(0,1))


pred1 <- predict(tempCoxph1,newdata=dataPheno)
pred2 <- predict(tempCoxph2,newdata=dataPheno)
pred3 <- predict(tempCoxph3,newdata=dataPheno)
pred4 <- predict(tempCoxph4,newdata=dataPheno)
pred5 <- predict(tempCoxph5,newdata=dataPheno)
pred6 <- predict(tempCoxph6,newdata=dataPheno)
pred7 <- predict(tempCoxph7,newdata=dataPheno)


out <- cbind(pred1,pred2,pred3, pred4, pred5, pred6, pred7)
colnames(out) <-  c("Stage","Mutation","Subtype", "Stage_Mut","Stage_Subtype","Mut_Subtype","Stage_Mut_Subtype")
outfile <- cbind(time=dataPheno$time,status=dataPheno$status,out)
outfile <- data.frame(outfile,check.names=FALSE) %>% mutate(Sample=dataPheno$SAMPLE_ID,.before=1)

outfilestr <-  "GIS_Genomic_Clinical_riskscore_prediction.txt"
write_tsv(outfile,outfilestr)

#############
infilestr <-  "GIS_Genomic_Clinical_riskscore_prediction.txt"
infile <- read_tsv(infilestr)
indata <- infile[,4:ncol(infile)]
infile$SURVIVAL <- Surv(time=infile$time,event=infile$status)
cIndex <- {}
for(i in 1:ncol(indata)){
  cIndex <- c(cIndex,  1-concordance(infile$SURVIVAL ~ pull(indata,i))$concordance)
}
outfile <- data.frame(Model=colnames(indata),cIndex=cIndex)
outfilestr <-"GIS_riskscore_prediction_cIndex.txt"
write_tsv(outfile,outfilestr)

