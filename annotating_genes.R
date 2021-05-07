#Annotating Lists
library(tidyverse)
library(svglite)
library(factoextra)
library(stringr)
library(stringi)
library(dplyr)





RIP<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/Copy of Supplementary Table S5_brief_ythdc1_ripseq_lui_nuc_acid_res_2020.csv")

x1<- read.csv(file="~/Desktop/Westmark Lab Diet Publications/GSE74397_YTHDC1_RIP_CuffDiff.csv")
write.csv(x1,file="~/Desktop/Westmark Lab Diet Publications/GSE74397_YTHDC1_RIP_CuffDiff_anno.csv")
x2<- read.csv(file="~/Desktop/Westmark Lab Diet Publications/GSE74397_Subcellular_RNA_Seq_genes.csv")
write.csv(x2,file="~/Desktop/Westmark Lab Diet Publications/GSE74397_Subcellular_RNA_Seq_genes_anno.csv")
x3<- read.csv(file="~/Desktop/Westmark Lab Diet Publications/GSE74397_Ythdc1KO_mESC_RNA-Seq_CuffDiff.csv")
write.csv(x3,file="~/Desktop/Westmark Lab Diet Publications/GSE74397_Ythdc1KO_mESC_RNA-Seq_CuffDiff_anno.csv")
mouse_genome<- read.csv(file="~/Desktop/Westmark Lab Diet Publications/mouse_annotation_list.csv", na.strings=c("","NA"))
mouse_genome$Gene<-mouse_genome$Approved.symbol
mouse_genome$gene_id<-mouse_genome$RefSeq.IDs
mouse_genome$tracking_id<-mouse_genome$RefSeq.IDs
#metabolism list
all_metab<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/all_metabolism_pathway_genes.csv")
list_metab <- c(all_metab[,2])

RIP2<-RIP%>%
  filter(gene_name %in% list_metab)

list_m<-mouse_genome%>%
 drop_na(RefSeq.IDs)%>%
  merge(all_metab, all.x=TRUE, all.y=TRUE, by="Gene")%>%
  drop_na(pathway)

x1_short<-x1_short%>%
  filter(status=="OK")%>%
  filter(significant=="yes")

x_all<-merge(x1_short, x2_short, all.y=TRUE, by="Gene.y")

write.csv(x_all,file="~/Desktop/Westmark Lab Diet Publications/GSE74397_Subcellular_RNA_Seq_genes_yRIP_filtered.csv")
