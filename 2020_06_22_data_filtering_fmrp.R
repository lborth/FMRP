
# Libraries ---------------------------------------------------------------

options(stringsAsFactors = FALSE)
library(svglite)
library(factoextra)
library(stringr)
library(stringi)
library(dplyr)
library(ggforce)
library(tidyverse)
library(rio)
#library(xlsx)
library(readxl)
library(data.table)

-------------------
#set wd
#  setwd("User/laura/Desktop/R script/FMRP")

#Import files ------------------------------------------
#rna
#RNAseq<- read.csv(file="~/Desktop/Westmark Lab Diet Publications/RNAseq of cytoplasmic and nuclear fractions_table5.csv")
RNAseq_ANOVA<- read.csv(file="~/Desktop/Westmark Lab Diet Publications/all stats/all_ANOVA2_rna.csv")
RNAseq<- read.csv(file="~/Desktop/Westmark Lab Diet Publications/RNAseq of cytoplasmic and nuclear fractions_table5.csv")
RNAseq_ANOVA<- read.csv(file="~/Desktop/Westmark Lab Diet Publications/all stats/all_ANOVA2_rna.csv")
#m6a
m6a_both2<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/unique_in_m6a.csv")
m6a_final<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/final_m6a.csv")
#original_m6a<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/combined_Table S3_m6a_peaks_lg.csv")
new_m6a<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/m6a_filtered_by_metab_4conditions_final_overlap.csv")
#protein
protein_mettl3<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/proteinlist_fmrp_mettl3.csv")
protein_fmrpki<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/proteinlist_fmrp_fmrp_mutant.csv")
#development
DEfraction<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/DEGenes by Fraction_adult_neonate_price_genomeres_2020.csv")
DEgene<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/DEGenes by GENE_adult_neonate_price_genomeres_2020.csv")

RNA_nuc<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/GSE121809_NUC_WT_vs_NUC_KO_Differential_Expression_edens_2019 cell reports.csv")
RNA_nuc$Gene.Symbol<-toupper(RNA_nuc$Gene.Symbol)
#metabolism list
all_metab<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/all_metabolism_pathway_genes.csv")
#darnell 2020 and 2011 overlap list
all_darnell<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/darnell_fmrp_targets_merged.csv")
  
--------------------------------------------------------------------------
##overlap and lists  
#  all_metab$tracking_id<-all_metab$Gene
#all_metab2<-left_join(all_metab, m6a_final, by="tracking_id")
#all_metab2<-all_metab2%>%
 # select(tracking_id, Discription, pathway)%>%
#  unique()

#ways to make lists and dataframes into filters; don't really know why some work better
#all_list$Genes<-toupper(a_list$Genes) 
#list_f<-unlist(apply(all_metab, 1, list), recursive = FALSE)
list_pro <- c(protein_all[,1])
list_metab <- c(all_metab[,2])
list_rna_anova<- c(RNAseq_ANOVA2final[,8])
#make a list-m6a hits in all 4 groups  
all_list<-inner_join(all_metab, m6a_both2, by = c( "Gene"= "tracking_id"))
list_all <- c(all_list[,2])

protein_all<-protein_all%>%
  filter(m6a.interplay.fmr1.mutant>1.2|fold.change.fmr1.mutant>1.2)
list_pro <- c(protein_all[,1])

#m_d_overlap<-all_metab%>%
#  merge(all_darnell)

DEfraction2<-DEfraction%>%
  filter(Symbol %in% list_metab)%>%
  filter(prenatal.FDR>0.05)%>%
  filter(adult.FDR<0.05)


DEGene4<-DEgene%>%
  filter(Symbol %in% list_metab)%>%
  filter(Nuc.FDR<0.05)%>%
  filter(Cyto.FDR<0.05)
DEGeness<-rbind(DEGene2,DEGene3,DEGene4)

RNA_nuc2<-RNA_nuc%>%
  filter(Gene.Symbol %in% list_metab)%>%
  filter(Significant=="yes")
write.csv(RNA_nuc2, file="~/Desktop/Westmark Lab Diet Publications/RNA_nuc_only_filtered_fdr05_metab.csv")

DEmerged<-DEGeness%>%
  unique()%>%
  merge(DEfraction2, by="Symbol", all.x=TRUE, all.y=TRUE)
write.csv(DEmerged, file="~/Desktop/Westmark Lab Diet Publications/DE_genes_filtered_multiple_fdr05_metab.csv")
m6a_2<-new_m6a%>%
  filter(Gene=="VCP")

RNAseq_ANOVA2final<-RNAseq_ANOVA%>%
  mutate(Padj.Genotype=pval_Factor2*359)%>%
  mutate(Padj.Cell=pval_Factor1*359)%>%
  mutate(Padj.pair=pval_pair*359)%>%
  mutate(Gene=tracking_id)%>%
  select(-tracking_id)%>%
  #filter(Padj.Genotype<0.05 & Padj.Cell<0.05)%>%
  filter(Gene %in% list_metab)
write.csv(RNAseq_ANOVA2final, file="~/Desktop/Westmark Lab Diet Publications/final data for paper/RNAseq_ANOVA2final.csv)


##
m6a_metab<- m6a_final%>%
  filter(tracking_id %in% list_all)%>%
  dplyr::rename(Gene=tracking_id)
write.csv(m6a_metab, file="~/Desktop/Westmark Lab Diet Publications/m6a_filtered_by_metab_4conditions.csv")

#number locations 
library(data.table) 
library(GenomicRanges)
x<-fread(file="~/Desktop/Westmark Lab Diet Publications/m6a_filtered_by_metab_4conditions.csv")
x<-rownames_to_column(x, var="overlap")
x$overlap<-as.numeric(x$overlap)
Granges<-x[, as.data.table(reduce(IRanges(start, end))), by = .(Gene)]
Granges2<-rownames_to_column(Granges, var="loc")

# convert data frame to GRanges
m6a.gr=makeGRangesFromDataFrame(x,
                                seqnames.field=c("Gene"),
                                start.field=c("start"),
                                end.field=c("end"))
Granges.gr=makeGRangesFromDataFrame(Granges2,
                                seqnames.field=c("Gene"),
                                start.field=c("start"),
                                end.field=c("end"))
#find overlap
hits<-findOverlaps(m6a.gr, Granges.gr)
#format
hit<-as.data.frame(hits)
hit<-hit%>%  
dplyr::rename(overlap = queryHits, group = subjectHits)

#attach overlap as new column
new_m6a<-merge(hit,x, by="overlap") 
#write.csv(new_m6a, file="~/Desktop/Westmark Lab Diet Publications/m6a_filtered_by_metab_4conditions_final_overlap.csv")


#match cases in 2 datasets
RNAseq$tracking_id<-toupper(RNAseq$tracking_id) 
df2<-m6a_final%>%
  filter(tracking_id %in% list)
#####################################################
#protein stuff
protein_all<-merge(protein_fmrpki, protein_mettl3, by="Gene", all = TRUE)
protein_all<-protein_all%>%
  filter(Gene %in% list_metab)
protein_overlap<-merge(protein_fmrpki, protein_mettl3, by="Gene")
protein_overlap<-protein_overlap%>%
  filter(Gene %in% list_all)

protein_graph<-protein_all%>%
  merge(all_metab, by="Gene")%>%
  #select(Gene, pathway, fold.change.fmr1.mutant, m6a.interplay.fmr1.mutant)%>%
  unite(Name, Gene, pathway)

#############################################################################################
#2way Anova
#format for anova
rna3<-RNAseq%>%
  gather(group, value, 2:17)%>%
  separate(group,into=c("Genotype", "Cell", "Sample"),sep = "[.]")%>%
  select(-"Sample")

m6a4<- df2%>%
  select(-Sample)
#make into dataframe and factors
m6a4<-as.data.frame(m6a4)
long_2 <-new_m6a%>%
  select(Gene, group, Cell, Genotype, value)%>%
  unite(ID, Gene, group)%>%
  unite(thing, Genotype, Cell)%>%
  group_by(ID)%>%
  summarise(count=n())%>%
  filter(count>4)
long_2<-as.data.frame(long_2)
list_short <- c(long_2[,1])
long_2 <-new_m6a%>%
  select(Gene, group, Cell, Genotype, value)%>%
  unite(ID, Gene, group)%>%
  filter(ID %in% list_short)
  
long_2<-rna3
long_2$Genotype<- as.factor(long_2$Genotype)
long_2$Cell<- as.factor(long_2$Cell)

# # Make into dataframe of lists
long_3 <- split(long_2, paste(long_2$ID))
tmpfn <- function(dat) {
  fit <- aov(value ~ Genotype + Cell  + Genotype:Cell, dat)
  sum_aov <- summary.aov(fit)
  pval_Genotype <- sum_aov[[1]][["Pr(>F)"]][1]
  pval_Cell <- sum_aov[[1]][["Pr(>F)"]][2]
  pval_pair <- sum_aov[[1]][["Pr(>F)"]][3]
  res <- (fit)$res
  data.frame(pval_Genotype,pval_Cell, pval_pair, res)
}

long_stat <- long_3 %>%
  map(tmpfn) %>%
  bind_rows(.id = "ID")

#filter 522 to 252 with p<0.5
long_stat2 <- long_stat %>%
  #select(-res) %>%
  unique()%>%
  filter(pval_Genotype<0.05|pval_Cell<0.05)

# long_interaction <- long_stat %>%
#   filter(pval_Cell <0.05)%>%
#   filter(pval_Genotype<0.05)
##################################################################################
# Plot the interaction 
# This script will make an interaction plot for each peptide that showed significant interaction between Factor1 and Factor2

#graph
long_stat_filter<-list<-c(long_interaction[,1])
# long_2_plot<-long_2%>%
#   filter(tracking_id %in% long_stat_filter)%>%
#   unite(ID, Gene, pathway)%>%
#   group_by(Genotype, Cell, ID) %>% 
#   summarise(value2 = mean(value, na.rm = TRUE),
#                sd = sd(value, na.rm = TRUE))

long_2$Gene<-toupper(long_2$tracking_id)
long_2_plot<-long_2%>%
  merge(m_d_overlap)%>%
  filter(Gene %in% list_rna_anova)%>%
  unite(ID, Gene, pathway)%>%
  group_by(Genotype, Cell, ID) %>% 
  summarise(value2 = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE))

pdf("08122020_rna_metab_sigx2.pdf")
for(i in 1:1)  
{
  print(ggplot(long_2_plot, aes(x=Genotype, y=value2, color=Cell, fill = Cell))+ 
      geom_col(position=position_dodge(), width=0.9) +
       geom_errorbar(aes(ymin=value2-sd, ymax=value2+sd), na.rm = TRUE, width=0.2, position=position_dodge(width=0.9))+ 
        facet_wrap_paginate(~ID,  nrow = 4, ncol = 2,scales = "free", page=i))
}
  dev.off()
  
#protein graph
#assemble and make long
   protein_graph2<-protein_graph%>%
     select(-Protein.names.x, -Protein.names.y, -pathway2, -11, -10)%>%
     rename(FC_fmr1/mut=fold.change.fmr1.mutant)%>%
     rename(m6a_fmr1/mut= m6a.interplay.fmr1.mutant)%>%
     rename(FC_mutant=foldChange_doxByNoDox.Mutant)%>%
     rename(FC_fmr1=foldChange_doxByNoDox.FMR1)%>%
     rename(m6a_FMR=m6a_foldChange_Dox.By.No.Dox.FMR_siMETTL3)%>%
     rename(m6a_mutant=m6a_foldChange_Dox.By.NoDox.Mutant_siMETTL3)%>%
     rename(m6a_si.ctl=m6a_foldChange_Dox.By.NoDox.siControl)%>%
   gather(key=group, value=fold.change, 2:8)
 
   protein_graph3<-protein_graph2%>% 
  mutate(exp2 = group)%>%
  separate(exp2, c("exp", "b"), sep = "_", extra="merge", remove=TRUE, convert=TRUE)%>%
     group_by(exp)%>%
     drop_na(fold.change)
   
   protein_graph3$newgroup<-stringr::str_wrap(protein_graph3$group, width = 7)
    
  pdf("07272020_protein_metab_all_justvalues_wide_high.pdf")
  for(i in 1:26)  
  {
    print(ggplot(protein_graph3, aes(x=newgroup, y=fold.change, color=exp, fill = exp))+ 
            geom_col(position=position_dodge(), width=0.3) +
            theme(legend.position="top")+
            ylim(0, 3)+
            geom_hline(yintercept=1)+
            facet_wrap_paginate(~Name,  nrow =4, ncol = 1,scales = "free", page=i))
  }
  dev.off()


#wrap long labels
  #aes(stringr::str_wrap(x, 15), y)
  
#line graph
aes(x = Factor1, y = value, color = Factor2) + 
  geom_line(aes(group = Factor2)) +
  geom_point() 


#graph m6a
#create filter from stats
stat_filter<-c(long_stat[,1])
short_stat<-c(long_stat2[,1])
#get Granges ready to merge
Granges_merge<-Granges2%>%
  unite(range, start, end, sep = "-")%>%
  select(loc, range)
#assemble using new_m6a
long_2_plot<-m6a_2%>%
  unite(ID, Gene, group)%>%
  #filter(ID %in% stat_filter)%>%
  separate(ID, c("Gene", "loc"), sep = "([_])")%>%
  group_by(Genotype, Cell, Gene, loc) %>% 
  summarise(avg.log2.enrich = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE))%>%
  #left_join(all_metab, by="Gene")%>%
  #filter(pathway2== "glucose")%>%
  #left_join(Granges_merge, by="loc")%>%
  mutate(name = Cell)%>%
  mutate(Name2=Genotype)%>%
  unite(Name, Name2, name)%>%
  #usually select by range not group and name
  select(Genotype, Cell, Gene, Name, range, loc, avg.log2.enrich, sd)%>%
  group_by(Gene)

pdf("08212020_m6a_vcp.pdf")
for(i in 1:1)  
{
  print(ggplot(long_2_plot, aes(x=range, y=avg.log2.enrich, color=Name, fill = Name))+ 
          geom_col(position=position_dodge(), width=0.9) +
          geom_errorbar(aes(ymin=avg.log2.enrich-sd, ymax=avg.log2.enrich+sd), na.rm = TRUE, width=0.2, position=position_dodge(width=0.9))+
          #geom_text(aes(label=Cell), colour="white", size=3.5, position=position_dodge(0.5))+
          facet_wrap_paginate(~Gene,  nrow = 5, ncol = 1,scales = "free_x", page=i))
}
dev.off()

geom_text(aes(label=Name), y=0.5, colour="white", size=3.5, position=position_dodge(0.9), check_overlap = T)+
#          facet_grid(.~Genotype)+

##################################################################################  
#stats
rna3_stats<-rna3%>%
group_by(Genotype,Cell, tracking_id)%>%
  summarise(rna_mean = mean(value, na.rm = TRUE),
            rna_sd = sd(value, na.rm = TRUE),
            count =n())%>%
  filter(tracking_id %in% list)

m6a_stats<-m6a4%>%
  group_by(Genotype,Cell, tracking_id)%>%
  summarise(m6a_mean = mean(value, na.rm = TRUE),
            m6a_sd = sd(value, na.rm = TRUE),
            count =n())%>%
  filter(tracking_id %in% list2)

#change to wide

new_m6a_short<-new_m6a%>%
  select(Gene)%>%
  unique()


final<-rna3_stats%>%
  unite(Sample2, Genotype, Cell, sep=".")%>%
  select(Sample2, tracking_id, rna_mean)%>%
  spread(key=Sample2, value=rna_mean)

final-m<-m6a_stats%>%
  unite(Sample2, Genotype, Cell, sep=".")%>%
  select(Sample2, tracking_id, rna_mean)%>%
  spread(key=Sample2, value=rna_mean)

final2<-rna3_stats%>%
  unite(Sample2, Genotype, Cell, sep=".")%>%
  select(Sample2, tracking_id, rna_sd)%>%
  spread(key=Sample2, value=rna_sd)

final3<-inner_join(final, final2, by="tracking_id")
final4<-inner_join(long_interaction, final3, by="tracking_id") 
##############################################################################3
#write file
write.csv(long_stat2_graph, file="~/Desktop/Westmark Lab Diet Publications/07152020_metab_rna_anova.csv")
#write results to xls
library("xlsx")
setwd("~/Desktop/Westmark Lab Diet Publications")

# Write the first data set in a new workbook
write.xlsx(glcy_in_rnaseq, file = "final_metabolism_fmrp_dataset.xlsx",
           sheetName = "glcy_in_rnaseq", append = FALSE)
# Add a second data set in a new worksheet
write.xlsx(metab_in_rnaseq, file = "metabolism_fmrp_dataset.xlsx", 
           sheetName="metab_in_rnaseq,", append=TRUE)
# Add a set in a new worksheet
write.xlsx(metab_in_m6a, file = "metabolism_fmrp_dataset.xlsx", 
           sheetName="metab_in_m6a,", append=TRUE)
# Add a  set in a new worksheet
write.xlsx(glcy_in_m6a, file = "metabolism_fmrp_dataset.xlsx", 
           sheetName="glyc_in_m6a,", append=TRUE)
#####################################################################################
###Notes
#to create final m6a
# 
#m6a <- excel_sheets("~/Desktop/Westmark Lab Diet Publications/Table S3_m6a_peaks.xlsx")
#make into database of lists
#list_all_example <- lapply(m6a, function(x) read_excel("~/Desktop/Westmark Lab Diet Publications/Table S3_m6a_peaks.xlsx", sheet = x))
#str(list_all_example)
#make into single database
#df2<-data.frame(Reduce(rbind, list_all_example))
#write combined file
#write.csv(df2, file="~/Desktop/Westmark Lab Diet Publications/combined_Table S3_m6a_peaks_lg.csv")
 df2_final <-long_m6a%>%
   drop_na(Symbol)%>% 
   separate(sample,into=c("Genotype", "Cell", "Sample"),sep = "[.]")%>%
   rename(gene=Symbol)%>%
   rename(value=R1.log2.Enrichment.)
 m6a_filtered<-inner_join(m6a_final, darnell_list, by="tracking_id")
write.csv(df2_filtered, file="~/Desktop/Westmark Lab Diet Publications/m6a_peaks_filtered_all_list_darnell.csv")
# df2_final_filtered<-df2_final%>%
#   filter(tracking_id %in% m6a_unique_list)
##to make unique list
# m6a_both<-m6a2%>%
#   select(-Sample)%>%
#   unite(yolo, Genotype, Cell, sep=".")%>%
#   group_by(yolo,tracking_id)%>%
#   summarise(protein_n=n())
# m6a_both2<-m6a_both%>%
#   group_by(tracking_id)%>%
#   summarise(protein_n=n())%>%
#   filter(protein_n>3)%>%
#   select(tracking_id)
#m6a_unique_list<-unlist(apply(m6a_both2, 1, list), recursive = FALSE)

m6a_stats<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/m6a_filtered_by_metab_4conditions_final_overlap.csv")

m6a<-m6a_stats%>%
  unite(ID, group, Gene)%>%
  unite(set, Genotype, Cell)%>%
  select(ID, set, value)%>%
    group_by(ID, set)%>%
    summarise(mean = mean(value, na.rm = TRUE))

m6a2<-m6a%>%
  spread(set, mean)%>%
  mutate(FC_KO_Cyto_Nuc = KO_Cyto/KO_Nuc, 
           FC_WT_Cyto_Nuc = WT_Cyto/WT_Nuc,
          FC_Cyto_KO_WT = KO_Cyto/WT_Cyto,
           FC_Nuc_KO_WT = KO_Nuc/WT_Nuc)%>%
  separate(ID, c("group", "Gene"), sep = "_", extra="merge", remove=TRUE, convert=TRUE)

write.csv(m6a2, file="~/Desktop/Westmark Lab Diet Publications/m6a_foldchange.csv")

