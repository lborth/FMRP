# Libraries ---------------------------------------------------------------

options(stringsAsFactors = FALSE)
library(tidyverse)
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
#RNA
#import data
RNAseq<- read.csv(file="~/Desktop/Westmark Lab Diet Publications/RNAseq of cytoplasmic and nuclear fractions_table5.csv")
RNAseq$tracking_id<-toupper(RNAseq$tracking_id) 
RNAseq_ANOVA<- read.csv(file="~/Desktop/Westmark Lab Diet Publications/all stats/all_ANOVA2_rna.csv")

#import metabolism list
all_metab<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/all_metabolism_pathway_genes.csv")
#make into a list
list_metab <- c(all_metab[,2])

#filter, adj pval
RNAseq_ANOVA2final<-RNAseq_ANOVA%>%
  mutate(Padj.Genotype=pval_Factor2*336)%>%
  mutate(Padj.Cell=pval_Factor1*336)%>%
  mutate(Padj.pair=pval_pair*336)%>%
  mutate(Gene=tracking_id)%>%
  select(-tracking_id)%>%
  filter(Gene %in% list_metab)%>%
  merge(all_metab, by="Gene")
rna3<-RNAseq%>%
  filter(tracking_id %in% list_metab)%>%
  gather(group, value, 2:17)%>%
  separate(group, into=c("Genotype", "Cell", "Sample"), sep = "[.]")%>%
  merge(all_metab, by="Gene")%>%
  select(-tracking_id)
rna3$Genotype<- as.factor(rna3$Genotype)
rna3$Cell<- as.factor(rna3$Cell)
rna_stats_final<-rna3%>%
  mutate(Gene2=Gene)%>%
  unite(ID, Gene2, pathway)%>%
  group_by(Genotype, Cell, ID, Gene) %>% 
  summarise(value2 = mean(value, na.rm = TRUE),
            n=n(),
            sd = sd(value, na.rm = TRUE),
            se=sd/sqrt(n))
rna_all<-merge(rna_stats_final, RNAseq_ANOVA2final, by="Gene")

#write file
write.csv(rna_all, file="~/Desktop/Westmark Lab Diet Publications/final data for paper/RNAseq_all_final.csv")
write.csv(RNAseq_ANOVA2final, file="~/Desktop/Westmark Lab Diet Publications/final data for paper/RNAseq_ANOVA_final.csv")

#Plot

pdf("~/Desktop/Westmark Lab Diet Publications/final data for paper/09032020_rna_metab_all_44.pdf")
for(i in 1:44)  
{
  print(ggplot(rna_stats_final, aes(x=Genotype, y=value2, color=Cell, fill = Cell))+ 
          geom_col(position=position_dodge(), width=0.9) +
          geom_errorbar(aes(ymin=value2-se, ymax=value2+se), na.rm = TRUE, width=0.2, position=position_dodge(width=0.9))+ 
          facet_wrap_paginate(~Gene,  nrow = 4, ncol = 2,scales = "free", page=i))
}
dev.off()

#Stats
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
long_2 <-rna3%>%
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
long_3 <- split(long_2, paste(long_2$tracking_id))
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

#filter with p<0.5
long_stat2 <- long_stat %>%
  #select(-res) %>%
  unique()%>%
  filter(pval_Genotype<0.05|pval_Cell<0.05)

#protein#
#merge
protein_all<-merge(protein_fmrpki, protein_mettl3, by="Gene", all = TRUE)
#filter by FC
protein_all<-protein_all%>%
  filter(m6a.interplay.fmr1.mutant>1.2|fold.change.fmr1.mutant>1.2)
list_pro <- c(protein_all[,1])
#rename columns
  names(protein_all)[3] <- "m6a_FMR1-I304N/nodox"
  names(protein_all)[4] <- "m6a_FMR1OE/nodox"
  names(protein_all)[5] <- "m6a_comp.FMR1OE/I304N"
  names(protein_all)[7] <- "no-m6a_FMR1OE/nodox"
  names(protein_all)[8] <- "no-m6a_FMR1-I304N/nodox"
  names(protein_all)[9] <- "no-m6a_Si.ctl/nodox"
  names(protein_all)[12] <- "no-m6a_comp.FMR1OE/I304N"
    #5rename(FC_fmr1/mut=fold.change.fmr1.mutant)%>%
   #12 rename(m6a_fmr1/mut= m6a.interplay.fmr1.mutant)%>%
   #3rename(FC_mutant=foldChange_doxByNoDox.Mutant)%>%
   #4 rename(FC_fmr1=foldChange_doxByNoDox.FMR1)%>%
   #7 rename(m6a_FMR=m6a_foldChange_Dox.By.No.Dox.FMR_siMETTL3)%>%
    #8rename(m6a_mutant=m6a_foldChange_Dox.By.NoDox.Mutant_siMETTL3)%>%
    #9rename(m6a_si.ctl=m6a_foldChange_Dox.By.NoDox.siControl)
#filter, make long
    protein_graph<-protein_all%>%
    filter(Gene %in% list_metab)%>%
    merge(all_metab, by="Gene")%>%
    unite(Name, Gene, pathway)%>%
    select(-Protein.names.x, -Protein.names.y, -pathway2, -11, -10)
    write.csv(protein_graph, file="~/Desktop/Westmark Lab Diet Publications/final data for paper/prot_filt_metab.csv")
    
    
    %>%
    gather(key=group, value=fold.change, 2:8)%>%
  mutate(exp2 = group)%>%
    separate(exp2, c("exp", "b"), sep = "_", extra="merge", remove=TRUE, convert=TRUE)%>%
    group_by(exp)%>%
    drop_na(fold.change)%>%
      rename(Translation.Rate_Fold.Change=fold.change)

pdf("~/Desktop/Westmark Lab Diet Publications/final data for paper/protein_filter_fc12_ymax3.pdf")
for(i in 1:8)  
{
  print(ggplot(protein_graph, aes(x=b, y=Translation.Rate_Fold.Change, color=exp, fill = exp))+ 
          geom_col(position=position_dodge(), width=0.3) +
          theme(legend.position="top")+
          ylim(0, 3)+
          geom_hline(yintercept=1)+
          geom_hline(yintercept=1.2, linetype=2)+
          facet_wrap_paginate(~Name,  nrow =4, ncol = 1,scales = "free", page=i))
}
dev.off()

##other nuclear RNA test

RNA_nuc2<-RNA_nuc%>%
  filter(Gene.Symbol %in% list_metab)%>%
  filter(Significant=="yes")%>%
  unite(Gene, Gene.Symbol, Strand, sep = ".")%>%
  select(1, 16:21)%>%
  gather(key=group, value=value, 2:7)%>%
  separate(group, c("a", "b", "c","group", "sample"), sep = "_", extra="merge", remove=TRUE, convert=TRUE)%>%
  group_by(Gene, group)%>%
  summarise(value2 = mean(value, na.rm = TRUE),
            n=n(),
            sd = sd(value, na.rm = TRUE),
            se=sd/sqrt(n))
  

pdf("~/Desktop/Westmark Lab Diet Publications/final data for paper/09102020_nuc_only_rna_metab_sig.pdf")
for(i in 1:34)  
{
  print(ggplot(RNA_nuc2, aes(x=group, y=value2, color=group, fill = group))+ 
          geom_col(position=position_dodge(), width=0.9) +
          geom_errorbar(aes(ymin=value2-se, ymax=value2+se), na.rm = TRUE, width=0.2, position=position_dodge(width=0.9))+ 
          facet_wrap_paginate(~Gene,  nrow = 4, ncol = 2,scales = "free", page=i))
}
dev.off()
