

#m6a
m6a_both2<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/unique_in_m6a.csv")
m6a_final<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/final_m6a.csv")
#original_m6a<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/combined_Table S3_m6a_peaks_lg.csv")
new_m6a<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/m6a_filtered_by_metab_4conditions_final_overlap.csv")
#metabolism list
all_metab<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/all_metabolism_pathway_genes.csv")
list_metab <- c(all_metab[,2])

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


#graph m6a
#create filter from stats
stat_filter<-c(long_stat[,1])
short_stat<-c(long_stat2[,1])
#get Granges ready to merge
Granges3<-read.csv(file="~/Desktop/Westmark Lab Diet Publications/Granges2.csv")
Granges_merge<-Granges3%>%
  unite(range, start, end, sep = "-")%>%
  select(loc, range)
Granges_merge$loc<-as.character(Granges_merge$loc)

#ready for plot
long_2_plot<-new_m6a%>%
  filter(Gene== "FASN")%>%
  unite(ID, Gene, group)%>%
  #filter(ID %in% stat_filter)%>%
  separate(ID, c("Gene", "loc"), sep = "([_])")%>%
  group_by(Genotype, Cell, Gene, loc) %>% 
  summarise(avg.log2.enrich = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE))%>%
 # left_join(all_metab, by="Gene")%>%
  left_join(Granges_merge, by="loc")%>%
  mutate(name = Cell)%>%
  mutate(Name2=Genotype)%>%
  unite(Name, Name2, name)

#%>%
  #usually select by range not group and name
#  select(Genotype, Cell, Gene, Name, range, loc, avg.log2.enrich, sd)%>%
 # group_by(Gene)

pdf("~/Desktop/Westmark Lab Diet Publications/final data for paper/09032020_m6a_pro_only2.pdf")
for(i in 1:13)  
{
  print(ggplot(long_2_plot, aes(x=range, y=avg.log2.enrich, color=Name, fill = Name))+ 
          geom_col(position=position_dodge(), width=0.9) +
          geom_errorbar(aes(ymin=avg.log2.enrich-sd, ymax=avg.log2.enrich+sd), na.rm = TRUE, width=0.2, position=position_dodge(width=0.9))+
          #geom_text(aes(label=Cell), colour="white", size=3.5, position=position_dodge(0.5))+
          facet_wrap_paginate(~Gene,  nrow = 2, ncol = 1,scales = "free_x", page=i))
}
dev.off()


ggplot(long_2_plot, aes(x=range, y=avg.log2.enrich, color=Name, fill = Name))+ 
          geom_col(position=position_dodge(), width=0.9) +
          geom_errorbar(aes(ymin=avg.log2.enrich-sd, ymax=avg.log2.enrich+sd), na.rm = TRUE, width=0.2, position=position_dodge(width=0.9))


