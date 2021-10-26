#####################################################
#####################################################
#rm(list=ls())
#Task 1 and 2
#load the the Rdata file
load("phyloseq_object_data_2.Rdata")

#extract abundance profiles 
OTUdata_2 <- abundances(phylo_obj_2)
#extract meta data
SampleData_2 <- meta(phylo_obj_2)
#extract taxonomy of microbiome
TAXAData_2 <- as.data.frame(tax_table(phylo_obj_2)@.Data)

#alpha diversity estimators:richness
Adiv_2 <- estimate_richness(phylo_obj_2,measures=c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson"))
#####################################################
#####################################################
#Task 3
Adiv_2$SampleID <- rownames(Adiv_2)
Adiv_DF_2 <- full_join(SampleData_2,Adiv_2,by="SampleID")#full_join function from dplyr

#normal distribution test
shapiro.test(Adiv_DF_2$Shannon)
#Shannon is normal distribution
##################################################
##################################################
#Task 4
index_1<-which(Adiv_DF_2$Condition=="Control-HFD")
index_2<-which(Adiv_DF_2$Condition=="Control-Chow")
value_1<-Adiv_DF_2$Shannon[index_1]
value_2<-Adiv_DF_2$Shannon[index_2]
wilcox.test(value_1,value_2)
#median(value_1)
#median(value_2)
##################################################
##################################################
#Task 5
cor.test(Adiv_DF_2$Shannon,Adiv_DF_2$Simpson,method="spearman")
##################################################
##################################################
#Task 6
#log10 transformation for the abundance data
phyobj_shift_2 <- microbiome::transform(phylo_obj_2,transform="log10")
#beta diversity based on euclidean distance 
Betdiv_2 <- as.matrix(phyloseq::distance(phyobj_shift_2, method="euclidean"))
#PermANOVA
adonis(Betdiv_2 ~ Condition, data=SampleData_2, permutations=999)

#if we only compare the beta diversity between Control-HFD and Control-Chow groups
index_1<-which(SampleData_2$Condition=="Control-HFD")
index_2<-which(SampleData_2$Condition=="Control-Chow")
index<-c(index_1,index_2)
SampleData_2<-SampleData_2[index,]

Betdiv_2<-Betdiv_2[SampleData_2$SampleID,SampleData_2$SampleID]
