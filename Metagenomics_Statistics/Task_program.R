#####################################################
#####################################################
#Task 1 
#load the the Rdata file
load("phyloseq_object_data_2.Rdata")

#extract abundance profiles 
OTUdata_2 <- abundances(phylo_obj_2)
#extract meta data
Metadata_2 <- meta(phylo_obj_2)
#extract taxonomy of microbiome
TAXtable_2 <- as.data.frame(tax_table(phylo_obj_2)@.Data)
#####################################################
#####################################################
#task 2
#alpha diversity estimators:richness
Adiv_2 <- estimate_richness(phylo_obj_2,measures=c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson"))
#####################################################
#####################################################
#Task 3
Adiv_2$SampleID <- rownames(Adiv_2)
Adiv_DF_2 <- full_join(Metadata_2,Adiv_2,by="SampleID")#full_join function from dplyr

#normal distribution test
shapiro.test(Adiv_DF_2$InvSimpson)
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
###
index_1<-which(Adiv_DF_2$Condition=="Blautia-HFD")
index_2<-which(Adiv_DF_2$Condition=="Control-Chow")
value_1<-Adiv_DF_2$InvSimpson[index_1]
value_2<-Adiv_DF_2$InvSimpson[index_2]
t.test(value_1,value_2)
