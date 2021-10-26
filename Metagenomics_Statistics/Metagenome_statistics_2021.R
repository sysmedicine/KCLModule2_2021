#Please click"Version control" and then "Pull branches"
#Go to the folder "Metagenomics_Statistics" and set this folder as working directory
#################################################################
#Clean all the variables in the environment
#rm(list=ls())
#################################################################
#Install the packages if you didn't installed before.Otherwise skip it
#install.packages("tidyr")#install tidyr
#install.package("dplyr")#install dplyr
#install.package("tibble")#install tibble

#install Vega
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Vega")

#install phyloseq
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("phyloseq")

#install microbiome
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("microbiome")
#################################################################
#load all the R packages
library(tidyr)
library(dplyr)
library(tibble)
library(vegan)
library(phyloseq)
library(microbiome)
#################################################################
#load the the Rdata file
load("phyloseq_object_data_1.Rdata")
################################################################
#extract abundance profiles 
OTUdata <- abundances(phylo_obj)
#extract meta data
SampleData <- meta(phylo_obj)
#extract taxonomy of microbiome
TAXAData <- as.data.frame(tax_table(phylo_obj)@.Data)
################################################################
#alpha diversity estimators:richness
Adiv <- estimate_richness(phylo_obj,measures=c("Observed","Chao1","ACE"))

#visualize ‘Adiv’
View(Adiv)

#############Task 1 here########################################
#############Task 1 here########################################
#alpha diversity estimators:richness and evenness
Adiv <- estimate_richness(phylo_obj, measures=c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson"))

#############Task 2 here########################################
#############Task 2 here########################################
#add a new column called 'SampleID' 
Adiv$SampleID <- rownames(Adiv)
Adiv_DF <- full_join(SampleData,Adiv,by="SampleID")#full_join function from dplyr

#normal distribution test:Shapiro-Wilk test
result<-shapiro.test(Adiv_DF$Chao1)
result
p_value=result$p.value
#Shannon is normally distributed. Chao1 is not.
#############Task 3 here########################################
#############Task 3 here########################################
#Difference analysis
#Compare shannon values between two genders (female vs male)
#Students' t-test
#get location (rows) of female
index_1<-which(Adiv_DF$Gender=="Female")
#get location (rows) of male
index_2<-which(Adiv_DF$Gender=="Male")

#extract the shannon values of female samples
shannon_1<-Adiv_DF$Shannon[index_1]
#extract the shannon values of male samples
shannon_2<-Adiv_DF$Shannon[index_2]

#t test
ttest_result<-t.test(shannon_1,shannon_2)
ttest_result
t_stat<-ttest_result$statistic
p_ttest<-ttest_result$p.value

#Wilcoxon rank sum test
wilcox.test(shannon_1,shannon_2)

#############Task 4 here#################
#############Task 4 here#################
#linear correlation
cor.test(Adiv_DF$Chao1,Adiv_DF$Simpson,method="spearman")
#############Task 5 here#################
#############Task 5 here#################
##Beta diversity
#log10 transformation for the abudance data
phyobj_shift <- microbiome::transform(phylo_obj,transform="log10")
#beta diversity based on eudlidean distance 
Betdiv <- as.matrix(phyloseq::distance(phyobj_shift, method="euclidean"))
#PermANOVA analysis
adonis(Betdiv ~ Treatment, data=SampleData, permutations=999)

#############Task 6 here#################
#############Task 6 here#################

