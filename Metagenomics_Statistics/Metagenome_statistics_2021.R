#clean all the variables in the environment
rm(list=ls())
#################################################################
#load all the R packages we need
library(tidyr)
library(dplyr)
library(tibble)
library(phyloseq)
library(microbiome)
library(vegan)
#################################################################
#load the the Rdata file
load("phyloseq_object_data_1.Rdata")
################################################################
#extract the data in phyloseq object
#extract abundance profiles 
OTUdata <- abundances(phylo_obj)
#extract meta data
Metadata <- meta(phylo_obj)
#extract taxonomy of microbiome
TAXtable <- as.data.frame(tax_table(phylo_obj)@.Data)
################################################################
#alpha diversity estimators:richness
Adiv <- estimate_richness(phylo_obj,measures=c("Observed","Chao1","ACE"))

#############Task 1 here########################################
#############Task 1 here########################################
#alpha diversity estimators:richness and eveness
Adiv <- estimate_richness(phylo_obj, measures=c("Shannon","Simpson","InvSimpson"))
Adiv <- estimate_richness(phylo_obj, measures=c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson"))

#############Task 2 here########################################
#############Task 2 here########################################
#add a new column called 'SampleID' 
Adiv$SampleID <- rownames(Adiv)
Adiv_DF <- full_join(Metadata,Adiv,by="SampleID")#full_join function from dplyr

#normal distribution test:Shapiro-Wilk test
result<-shapiro.test(Adiv_DF$Shannon)
result
p_value=result$p.value
#Shannon is normally distributed. Chao1 is not.
#############Task 3 here########################################
#############Task 3 here########################################
#Difference analysis
#Compare shannon values between two groups (Treatment vs NoTreatment)
#Students' t-test
#get location (rows) of treatment
index_1<-which(Adiv_DF$Treatment=="Treatment")
#get location (rows) of NoTreatment
index_2<-which(Adiv_DF$Treatment=="NoTreatment")

#extract the shannon values of comorbidity samples
shannon_1<-Adiv_DF$Shannon[index_1]
#extract the shannon values of non-comorbidity samples
shannon_2<-Adiv_DF$Shannon[index_2]

#t test
ttest_result<-t.test(shannon_1,shannon_2)
ttest_result
t_stat<-ttest_result$statistic
p_ttest<-ttest_result$p.value

#Wilcoxon rank sum test
wilcox.test(shannon_1,shannon_2)

#############Task here#################
#############Task here#################

#build linear model for two estimators
#Continuous data
glmShannonCRPLevel =glm(CRPLevel ~ Shannon, data=Adiv_DF,family="quasipoisson")
summary(glmShannonCRPLevel)

#result<-summary(glmShannonCRPLevel)
#result$coefficients
#intercept<-result$coefficients[1,1]
#coef<-result$coefficients[2,1]
#p_value<-result$coefficients[2,"Pr(>|t|)"]

#plot(Adiv_DF$Shannon,Adiv_DF$CRPLevel)
#abline(glmShannonCRPLevel)

#############Task 4 here#################
#############Task 4 here#################
#correlation
cor.test(Adiv_DF$CRPLevel,Adiv_DF$Shannon,method="pearson")

##Beta diversity
phyobj_shift <- microbiome::transform(phylo_obj, transform="shift", shift=1)
phylo_Betdiv <- as.matrix(phyloseq::distance(phyobj_shift, method="bray"))
#head(phylo_Betdiv)
RE_MetaData_2 <- meta(phyobj_shift)

#PermANOVA analysis
adonis(phylo_Betdiv ~ Comorbidity, data=RE_MetaData_2, permutations=999)

#############Task here#################
#############Task here#################



#######################
#############################ddisease<-main$ER
p_cutoff<-0.05
health<-main[,3:8]#health samples' abundance
disease<-main[,3:8]#ldisease samples' abundance

t_list<-NULL
p_list<-NULL

for (i in 1:dim(disease)[1]){
  result=t.test(disease[i,],health[i,])
  t<-result$statistic
  p<-result$p.value
  t_list<-rbind(t_list,t)
  p_list<-rbind(p_list,p)
}

index_sig<-which(t_list>0&p_list<p_cutoff)

bac2<-as.matrix(unlist(main[62,3:8]))
SampleData_Rep$bac2<-bac2

#####################
library(ggplot2)
ggsave("plot.pdf",width=30,height=15, dpi=300)


pdf("heatmap.3.pdf",width = 20,height = 20)
heatmap2()
dev.off()