install.packages(c('readxl', 'tidyverse', 'dplyr', 'vegan', 'ggpubr'))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c('phyloseq', ask = F, update = F))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c('microbiome'))


library(readxl)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(phyloseq)
library(microbiome)
library(vegan)

main=read.table('merged_estimated_number_read.txt',header=TRUE,sep="\t")
t = separate(main,X.clade_name, into = c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species"), sep="\\|")
t=t[!is.na(t$Species), ]
t=t[!duplicated(t$Species),]

x=1:dim(t)[1]
OTUs= paste("OTU",x)
t=add_column(t, OTUs = OTUs, .before = "clade_taxid")

OTU =as.data.frame(t[,-which(names(t) %in% c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Type","clade_taxid"))])
row.names(OTU) = OTU$OTUs
OTU =as.data.frame(OTU[,-which(names(OTU) %in% c("OTUs","Species"))])
OTU[is.na(OTU)] = 0
names(OTU) = gsub(pattern = "_profile", replacement = "", x = names(OTU))
CountMatrix = OTU %>% as.matrix()
mode(CountMatrix) <- 'integer'


TAX =as.data.frame(t[,which(names(t) %in% c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species","OTUs"))])
row.names(TAX) = TAX$OTUs
TAX =as.data.frame(TAX[,-which(names(TAX) %in% c("OTUs"))])
TaxaMatrix <- TAX %>% as.matrix()

Metadata <- read_xlsx("KCL_META.xlsx", sheet=1) %>% as.data.frame()
rownames(Metadata) <- Metadata$SampleID



otuTABLE <- otu_table(CountMatrix, taxa_are_rows = TRUE)
taxTABLE <- tax_table(TaxaMatrix)
sampleDATA <- sample_data(Metadata)


phylo_obj <- phyloseq(otuTABLE, taxTABLE, sampleDATA)

### Extracting data from phyloseq object
OTUdata <- abundances(phylo_obj)
SampleData <- meta(phylo_obj)
TAXAData <- as.data.frame(tax_table(phylo_obj)@.Data)
