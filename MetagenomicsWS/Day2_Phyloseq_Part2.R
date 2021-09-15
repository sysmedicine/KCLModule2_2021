library(readxl)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(phyloseq)
library(microbiome)
library(vegan)

main=read.table('Data/merged_abundance_table.txt',header=TRUE,sep="\t")
t = separate(main, ID, into = c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species","Type"), sep="\\|")
t=t[!is.na(t$Species), ]
t=t[!duplicated(t$Species),]

x=1:dim(t)[1]
OTUs= paste("OTU",x)
t=add_column(t, OTUs = OTUs, .before = "s1")

OTU =as.data.frame(t[,-which(names(t) %in% c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Type"))])
row.names(OTU) = OTU$OTUs
OTU =as.data.frame(OTU[,-which(names(OTU) %in% c("OTUs","Species"))])
CountMatrix = OTU %>% as.matrix()
mode(CountMatrix) <- 'integer'

TAX =as.data.frame(t[,which(names(t) %in% c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species","OTUs"))])
row.names(TAX) = TAX$OTUs
TAX =as.data.frame(TAX[,-which(names(TAX) %in% c("OTUs"))])
TaxaMatrix <- TAX %>% as.matrix()


Metadata <- read.table('Data/metadata_Task2.txt',header=TRUE,sep="\t")
rownames(Metadata) <- Metadata$SampleID


otuTABLE <- otu_table(CountMatrix, taxa_are_rows = TRUE)
taxTABLE <- tax_table(TaxaMatrix)
sampleDATA <- sample_data(Metadata)


phylo_obj <- phyloseq(otuTABLE, taxTABLE, sampleDATA)

### Extracting data from phyloseq object
OTUdata <- abundances(phylo_obj)
SampleData <- meta(phylo_obj)
TAXAData <- as.data.frame(tax_table(phylo_obj)@.Data)

### Abundance bar plot
gp = subset_taxa(phylo_obj)
plot_bar(gp, fill="Species")

save(list = ls(all.names = TRUE),file = 'Results/Day2_Part2.Rdata')
