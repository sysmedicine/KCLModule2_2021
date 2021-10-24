#########################Set the path and working environement########################
path='directory for data'
setwd("directory for data")

########################install mutiple packages at once##############################
list.of.packages <- c("ggplot2", "tidyverse", "reshape2","wesanderson","gplots","ggpubr",
                      "readxl","phyloseq","microbiome","dplyr","vegan","tibble","tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#################################Load packages########################################
library(readxl)
library(tidyr)
library(phyloseq)
library(ggplot2)
library(microbiome)
library(dplyr)
library(vegan)
library(tibble)
library(ggpubr)
library(wesanderson)
library(tidyverse)
library(reshape2)
library(gplots)
#################################pre-processing########################################
load(file = "phylo_obj.Rdata")

### Abundance bar plot
gp = subset_taxa(phylo_obj)
plot_bar(gp, fill = "Genus")

###################################Data visualization##################################
############################Prepare data for plotting##################################
pd <- psmelt(gp)
############################Barplot####################################################
# create basic bar plot
ggplot(data=pd, aes(x=SampleID, y=CRPLevel)) +
  geom_bar(stat="identity")

# change the width of the bar
ggplot(data=pd, aes(x=SampleID, y=CRPLevel)) +
  geom_bar(stat="identity", width = 0.8)

# change the color of the bar
ggplot(data=pd, aes(x=SampleID, y=CRPLevel)) +
  geom_bar(stat="identity", fill="steelblue")

# Visualization data using bar plot by group
ggplot(data=pd, aes(x=SampleID, y=CRPLevel, fill = Group)) +
  geom_bar(stat="identity")

# group samples by gender
ggplot(data=pd, aes(x=Gender, y=CRPLevel, fill = Group)) +
  geom_bar(stat="identity") 

# mimic plot in fig2b in https://www.nature.com/articles/s41586-019-1560-1
ggplot(data=pd, aes(x=Gender, y=Abundance, fill = Group)) + # define which columns you will display in plot
  geom_bar(stat="identity") + # plotting bar plot
  labs(x = "Gender groups", y = "Relative mean abundance (%)",
       title = "Bar plot of abundance by age groups", fill = "Age groups") + # define names for each axis and title
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) # define the style of plot title. 


#position = "dodge" in order put the bars side by side (instead of on top of each other)
ggplot(data=pd, aes(x=Gender, y=CRPLevel, fill = Group)) +
  geom_bar(stat="identity", position = "dodge")

# mimic plot in fig2a in https://www.nature.com/articles/s41586-019-1560-1
ggplot(data=pd, aes(x=Gender, y=CRPLevel, fill = Group)) + # define which columns you will display in plot
  geom_bar(stat="identity", position = "dodge") + # plotting bar plot
  labs(x = "Gender groups", y = "Relative mean abundance (%)",
       title = "Bar plot of abundance by age groups", fill = "Age groups") + # define names for each axis and title
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) + # define the style of plot title. 
  coord_flip() # Horonital the plot

# let us make a decent plot like Fig2 in https://www.nature.com/articles/s41586-019-1560-1
# ggarrange(fig2a, fig2b, nrow = 2, labels = c("A","B"), common.legend = TRUE)

############################Boxplot####################################################
# Basic boxplot
ggboxplot(SampleData, x = "Treatment", y = "CRPLevel")
# fill color for boxplot
ggboxplot(SampleData, x = "Treatment", y = "CRPLevel", fill = "Treatment")
# change color to what you like
ggboxplot(SampleData, x = "Treatment", y = "CRPLevel", fill = "Treatment",
          palette =c("#00AFBB", "#E7B800"))
# add dots
# Dots can be added to a box plot using the parameter add = "jitter" :
ggboxplot(pd_filter, x = "Phylum", y = "Abundance", 
          color = "Treatment", add = "jitter")+
  stat_compare_means()
# Add p-values comparing groups
ggboxplot(Adiv_DF, x = "Treatment", y = "Chao1", notch = TRUE,
          color = "Treatment",palette =c("#00AFBB", "#E7B800"), add = "jitter") +
  stat_compare_means(method = "t.test", label.x = 1.3, label = "p")
############################volinplot####################################################
pd_filter <- read.table(file = 'pd.txt', sep = "\t",header = TRUE)

my_comparisons <- list( c("day 4", "day 7"), c("day 4", "day 21"), c("day 7", "infancy") )
ggviolin(pd_filter, x = "Group", y = "CRPLevel", fill = "Group",
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")
# change x axis
#pd_filter <- pd %>% filter(Abundance == 0)
ggboxplot(pd_filter, x = "Phylum", y = "Abundance", 
          color = "Group")

############################Heatmap####################################################
heatmap.2(OTUdata)
# add parameters to control the outlook of plots.
# read the help() information for heatmap.2 function, which has a detailed 
# description of each parameters
#add whatever you want, for example, change the color, add title.
heatmap.2(OTUdata,scale = "none",col = bluered(100), 
          trace = "none", density.info = "none")


########################Change fill color of your plot##############################
ggboxplot(SampleData, x = "Treatment", y = "CRPLevel", color = "Treatment",add = "jitter")   
#In R, colors can be specified either by name (e.g col = “red”)
ggboxplot(SampleData, x = "Treatment", y = "CRPLevel", color = "Treatment",add = "jitter",
          palette = c("red", "blue"))  
#as a hexadecimal RGB triplet (such as col = “#FFCC00”). 
ggboxplot(SampleData, x = "Treatment", y = "CRPLevel", color = "Treatment",add = "jitter",
          palette = c("#00AFBB", "#E7B800"))  
#You can also use other color systems such as ones taken from the RColorBrewer package.
# First, how many colors you need for your data
# for example, when we plot boxplot for group, we need 4 colors
length(levels(as.factor(pd_filter$Group)))
#4
# list the name of palattes and pick 4 colors from one palatte, and then save.
names(wes_palettes)
colors <- wes_palette(n=4, name="Darjeeling1")

ggboxplot(pd_filter, x = "Phylum", y = "Abundance", 
          color = "Group", palette = colors)

#When the number of colors you need are more that the number of palattes can provide
#for example, when you want to plot Genus, you need
length(levels(as.factor(TAXAData$Genus)))
# 45
colors <- wes_palette(n=8, name="Darjeeling1")
# What kind of error you get?
# add the parameter: type = "continuous" and try again.
colors <- wes_palette(n=45, name="Darjeeling1", type = "continuous")
ggplot(data=pd, aes(x=Gender, y=Abundance, fill = Genus)) +
  geom_bar(stat="identity", position = "dodge") +
  scale_fill_manual(values = colors)
# for stacked bar plot
ggplot(data=pd, aes(x=Gender, y=Abundance, fill = Genus)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = colors)

save.image(file = "Metaboime_study.Rdata")
