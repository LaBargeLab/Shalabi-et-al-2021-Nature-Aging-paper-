
### This file provides the R code for producing the results in the cKit progenitor differentiation assays.

### Load file of the quantified mean KRT19 and KRT14 signals from the fixed and stained cKit progenitor cells at 2 or 7 days. 
### This file is provided as a source data file for figure 3. 

###set working directroy 
### load output file from CellProfiler 

Data <- data.frame(read.csv("ckit_analysis_final_07_14_2021.csv"))

##Load needed packages
library(ggpointdensity)
library(ggplot2)
library(dplyr)

#subset data by taking out empty images, false positive cells (small area, no intensity)
Data = subset(Data, 
                 Intensity_MeanIntensity_K19 > 0 &
                 Intensity_MeanIntensity_K14 > 0  );
###The file provided was saved after this step 

###### To draw histograms 
library(scales)
Data$norm_K19 <- rescale(Data$Intensity_MeanIntensity_K19)
Data$norm_K14 <- rescale(Data$Intensity_MeanIntensity_K14)
Data$k19_k14 <- log2(Data$norm_K14/Data$norm_K19)
### get the log2(k14/k19)
Data$k14_k19 <- log2(Data$Intensity_MeanIntensity_K14/Data$Intensity_MeanIntensity_K19)

###subset data into a 2day and 7 day time points 
x <- c("7DAYS", "7days")
Data_7days <- subset(Data, grepl(paste(x, collapse= "|"), Data$Metadata_strain_name));
Data_48H <- subset(Data, grepl(paste("48H", collapse= "|"), Data$Metadata_strain_name));

h <- ggplot(Data_48H)+ aes(x=Data48H$k19_k14, y = ..count../sum(..count..)*100, fill= Group)+
  geom_histogram(col="black", breaks=seq(-7.5, 7.5, by=0.25))+ 
  theme(legend.position="top")+
  ylim(c(0,25))+ ylab("Percentage of cells")+
  xlab("log2(K14/K19)")+
  scale_fill_manual(values = c("red" = "red", "green" = "green", "yellow" = "yellow"))

##### to remove gray grid and background 
h + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

### to draw density contour plots 
library(tidyverse)
library(RColorBrewer)
library(viridis)
display.brewer.all()

### for the 2 day timepoint 
ggplot(Data_48H, aes(Intensity_MeanIntensity_K14, Intensity_MeanIntensity_K19)) +
  geom_density_2d_filled(contour_var = "ndensity", h = 0.25) +
  xlab("K14") +
  ylab("K19") +
   facet_wrap(~ Mutation)

### for the 7 day timepoint 
ggplot(Data_7days, aes(Intensity_MeanIntensity_K14, Intensity_MeanIntensity_K19)) +
  geom_density_2d_filled(contour_var = "ndensity", h = 0.25) +
  xlab("K14") +
  ylab("K19") +
  facet_wrap(~ Mutation)

