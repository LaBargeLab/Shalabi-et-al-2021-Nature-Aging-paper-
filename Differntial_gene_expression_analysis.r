#### set working directory to RNA seq data

directory <- 'Data'
setwd(directory)

### Load sample data and count data from the link we deposited our data at: XXXXX 

ShalabiRNASeq_LEP.MEP.sampleData<- data.frame(read.csv("ShalabiRNASeq_LEP.MEP.sampleData.csv"))
load(file=file.path(getwd(), 'ShalabiRNASeq_LEP.MEP.countData.RData'))

######define CellType. Risk 
ShalabiRNASeq_LEP.MEP.sampleData$CellType.Risk <- NA
ShalabiRNASeq_LEP.MEP.sampleData[ShalabiRNASeq_LEP.MEP.sampleData$CellType == "LEP" & 
                                   ShalabiRNASeq_LEP.MEP.sampleData$Risk == "High Risk", ]$CellType.Risk <- "LEP.High Risk" 
ShalabiRNASeq_LEP.MEP.sampleData[ShalabiRNASeq_LEP.MEP.sampleData$CellType == "LEP" & 
                                   ShalabiRNASeq_LEP.MEP.sampleData$Risk == "Normal Risk", ]$CellType.Risk <- "LEP.Normal Risk"
ShalabiRNASeq_LEP.MEP.sampleData[ShalabiRNASeq_LEP.MEP.sampleData$CellType == "MEP" & 
                                   ShalabiRNASeq_LEP.MEP.sampleData$Risk == "High Risk", ]$CellType.Risk <- "MEP.High Risk"
ShalabiRNASeq_LEP.MEP.sampleData[ShalabiRNASeq_LEP.MEP.sampleData$CellType == "MEP" & 
                                   ShalabiRNASeq_LEP.MEP.sampleData$Risk == "Normal Risk", ]$CellType.Risk <- "MEP.Normal Risk"

table(ShalabiRNASeq_LEP.MEP.sampleData$CellType.Risk)

# LEP.High Risk LEP.Normal Risk   MEP.High Risk MEP.Normal Risk 
# 11              10              11              10 

########## Remove all the genes or counts with a zero value 
dim(ShalabiRNASeq_LEP.MEP.countData)
# [1] 34623    42

AllZero.indx <- rowSums(ShalabiRNASeq_LEP.MEP.countData) == 0
length(which(AllZero.indx == TRUE))
# [1] 6824

countdata_clean <- ShalabiRNASeq_LEP.MEP.countData[!AllZero.indx, ]
dim(countdata_clean)
# [1] 27799    42

sampledata <-ShalabiRNASeq_LEP.MEP.sampleData
dim(sampledata)

####### make CellType.Risk into a factor with levels
sampledata$CellType.Risk <- factor(sampledata$CellType.Risk, levels=c("LEP.High Risk", "LEP.Normal Risk", "MEP.High Risk", "MEP.Normal Risk"))
sampledata$CellType.Risk

############## Differntial gene expression analysis in LEPs comapring high-risk vs average (normal) risk in DESeq2####

##### Create the dds object 

library("DESeq2")

ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countdata_clean,
                                            colData = sampledata,
                                            design = ~ CellType.Risk)


###### filter the genes with low counts 

keep <- rowSums(counts(ddsFullCountTable)) >= 10
dds <- ddsFullCountTable[keep,]


###### Now do the differntial gene expression analysis between HR and AR LEPs 
###########by default alpha or false discovery rate or adjusted p values is set to 0.1, so adjust alpha to 0.05 

res05 <- results(dds, contrast=c("CellType.Risk","LEP.High Risk", "LEP.Normal Risk"), alpha=0.05)
summary(res05)
res05

#### what is the number of differntially expressed genes with an adjusted p value <0.05?
sum(res05$padj < 0.05, na.rm=TRUE)

####save that list in the same directory or change the directory to the results folder 
DE1 <- subset(res05, res05$padj < 0.05)
write.csv(DE1, "DEgenes_LEPs_HR_VS_NR_res05_padj0.05.csv") 

### how many differntially expressed genes are with a log2 foldchange of 1 and -1 
sum(DE1$log2FoldChange>=1 | DE1$log2FoldChange<=-1)

###### to draw a volcano plot for the differnially expresed genes 
###
#### load the needed packages 
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

#### define the colors 
keyvals <- ifelse(
  res05$log2FoldChange < -1 & res05$padj <0.05, '#008B8B',
  ifelse(res05$log2FoldChange > 1 & res05$padj <0.05, '#FA8072',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == '#FA8072'] <- 'Upregulated'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == '#008B8B'] <- 'Downregulated'

#### now draw the volcano plot 
EnhancedVolcano(res05,
                lab = gsub(".*_","",rownames(res05)),       
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 4,
                ylim = c(0,6.2),
                xlim = c(-4,4),
                colCustom = keyvals
                )
                



