
### This file shows the code for the geneoverlap analysis done in our paper 

### load needed packages 
library("GeneOverlap")

###
####Load MEP specific genes and LEP specifc genes 
####These are supplementary tables #4 and #11 respectively in the paper. 
MEP_Specific_Genes <- read.csv("MEP_SPecific_Genes.csv")
LEP_Specific_Genes <- read.csv("LEP_SPecific_Genes.csv")

#### give rownames of these gene lists based on their Ensembl.ID_Gene.Symbol 
rownames(MEP_Specific_Genes) <- MEP_Specific_Genes$Ensembl.ID_Gene.Symbol
rownames(LEP_Specific_Genes) <- LEP_Specific_Genes$Ensembl.ID_Gene.Symbol

### To test the gene overlap between MEP specific genes and the genes upregulated in high-risk (HR) LEPs comapred to average-risk (AR) LEPs 
#### take the genes upregulated in HR LEPs comapred to AR LEPs from the Differntially expressed genes from DESeq2 output 
HR_LEP_vs_NR_LEP_sig_upregulated <- subset(res05, padj < 0.05 & log2FoldChange >0 )

###now test for the gene overlap 
##this function shows the number of genes that are common between the two lists and creates the overalp object 
Overlap_object_LEP_general <- newGeneOverlap(rownames(MEP_Specific_Genes),rownames(HR_LEP_vs_NR_LEP_sig_upregulated))
###this function performs the statistical test of the gene overlap. It returns the Odd's ratio of the overlap and the p value 
Overlap_object_LEP_stat_test_general <- testGeneOverlap(Overlap_object_LEP_general)
Overlap_object_LEP_stat_test_general

#### load the LEP and MEP aging signatures. These are supplementary tables #6 and #7 respectively in the paper 
#### load the the senescence genes from the molecular signature data base 
#### select the genes upregualted in the senescence list by selecting the genes with log2foldchange>0 
Sen_up_Genes <- read.csv("sensecence_upregulated_geneset.csv")
Aging_genes_LEPs_upregulated <- read.csv("Aging_genes_in_LEP_upregulated.csv")

### to test for the overlap between genes upregulates in HR LEPs comapred to AR LEPs with the senescence gene signature 
Overlap_sen_genes_LEP <- newGeneOverlap(rownames(Sen_up_Genes),rownames(HR_LEP_vs_NR_LEP_sig_upregulated))
Overlap_sen_genes_LEP_stat_test <- testGeneOverlap(Overlap_sen_genes_LEP)

#### ### to test for the overlap between genes upregulates in HR LEPs comapred to AR LEPs with the aging signature in LEPs
Overlap_aging_genes_LEP <- newGeneOverlap(rownames(Aging_genes_LEPs_upregulated),rownames(HR_LEP_vs_NR_LEP_sig_upregulated))
Overlap_aging_genes_LEP_stat_test <- testGeneOverlap(Overlap_aging_genes_LEP)
