##### This file here shows how the gene set enrichment analysis was done in our paper. 

### load the needed packages 
library("clusterProfiler")
library("enrichplot")
library("msigdbr")
library("dplyr")

####For human annotaion
library("org.Hs.eg.db", character.only = TRUE)

######### prepare input 
###keep only gene symbols and remove Ensemble_ID

####use the output we got from DESeq2 of the differntially expressed genes between high-risk and average-risk (normal risk) LEPs
res05_GeneSYM <- gsub(".*_","",rownames(res05))
rownames(res05) <- res05_GeneSYM
head(res05)
write.csv("res05_genesymb_for_gsea.csv")
df = read.csv("res05_genesymb_for_gsea.csv", header=TRUE)
df <- na.omit(df)
####need to vonvert into a vector 
rownames(df) <- df$X
original_gene_list <- as.vector(df$log2FoldChange)
names(original_gene_list) <- df$X
gene_list_LEP = sort(original_gene_list, decreasing = TRUE)

### use molecular signature database 
msigdbr_show_species()
msigdb_human = msigdbr(species = "Homo sapiens")
head(msigdb_human)
###choose hallmark MSigDB pathways 
msigdb_human_H_category = msigdbr(species = "Homo sapiens", category = "H")

msigdb_human_H_category_ready = msigdb_human_H_category %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

### now run GSEA 

LEP_HR_VS_NR_msigdb <- GSEA(gene_list_LEP, TERM2GENE = msigdb_human_H_category_ready)
head(LEP_HR_VS_NR_msigdb)
