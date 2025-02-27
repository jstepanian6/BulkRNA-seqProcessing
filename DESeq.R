#Install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
install.packages("Hmisc")

#Load packages
library(EnhancedVolcano)
library(DESeq2)
library(ggplot2)
library(rgl)
library(plotly)

#set working directory 
setwd("/home/jstepanian/Documents/Fertilidad/DESeq_genes")

#Load files 
counts <- read.csv("geneCountMatrix_hg38_genes.csv", sep = ",", row.names = 1)
head(counts)
metadata<- read.csv("metadata.tsv", sep="\t", row.names = 1)
all(colnames(counts)==rownames(metadata))
ncol(counts) == nrow(colData) 

#Create a DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata, 
                              design=~Categoría)
vsd <-  vst(dds, blind=F)
plotPCA(vsd, intgroup=c("Categoría", "Paciente.1"))
plotPCA(vsd, intgroup=c("Otras.características"))

rld <- rlog(dds)
plotPCA(rld, intgroup=c("Categoría", "Paciente.1"))

# mean vs variance 
vstNormalized <-(assay(vsd))

row_means <- rowMeans(vstNormalized)  
row_var <- rowVars(vstNormalized)

df <- data.frame(rownames(vstNormalized),row_means, row_var)
ggplot(df, aes(x = row_means, y = row_var)) +
  geom_point() +  labs(x = "Mean counts", y = "Variance counts") +
  theme_minimal()  

##What genes have huge variance (>10)?
which(row_var>10) 
gene <- (data.frame(vstNormalized["gene-LIF",]))
colnames(gene) <- "LIF"
ggplot(gene, aes(x= rownames(gene) , y=LIF))+ geom_point()

plot(vstNormalized["gene-DEPP1",])

##Filtering zero count genes from the matrix 
filtered_counts <- counts[rowSums(counts != 0) > 0, ]
dds <- DESeqDataSetFromMatrix(countData=filtered_counts, 
                              colData=metadata, 
                              design=~Categoría)
dds$Group <- relevel(dds$Categoría, ref = "Control")
vsd <-  vst(dds)
plotPCA(vsd, intgroup="Categoría")
#Export PCA coordinates
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(vsd)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3=pca$x[,3])
ggplot(data = d, aes_string(x = "PC1", y = "PC2")) +  geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] *100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] *100), "% variance"))

plot_ly(data = d, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "markers") %>%
  layout(scene = list(xaxis = list(title = paste0("PC1: ", round(percentVar[1] * 100), "% variance")),
                      yaxis = list(title = paste0("PC2: ", round(percentVar[2] * 100), "% variance")),
                      zaxis = list(title = paste0("PC3: ", round(percentVar[3] * 100), "% variance"))))

ggplot(data = d, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  scale_color_manual(values = rainbow(k)) +  # Change colors if needed
  theme_minimal()

#Differential expresion analysis
dds <- DESeq(dds)
resultsNames(dds)
plotDispEsts(dds)
res <- results(dds)
EnhancedVolcano(res,
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                #pCutoff = 10e-32,
                FCcutoff = 15.0,
                pointSize = 1.0,
                labSize = 5.0,
                colAlpha = 1,
                legendPosition = 'bottom',
                subtitle = "")
write.csv(res, "controlvspatient.csv")
