#install.packages("BiocManager")
#BiocManager::install("WGCNA")
library(DESeq2) 
library(WGCNA)
library(genefilter)
library(tidyverse)

# Borra todos los objetos en la sesión de R
rm(list = ls())
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "/home/jstepanian/Documents/PostMaster/Pvulgaris_Valerie/Yang_ORourke_Astudillo_Vlasova_JefferyNoSoaked_fixed_name"
setwd(workingDir); 
## Load in data
data <- read.delim("geneCountMatrix_Yang_ORourke_Astudillo_Vlasova_JefferyNoSoaked_fixed_name.txt", sep = ",",header=TRUE, row.names=1, check.names=F) 
meta <- read.table("sra_Yang_ORourke_Astudillo_Vlasova_JefferyNoSoaked_fixed_name.txt",sep = "\t", header=T, row.names=1)

### Check that sample names match in both files
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ tissue)
##select ref level by comparation
#dds$type <- relevel(dds$type,ref="22DAP") 

dds <- DESeq(dds)

# remove  genes  without  any  counts
dds <- dds[ rowSums(counts(dds)) > 0, ]
###estimate size actors
dds <- estimateSizeFactors(dds)

vsd <- varianceStabilizingTransformation(dds)

expr_normalized <- getVarianceStabilizedData(dds)
head(expr_normalized)
dim(expr_normalized)
summary(expr_normalized)
write.csv(as.data.frame(expr_normalized), file="Yang_ORourke_Astudillo_Vlasova_JefferyNoSoaked_normalizated.csv" )
##PCA
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(vsd)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3=pca$x[,3])
ggplot(data = d, aes_string(x = "PC1", y = "PC2")) +  geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] *100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] *100), "% variance"))
pdf(file = "PCA_Yang_ORourke_Astudillo_Vlasova_JefferyNoSoaked.pdf", width = 12, height = 9);
plotPCA(vsd, intgroup="tissue")
dev.off()
################### 
#Establecer estos parametros es necesario para el funcionamiento del paquete
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

datExpr0 = as.data.frame(t(expr_normalized[, -c(1)]));
#We first check for genes and samples with too many missing values:
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples
#from the data:

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Next we cluster the samples (in contrast to clustering genes that will come later) 
# to see if there are any obvious outliers.

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "sampleClustering_Yang_ORourke_Astudillo_Vlasova_JefferyNoSoaked.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off();


#1.c Loading oil trait data
#traitData = read.csv("16DAP_FAMAT.csv");
traitData <- meta["tissue"]

dim(traitData)
names(traitData)
rownames(traitData)
allTrait=traitData
dim(allTrait)

# Form a data frame analogous to expression data that will hold the traits.
sample = rownames(datExpr0);
traitRows = match(sample,  rownames(meta));
#datTraits = traitData[traitRows,-1];
datTraits = traitRows
names(datTraits) = traitData[traitRows, 1];
#rownames(datTraits)= traitData[traitRows, 1] 
head(datTraits)

collectGarbage()

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed=FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = "Tissue",
                    main = "Sample dendrogram and Tissue heatmap")
sizeGrWindow(13,10)
pdf(file = "Dendograma_heatmap.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and Tissue heatmap heatmap")
dev.off()

#The last step is to save the relevant expression and trait data for use in the next steps 
save(datExpr0, datTraits, file = "Yang_ORourke_Astudillo_Vlasova_JefferyNoSoaked_dataInput.RData")
dim(datTraits)
dim(datExpr0)
dim(datTraits)

