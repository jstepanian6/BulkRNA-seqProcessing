###################### II Part ######################
#cargar las librerias necesarias
library(ggplot2)
#library(dplyr)
#library(lattice)
#library(DESeq2)
library(WGCNA)

# Borra todos los objetos en la sesión de R
rm(list = ls())
# Display the current working directory
getwd();
dir= "/hpcfs/home/ing_sistemas/j.stepanian/PostMaster/PVulgaris_IteracionValerieyJohanna"
#dir = "/home/jstepanian/Documents/PostMaster/Pvulgaris_Valerie/Yang_ORourke_Astudillo_Vlasova_JefferyNoSoaked_fixed_name" #cambiar a la carpeta donde descargaron los archivos de entrada.
setwd(dir)

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

enableWGCNAThreads()

# Load the data saved in the first part
lnames = load(file = "Yang_ORourke_Astudillo_Vlasova_JefferyNoSoaked_dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 11, to=20, by=1))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file = "Vlasova_ORourke_Yang_Scale-free_topology.pdf", wi = 9, he = 6)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], cex=cex1,col="red")
dev.off()

#now calculate the adjacencies, using the soft thresholding power
softPower = 17;
adjacency = adjacency(datExpr0, power = softPower)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency, TOMType = "signed", verbose = 3);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
pdf(file = "GeneClustering.pdf")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
table(dynamicMods)


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
pdf(file = "geneDendro.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
pdf(file = "geneDendro2.pdf", wi = 9, he = 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
dev.off()

sizeGrWindow(12, 9)
pdf(file = "geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "networkConstruction-stepBySte.RData")
