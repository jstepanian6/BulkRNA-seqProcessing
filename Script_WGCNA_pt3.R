###################### III Part ######################
# Borra todos los objetos en la sesión de R
rm(list = ls())

#cargar las librerias necesarias
library(ggplot2)
library(dplyr)
library(lattice)
#library(DESeq2)
library(WGCNA)
# Display the current working directory
getwd();
dir = "/hpcfs/home/ing_sistemas/j.stepanian/PostMaster/PVulgaris_IteracionValerieyJohanna"
setwd(dir)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
enableWGCNAThreads()
# Load the expression and trait data saved in the first part
lnames = load(file = "Yang_ORourke_Astudillo_Vlasova_JefferyNoSoaked_dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames

# Load network data saved in the second part.
lnames = load(file = "networkConstruction-stepBySte.RData");
lnames


# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
head(MEs)
moduleTraitCor = cor(MEs, datTraits, use = "p");
head(moduleTraitCor)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

write.csv(MEs, "MEs_Vlasova_ORourke_Yang.csv")

sizeGrWindow(2,15)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), " - (",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

names(datTraits)
dim(datTraits)
head(moduleTraitCor)


# Display the correlation values within a heatmap plot
#labeledHeatmap(Matrix = moduleTraitCor,
#               xLabels = "Tissue",
#              yLabels = names(MEs),
#              ySymbols = names(MEs),
#               colorLabels = FALSE,
#               colors = blueWhiteRed(50),
#               textMatrix = textMatrix,
#               setStdMargins = FALSE,
#               cex.text = 0.5,
#               zlim = c(-1,1),
#               main = paste("Module-trait relationships"))




# Define variable weight containing the weight column of datTrait
tissue = as.data.frame(datTraits);
names(tissue) = "Tissue"
# names (colors) of the modules
modNames = substring(names(MEs), 3)


geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr0, tissue, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
head(GSPvalue)

names(geneTraitSignificance) = paste("GS.", names(tissue), sep="");
names(GSPvalue) = paste("p.GS.", names(tissue), sep="");


#####################
propVarExplained(datExpr0,moduleColors ,MEs,corFnc = "cor")
dim(datExpr0)
head(datExpr0)


module = "darkolivegreen"
module="yellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for tissue",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



names(datExpr0)
moduleIDGenes= as.data.frame(names(datExpr0)[moduleColors==module])
dim(moduleIDGenes)

##
annot = read.csv(file = "annotation_Pvulgaris442.tsv", sep="\t")
dim(annot)
names(annot)
probes = names(datExpr0)
head(probes)
probes2annot =match(probes,annot$probe)


# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.


# Create the starting data frame
geneInfo0 = data.frame(geneSymbol =probes,
                       at =annot$probe[probes2annot],
                       anat=annot$locusName[probes2annot],
                       bestHitArabname=annot$Best.hit.arabi.name[probes2annot], 
                       arabi.symbol=annot$arabi.symbol[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for tissue
modOrder = order(-abs(cor(MEs, tissue, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance

geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Tissue));
geneInfo = geneInfo0[geneOrder, ]


write.csv(geneInfo, file = "geneInfo_Vlasova_ORourke_Yang_byTissue.csv")
##########
#Create input cytoscape
# Define the module of interest
modules = c("brown", "cyan", "black", "mediumorchid", "purple")  # Add more modules as needed

# Extract the topological overlap matrix (TOM)
TOM = TOMsimilarityFromExpr(datExpr0, power = 6)  # Adjust power if necessary

# Select the probes
probes = names(datExpr0)

# Loop through each module
for (module in modules) {
  # Select the genes in the current module
  moduleGenes = moduleColors == module

  # Get the module probes and annotations
  moduleProbes = probes[moduleGenes]
  moduleGenesAnnot = annot$gene_symbol[match(moduleProbes, annot$probe)]
  moduleTOM = TOM[moduleGenes, moduleGenes]

  # Set the dimnames for the TOM
  dimnames(moduleTOM) = list(moduleProbes, moduleProbes)

  # Export the network to Cytoscape
  cyt = exportNetworkToCytoscape(moduleTOM,
                                 edgeFile = paste("CytoscapeInput-edges-",
                                                  module, ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-",
                                                  module, ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = moduleProbes,
                                 altNodeNames = moduleGenesAnnot,
                                 nodeAttr = moduleColors[moduleGenes])
}

########## 

#hubGeneSignificance(datKME = datExpr0,GS = geneModuleMembership)
#head(GSPvalue)
#head(geneTraitSignificance)

#MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
#moduleGenes()

#automaticNetworkScreeningGS()




