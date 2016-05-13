options(stringsAsFactors = FALSE)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(wgcna_exp, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow =r c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");

plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

bwnet = blockwiseModules(wgcna_exp, maxBlockSize = 10000,
                           power = 4, minModuleSize = 10,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE,
                           saveTOMs = TRUE,
                           verbose = 3)

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(bwnet$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE,
addGuide = TRUE)



adjacency=adjacency(wgcna_exp,power=4)
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
minModuleSize = 15
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = F, minClusterSize=minModuleSize)

table(dynamicMods)


###SPLS
library("spls")
  
  setwd("F:/Job/project/omics/109/variation")


  set.seed(1)
  cv <- cv.spls( mice$x, mice$y, eta = seq(0.1,0.9,0.1), K = c(1:5) )

  #SPLS fits are obtained as below.
  f <- spls( mice$x, mice$y, eta = cv$eta.opt, K = cv$K.opt )
  print(f)

  gene_facors<-rownames(f$projection)








