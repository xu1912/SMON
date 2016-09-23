##Read the pre-processed data
options(stringsAsFactors = FALSE)
library(survival)
library(sqldf)
library(qvalue)
library(glmnet)

dc=read.table("D:/TCGA/mRNA_ans/phe.csv",sep=",",header=T,na.strings="NA")
cdd=read.table("D:/TCGA/mRNA_ans/complete_data.txt",sep=",",header=T,row.names=1)

tdd=sqldf("select a.* from cdd a, dc b where a.V1=b.bcr_patient_barcode")

ty[,1]=ty[,1]/365.0
colnames(ty)[1]="time"
tx=as.matrix(cbind(tdd[,-1], dc$age_at_initial_pathologic_diagnosis))
tx=apply(tx, 2,as.numeric)
dc$age2=dc$age_at_initial_pathologic_diagnosis^2

my.surv <- Surv(ty[,1], ty[,2])
gene_ln=colnames(tx)[-17815]
p_n=c()
for (i in 1:length(gene_ln)){

	g_1=coxph(my.surv~age_at_initial_pathologic_diagnosis+age2+as.matrix(tx[,i]),data=dc)
	p_n[i]=summary(g_1)$coefficients[3,5]

}

q_n=qvalue(p_n)

y=as.matrix(y)

cv_r=cv.glmnet(x,y,family="cox")
r_beta=cv_r$glmnet.fit$beta
r_beta_m=as.matrix(r_beta)

wgcna_exp=x[,which(rowSums(r_beta_m)!=0)]
wgcna_exp=data.frame(wgcna_exp[,-length(wgcna_exp[1,])])


powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(wgcna_exp, powerVector = powers, verbose = 5)

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
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = F, minClusterSize=minModuleSize)

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

