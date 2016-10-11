##Read the pre-processed data
options(stringsAsFactors = FALSE)
library("survival")
library("sqldf")
library("qvalue")
library("glmnet")
library("WGCNA")


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

y=as.matrix(ty)

library("Coxnet")

#cv_r=cv.glmnet(x,y,family="cox")
#r_beta=cv_r$glmnet.fit$beta
#r_beta_m=as.matrix(r_beta)

res=Coxnet(x, y, penalty = "Enet",alpha = 0.1,nlambda=20,nfolds=10)
res=Coxnet(x, y, penalty = "Enet",alpha = 0.1,lambda=res$fit0$lambda)

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


library("spls")
  
  #setwd("F:/Job/project/omics/109/variation")

  #mRNA=read.table("D:/TCGA/RM_outlier/mrna_exp_rm_outlier.txt",sep=",")
  #miRNA=read.table("D:/TCGA/RM_outlier/mirna_exp_rm_outlier.txt",sep=",")
  #meth=read.table("D:/TCGA/RM_outlier/methy_combined_Level_3_JHU_USC__HumanMethylation27_rm_outlier.txt",sep=",")
  
  
   mRNA=read.table("/home/jzhang9/TCGA/mrna_exp_rm_outlier.txt",sep=",")
   miRNA=read.table("/home/jzhang9/TCGA/mirna_exp_rm_outlier.txt",sep=",")
   meth=read.table("/home/jzhang9/TCGA/methy_combined_Level_3_JHU_USC__HumanMethylation27_rm_outlier.txt",sep=",")
  
   
  
  lable_mRNA=NULL

  for (i in 1:(length(table(dynamicMods_0_02))-1)){
  
	subject_lable=c(rownames(meth),rownames(miRNA))[duplicated(c(rownames(meth),rownames(miRNA)))]
	lable_mRNA=lable_gene[which(dynamicMods_0_02==i),1]
	predict_mm=scale(cbind(miRNA[subject_lable,],meth[subject_lable,]))

  	set.seed(i)
  	cv <- cv.spls(predict_mm,scale(mRNA[subject_lable,lable_mRNA]), eta = seq(0.1,0.9,0.1), K = c(1:10) )

  	#SPLS fits are obtained as below.
  	f <- spls( predict_mm,scale(mRNA[subject_lable,lable_mRNA]), eta = cv$eta.opt, K = cv$K.opt )
  	#print(f)

  	gene_facors<-rownames(f$projection)

  	write.list(list(colnames(mRNA[subject_lable,lable_mRNA]),gene_facors),output=paste(i,"_module.txt",sep=""))

   }
