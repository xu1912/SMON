options(stringsAsFactors = FALSE)
library("survival")
library("sqldf")
library("qvalue")
library("glmnet")
library("WGCNA")

#### 1.Read the pre-processed data
## Read phenotype data in csv format. Header included but can be customized.
dc=read.table("D:/TCGA/mRNA_ans/phe.csv",sep=",",header=T,na.strings="NA")

## Read the mRNA expression data in txt format. Each row represents a subject. Each column represents a gene. Header included but can be customized.
cdd=read.table("D:/TCGA/mRNA_ans/complete_data.txt",sep=",",header=T,row.names=1)

## Combine phenotype and expression data by the common id.
tdd=sqldf("select a.* from cdd a, dc b where a.V1=b.bcr_patient_barcode")

## In our manuscript, we used Cox survival model. You may change other model, lm or glm.
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


#### 2.Elastic net penalized regression for Cox model. Use glmnet for lm and glm model.
library("Coxnet")
y=as.matrix(ty)
x=as.matrix(tdd[,-1])

## Here, alpha=0.1 was also resulted from a tunning process.
res=Coxnet(x, y, penalty = "Enet",alpha = 0.1,nlambda=20,nfolds=10)
res=Coxnet(x, y, penalty = "Enet",alpha = 0.1,lambda=res$fit0$lambda)

#cv_r=cv.glmnet(x,y)
#r_beta=cv_r$glmnet.fit$beta
#r_beta_m=as.matrix(r_beta)

#### 3.WGCNA to identify modules in the trait-associated mRNAs.
wgcna_exp=x[,which(rowSums(r_beta_m)!=0)]
wgcna_exp=data.frame(wgcna_exp[,-length(wgcna_exp[1,])])

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(wgcna_exp, powerVector = powers, verbose = 5)

## Power=4 was resulted following WGCNA instructions.
## We tried two ways provided by the WGCNA to get the modules.
# 3.1
bwnet = blockwiseModules(wgcna_exp, maxBlockSize = 10000,
                           power = 4, minModuleSize = 10,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE,
                           saveTOMs = TRUE,
                           verbose = 3)
# 3.2
adjacency=adjacency(wgcna_exp,power=4)
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
minModuleSize = 15
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = F, minClusterSize=minModuleSize)

table(dynamicMods)

#### 4.SPLS to get module gene-associated regulators.
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
