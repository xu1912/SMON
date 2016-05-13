##Read the pre-processed data
comp_d=(as.matrix(read.table("D:/TCGA/RM_outlier/mrna_exp_rm_outlier.txt",sep=",",header=T,row.names=1)))

pd=read.table("D:/TCGA/mRNA_ans/phe.csv",sep=",",header=T,stringsAsFactor=F,na.strings="NA")

pd[which(is.na(pd$days_to_death)),9]=pd[which(is.na(pd$days_to_death)),8]
pd$status=0
pd[which(pd$days_to_last_followup>=pd$days_to_death & pd$vital_status=="Dead"),11]=1

library(survival)

dc=pd[match(rownames(comp_d), pd$bcr_patient_barcode),c(-3,-5)]

ty=dc[,c(6,9)]
ty[,1]=ty[,1]/365.0
colnames(ty)[1]="time"
tx=as.matrix(cbind(comp_d, dc$age_at_initial_pathologic_diagnosis))
tx=apply(tx, 2,as.numeric)


colnames(tx)
gene_in=read.table("D:/TCGA/mRNA_ans/all_pre_module_gene.txt",header=F,stringsAsFactors=F)
for (i in 1:length(gene_in$V1)){

	tx[,c(which(colnames(tx)==gene_in$V1[i]),17815)]

}

my.surv <- Surv(ty[,1], ty[,2])
gene_ln=colnames(tx)[-17815]
p_n=c()
for (i in 1:length(gene_in$V1)){

	g_1=coxph(my.surv~age_at_initial_pathologic_diagnosis+as.matrix(tx[,which(colnames(tx)==gene_in$V1[i])]),data=dc)
	p_n[i]=summary(g_1)$coefficients[2,5]

}
write.table(data.frame(gene_in$V1,(p_n)), "D:/TCGA/RM_outlier/single_gene_217.txt", quote=F,row.names=F,sep="\t")


module1_gene=read.table("D:/TCGA/RM_outlier/module1.txt",header=F,stringsAsFactors=F)
module2_gene=read.table("D:/TCGA/RM_outlier/module2.txt",header=F,stringsAsFactors=F)
module3_gene=read.table("D:/TCGA/RM_outlier/module3.txt",header=F,stringsAsFactors=F)

m_gene=c(module1_gene$V1,module2_gene$V1,module3_gene$V1)

m_gene_exp=matrix(0,nrow=length(m_gene),ncol=512)
for(i in 1:length(m_gene)){

	m_gene_exp[i,]=tx[,which(colnames(tx)==m_gene[i])]

}



library(ggplot2)
library(reshape2)
qplot(x=Var1,y=Var2,data=melt(cor(t(m_gene_exp))),fill=value,geom="tile")

m_gene_hm=melt(abs(cor(t(m_gene_exp))))
colnames(m_gene_hm)=c("Gene1","Gene2","Correlation")

ggsave("D:/TCGA/RM_outlier/heatmap.pdf",dpi=600)
ggplot(m_gene_hm,aes(Gene1, Gene2)) + geom_tile(aes(fill = Correlation), color="white") + scale_fill_gradient(low="green",high="red")

pdf("D:/TCGA/RM_outlier/heatmap.pdf",dpi=600)

heatmap.2(abs(cor(t(m_gene_exp))), dendrogram="none",na.rm=TRUE,Rowv=F,symm=T)
