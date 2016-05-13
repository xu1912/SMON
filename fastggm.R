install.packages(pkgs = "C:/Users/biostatadmin/Downloads/FastGGM.tar.gz", repos = NULL, type = "source")
library(Rcpp)
Sys.setenv("PKG_CPPFLAGS"="-std=c++0x")
library(FastGGM)


d=read.table("D:/TCGA/RM_outlier/mrna_exp_rm_outlier.txt",row.names=1,header=T,sep=",")
d_mi=read.table("D:/TCGA/RM_outlier/mirna_exp_rm_outlier.txt",row.names=1,header=T,sep=",")
d_me=read.table("D:/TCGA/RM_outlier/methy_combined_Level_3_JHU_USC__HumanMethylation27_rm_outlier.txt",row.names=1,header=T,sep=",")

genl=read.table("D:/TCGA/RM_outlier/module1.txt")
mirnal=read.table("D:/TCGA/RM_outlier/mirna_in_module1.txt")

module=read.table("1_module.txt",header=F,sep="\t",stringsAsFactors=F)

m3_gene=m3[which(m3[,2]=="gene"),1]
m3_mirna=m3[which(m3[,2]=="mirna"),1]
m3_methy=m3[which(m3[,2]=="methy"),1]

d_mi_ins=d[match(rownames(d_mi),rownames(d)),]

dg=matrix(0, nrow=length(d_mi_ins[,1]), ncol=length(genl$V1))
dg_mi=matrix(0, nrow=length(d_mi_ins[,1]), ncol=length(mirnal$V1))

for(i in 1:length(genl$V1)){
	dg[,i]=d_mi_ins[,which(colnames(d_mi_ins)==genl$V1[i])]
}

colnames(dg)=genl$V1

for(i in 1:length(mirnal$V1)){
	dg_mi[,i]=d_mi[,which(colnames(d_mi)==mirnal$V1[i])]
}

colnames(dg_mi)=mirnal$V1
cor(dg,dg_mi)


outlist1 <- FastGGM(dg)

pq=qvalue((outlist1$p_partialCor))

nnc=c()
nnp=c()
n1=c()
n2=c()
for (i in 1:(length(outlist1$p_partialCor[1,])-1)){

	for (j in (i+1):length(outlist1$p_partialCor[1,])){
			n1=append(n1,colnames(dg)[i])
			n2=append(n2,colnames(dg)[j])
			nnc=append(nnc,cor(dg[,i],dg[,j]))
			nnp=append(nnp,outlist1$p_partialCor[i,j])
	}

}
module_eg=data.frame(n1,n2,nnc,nnp)
pq=qvalue((nnp))

module_eg[which(p.adjust(nnp,method="BH")<0.05),1:3]
write.table(module_eg[which(p.adjust(nnp,method="BH")<0.05),1:3],"D:/TCGA/RM_outlier/module3_parcor_qvalue.txt",quote=F,sep="\t",row.names=F)


nnc=c()
n1=c()
n2=c()
for (i in 1:(length(outlist1$p_partialCor[1,])-1)){

	for (j in (i+1):length(outlist1$p_partialCor[1,])){
		if (pq$qvalues[i,j]<0.05){
			n1=append(n1,colnames(dg)[i])
			n2=append(n2,colnames(dg)[j])
			nnc=append(nnc,cor(dg[,i],dg[,j]))
		}
	}

}

module1_eg=data.frame(n1,n2,nnc)
write.table(module1_eg,"D:/TCGA/RM_outlier/module1_parcor_qvalue.txt",quote=F,sep="\t",row.names=F)


nnc=c()
n1=c()
n2=c()
for (i in 1:(length(outlist1$p_partialCor[1,])-1)){

	for (j in (i+1):length(outlist1$p_partialCor[1,])){
		if (outlist1$p_partialCor[i,j]<0.05){
			n1=append(n1,colnames(dg)[i])
			n2=append(n2,colnames(dg)[j])
			nnc=append(nnc,cor(dg[,i],dg[,j]))
		}
	}

}

module3_eg=data.frame(n1,n2,nnc)
write.table(module3_eg,"D:/TCGA/RM_outlier/module3_parcor.txt",quote=F,sep="\t",row.names=F)
