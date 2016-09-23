##Read the pre-processed data
library(survival)
library(sqldf)
library(qvalue)
library(glmnet)

dc=read.table("D:/TCGA/mRNA_ans/phe.csv",sep=",",header=T,stringsAsFactor=F,na.strings="NA")
cdd=read.table("D:/TCGA/mRNA_ans/complete_data.txt",sep=",",header=T,row.names=1,stringsAsFactors=F)

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

write.table(wgcna_exp,"wgcna_exp_input.txt",quote=F,sep="\t",col.names=T,row.names=F)

idx=numeric()
j=0
for(i in 1:length(r_beta_m[,1])){

	if(length(which(r_beta_m[i,]!=0))>0) {
		j=j+1
		idx[j]=i
	}

}


