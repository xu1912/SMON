##Read the pre-processed data
comp_d=(as.matrix(read.table("D:/TCGA/mRNA_ans/complete_data.txt",sep=",",header=T,row.names=1)))

pd=read.table("D:/TCGA/mRNA_ans/phe.csv",sep=",",header=T,stringsAsFactor=F,na.strings="NA")

pd[which(is.na(pd$days_to_death)),9]=pd[which(is.na(pd$days_to_death)),8]
pd$status=0
pd[which(pd$days_to_last_followup>=pd$days_to_death & pd$vital_status=="Dead"),11]=1

library(survival)

dc=pd[complete.cases(pd[,c(-3,-5)]),c(-3,-5)]

cdd=data.frame(cbind(row.names(comp_d),comp_d),stringsAsFactors=F)
cdd[,1]=gsub("\\.","-",cdd[,1])
colnames(cdd)[13112]="Rgr_small"

library(sqldf)

tdd=sqldf("select a.* from cdd a, dc b where a.V1=b.bcr_patient_barcode")

ty=dc[,c(6,9)]
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

library(qvalue)

q_n=qvalue(p_n)




library(glmnet)

y=ty[which(ty[,1]>0),]
y=as.matrix(y)
x=tx[which(ty[,1]>0),]

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


