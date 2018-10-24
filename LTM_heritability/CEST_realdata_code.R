###################################################
########### CEST real data analysis ###############
###################################################

##################
#### LAM (n4)
##################

#### Load Matching index & grouping
setwd("/home2/wjkim/project/LAM/final_analysis/")
phe <- read.table("0.data/phenofile_withage.txt",head=T,stringsAsFactor=F)
rm.list <- read.table("3.clogit/1.revision/0.data/outliers.txt",head=F,stringsAsFactor=F)
phe$pair_age[which(phe$pair_age%in%phe[phe$IID%in%rm.list$V2,'pair_age'])] <- NA
valid.idx <- which(!is.na(phe$pair_age))
new.phe <- phe[valid.idx,]

#### Generate relationship matrix 
## 1) Corr = 1
mat <- matrix(0,nrow(new.phe),nrow(new.phe))
for(i in unique(new.phe$pair_age)){
	print(i)
	tmp.idx <- which(new.phe$pair_age==i)
	mat[tmp.idx,tmp.idx] <- 1
}
mat <- rbind(new.phe$IID,mat)
setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/0.data")
write.table(mat,"Corr_1.cor",row.names=F,col.names=F,quote=F)

## 2) GRM
setwd("/home2/wjkim/project/LAM/final_analysis/")
system("wisard --bed 0.data/final_clean_withage --medcor --makecor --out /home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/0.data/GRM --thread 20")


#### Load MAC data
#setwd("/home2/wjkim/project/LAM/final_analysis/0.data/")
#system("plink --bfile final_clean_withage --snps rs4544201,rs2006950 --recodeA --out /home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/0.data/sigSNPs")


#### Do CEST (Corr = 1)
##  Sig SNPs
setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/")
raw <- read.table("0.data/sigSNPs.raw",head=T,stringsAsFactor=F)

setwd("/home2/wjkim/project/LAM/final_analysis/")
phe <- read.table("0.data/phenofile_withage.txt",head=T,stringsAsFactor=F)
rm.list <- read.table("3.clogit/1.revision/0.data/outliers.txt",head=F,stringsAsFactor=F)
phe$pair_age[which(phe$pair_age%in%phe[phe$IID%in%rm.list$V2,'pair_age'])] <- NA
valid.idx <- which(!is.na(phe$pair_age))
new.phe <- phe[valid.idx,]

setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/0.data")
VV <- as.matrix(read.table("Corr_1.cor",head=T,stringsAsFactor=F))
out <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/1.sigSNPs/sigSNPs_Corr_1_181024.txt"
write.table(data.frame('SNP','Score','var_Score_ver1','var_Score_ver2','Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)

source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
DoCEST.beta <- function(ii){
	dataset <- cbind(phe,SNP=raw[,ii+6])[valid.idx,]
	dataset <- dataset[!is.na(dataset$SNP),]
	dataset$stdSNP <- (dataset$SNP-mean(dataset$SNP))/sd(dataset$SNP)
	model <- pheno01~stdSNP
	test.beta <- 'stdSNP'
	
	Y <- dataset[,as.character(model)[2],drop=F]
	full.X <- model.matrix(model,dataset)
	test.beta.idx <- which(colnames(full.X)%in%test.beta)
	X <- full.X[,-test.beta.idx,drop=F]
	
	famid <- dataset$pair_age
	ord.idx <- match(dataset$IID,new.phe$IID)
	V <- VV[ord.idx,ord.idx]
	colnames(V) <- NULL
	
	# Testing
	init_beta <- matrix(0,1,1)
	init_h2 <- 0.05
	n.cores <- 20 
	prev <- 0.00001
	Testing.res <- Testing.beta(init_beta,init_h2,V,famid,prev,Y,X,full.X,test.beta.idx,max.iter=100,max.sub.iter=50,n.cores)
	res <- cbind(SNP=gsub('_[ACGT]','',colnames(raw)[ii+6]),Testing.res)
	write.table(res,out,col.names=F,row.names=F,quote=F,append=TRUE)
}

for(jj in 1:2) DoCEST.beta(jj)



#### GRM
setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/")
raw <- read.table("0.data/sigSNPs.raw",head=T,stringsAsFactor=F)

setwd("/home2/wjkim/project/LAM/final_analysis/")
phe <- read.table("0.data/phenofile_withage.txt",head=T,stringsAsFactor=F)
rm.list <- read.table("3.clogit/1.revision/0.data/outliers.txt",head=F,stringsAsFactor=F)
phe$pair_age[which(phe$pair_age%in%phe[phe$IID%in%rm.list$V2,'pair_age'])] <- NA
valid.idx <- which(!is.na(phe$pair_age))
new.phe <- phe[valid.idx,]

setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/0.data")
VV <- as.matrix(read.table("GRM.empi.med.cor",head=T,stringsAsFactor=F))
out <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/1.sigSNPs/sigSNPs_GRM_181024.txt"
write.table(data.frame('SNP','Score','var_Score_ver1','var_Score_ver2','Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)

source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
DoCEST.beta <- function(ii){
	dataset <- cbind(phe,SNP=raw[,ii+6])[valid.idx,]
	dataset <- dataset[!is.na(dataset$SNP),]
	dataset$stdSNP <- (dataset$SNP-mean(dataset$SNP))/sd(dataset$SNP)
	model <- pheno01~stdSNP
	test.beta <- 'stdSNP'
	
	Y <- dataset[,as.character(model)[2],drop=F]
	full.X <- model.matrix(model,dataset)
	test.beta.idx <- which(colnames(full.X)%in%test.beta)
	X <- full.X[,-test.beta.idx,drop=F]
	
	famid <- dataset$pair_age
	ord.idx <- match(dataset$IID,new.phe$IID)
	V <- VV[ord.idx,ord.idx]
	colnames(V) <- NULL
	
	# Testing
	init_beta <- matrix(0,1,1)
	init_h2 <- 0.05
	n.cores <- 20 
	prev <- 0.00001
	Testing.res <- Testing.beta(init_beta,init_h2,V,famid,prev,Y,X,full.X,test.beta.idx,max.iter=100,max.sub.iter=50,n.cores)
	res <- cbind(SNP=gsub('_[ACGT]','',colnames(raw)[ii+6]),Testing.res)
	write.table(res,out,col.names=F,row.names=F,quote=F,append=TRUE)
}

for(jj in 1:2) DoCEST.beta(jj)






