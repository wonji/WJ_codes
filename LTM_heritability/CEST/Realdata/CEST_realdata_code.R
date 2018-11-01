###################################################
########### CEST real data analysis ###############
###################################################

##################
#### LAM : Do CEST for two significant SNPs (n4)
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


#### Calculate MAC data for two significant SNPs
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


######### Pigmentation SNPs
#setwd("/home2/wjkim/project/LAM/final_analysis/0.data/")
#system("plink --bfile final_clean_withage --snps rs1129038,rs12913832 --recodeA --out /home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/0.data/pigSNPs")

#### Do CEST (Corr = 1)
##  Sig SNPs
setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/")
raw <- read.table("0.data/pigSNPs.raw",head=T,stringsAsFactor=F)

setwd("/home2/wjkim/project/LAM/final_analysis/")
phe <- read.table("0.data/phenofile_withage.txt",head=T,stringsAsFactor=F)
rm.list <- read.table("3.clogit/1.revision/0.data/outliers.txt",head=F,stringsAsFactor=F)
phe$pair_age[which(phe$pair_age%in%phe[phe$IID%in%rm.list$V2,'pair_age'])] <- NA
valid.idx <- which(!is.na(phe$pair_age))
new.phe <- phe[valid.idx,]

setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/0.data")
VV <- as.matrix(read.table("Corr_1.cor",head=T,stringsAsFactor=F))
out <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/1.sigSNPs/pigSNPs_Corr_1_181029.txt"
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
  n.cores <- 1
  prev <- 0.00001
  Testing.res <- Testing.beta(init_beta,init_h2,V,famid,prev,Y,X,full.X,test.beta.idx,max.iter=100,max.sub.iter=50,n.cores)
  res <- cbind(SNP=gsub('_[ACGT]','',colnames(raw)[ii+6]),Testing.res)
  write.table(res,out,col.names=F,row.names=F,quote=F,append=TRUE)
}

for(jj in 1:2) DoCEST.beta(jj)



#### GRM
setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/")
raw <- read.table("0.data/pigSNPs.raw",head=T,stringsAsFactor=F)

setwd("/home2/wjkim/project/LAM/final_analysis/")
phe <- read.table("0.data/phenofile_withage.txt",head=T,stringsAsFactor=F)
rm.list <- read.table("3.clogit/1.revision/0.data/outliers.txt",head=F,stringsAsFactor=F)
phe$pair_age[which(phe$pair_age%in%phe[phe$IID%in%rm.list$V2,'pair_age'])] <- NA
valid.idx <- which(!is.na(phe$pair_age))
new.phe <- phe[valid.idx,]

setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/0.data")
VV <- as.matrix(read.table("GRM.empi.med.cor",head=T,stringsAsFactor=F))
out <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/1.sigSNPs/pigSNPs_GRM_181029.txt"
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
  n.cores <- 1
  prev <- 0.00001
  Testing.res <- Testing.beta(init_beta,init_h2,V,famid,prev,Y,X,full.X,test.beta.idx,max.iter=100,max.sub.iter=50,n.cores)
  res <- cbind(SNP=gsub('_[ACGT]','',colnames(raw)[ii+6]),Testing.res)
  write.table(res,out,col.names=F,row.names=F,quote=F,append=TRUE)
}

for(jj in 1:2) DoCEST.beta(jj)





##################
#### LAM : Do CEST for whole chromosome
##################

#### 1. Make a CEST function for each chromosome
## Please find the function in CEST.R ('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')


#### 2. Link raw files to the working directory
IN <- "/home2/wjkim/project/LAM/final_analysis/3.clogit/"
OUT <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/"
setwd(OUT)
for(chr in 1:22){
  setwd(OUT)
  system(paste0('mkdir chr',chr))
  setwd(paste0(OUT,"chr",chr))
  system(paste0("ln -s ",IN,"chr",chr,"/lam_chr",chr,".raw chr",chr,".raw"))
}    
  

#### 3. Do CEST
## 3.1 GRM
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')

model <- pheno01 ~ 1
FID <- 'pair_age'
IID <- 'IID'
working_dir <- '/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/'
# pheno file
phe <- read.table("/home2/wjkim/project/LAM/final_analysis/0.data/phenofile_withage.txt",head=T,stringsAsFactor=F)
rm.list <- read.table("/home2/wjkim/project/LAM/final_analysis/3.clogit/1.revision/0.data/outliers.txt",head=F,stringsAsFactor=F)
phe$pair_age[which(phe$pair_age%in%phe[phe$IID%in%rm.list$V2,'pair_age'])] <- NA
valid.idx <- which(!is.na(phe$pair_age))

# relationship matrix (GRM)
setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/0.data")
VV <- as.matrix(read.table("GRM.empi.med.cor",head=T,stringsAsFactor=F))
V.name <- phe[valid.idx,'IID']
prev <- 0.00001
init_beta <- matrix(0,1,1)
init_h2 <- 0.1

# n5 24 cores, chr 1
n.cores <- 24 
for(chr in c(1)){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  tmp.DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}

# n14 32 cores, chr 7,18
n.cores <- 32 
for(chr in c(7,18)){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  tmp.DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}

# n6 24 cores, chr 4,10,14
n.cores <- 24 
for(chr in c(4,10,14)){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  tmp.DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}

# n10 20 cores, chr 2
n.cores <- 20
for(chr in c(2)){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  tmp.DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}

# n2 20 cores, chr 17,19~22
n.cores <- 20
for(chr in c(17,19:22)){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}

# n11 20 cores, chr 5,6
n.cores <- 20
for(chr in c(5,6)){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  if(chr==5) tmp.DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
  if(chr==6) DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}

# n12 50 cores, chr 11~13
n.cores <- 50
for(chr in c(11:13)){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  if(chr==11) tmp.DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
  if(chr!=11) DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}

# n8 24 cores, chr 3,9
n.cores <- 24
for(chr in c(3,9)){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  if(chr==3) tmp.DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
  if(chr==9) DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}

# n9 24 cores, chr 8,15
n.cores <- 24
for(chr in c(8,15)){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  if(chr==8) tmp.DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
  if(chr==15) DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}


setwd('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/')
a <- sapply(1:22,function(chr) return(ifelse('GRM_181024.txt'%in%system(paste0('ls chr',chr),intern=T),chr,'')))
aa <- a[a!='']
dat <- c()
for(i in aa){
	bb <- read.table(paste0('chr',i,'/GRM_181024.txt'),head=T)
	dat <- rbind(dat,bb)
	print(paste0('chr',i,', ',nrow(bb),' SNPs'))
}
plotdat <- dat[,c('CHR','SNP','Pvalue')]
colnames(plotdat) <- c('CHR','GENE','p.val')
write.table(plotdat,'temp_plotdata.txt',row.names=F,quote=F)

system('Rscript ~/Output_plot.r temp_plotdata.txt temp_GWAS')

##################
#### TWIN : Do CEST for testing H0: h2=0 & estimating
##################

#### Load Matching index & grouping
phe <- read.delim("/data/twin/TWIN_1801_phe.txt",head=T,stringsAsFactor=F)

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








