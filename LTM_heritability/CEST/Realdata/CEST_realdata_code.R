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
init_h2 <- 0

# n10 24 cores, chr 13,14
n.cores <- 24 
for(chr in c(13,14)){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  tmp.DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}

# n11 24 cores, chr 16,18
n.cores <- 24 
for(chr in c(16,18)){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  tmp.DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}

# n8 24 cores, chr 9
n.cores <- 24
for(chr in c(9)){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  tmp.DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}

# n15 32 cores, chr 6,7
n.cores <- 32 
for(chr in c(6,7)){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  tmp.DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}

# n9 24 cores, chr 6
n.cores <- 24
for(chr in 6){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  tmp.DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}

# n4 20 cores, chr 7
n.cores <- 20
for(chr in c(7)){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  tmp.DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}

# n14 32 cores, chr 12
n.cores <- 32
for(chr in 12){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}

# n7 24 cores, chr 11
n.cores <- 24 
for(chr in c(11)){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  tmp.DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
}

# n2 20 cores, chr 10
n.cores <- 20
for(chr in c(10)){
  out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/chr",chr,"/GRM_181024.txt")
  tmp.DoCEST.beta.WGS(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx,out)
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
write.table(plotdat,'LAM_plotdata.txt',row.names=F,quote=F)

system('Rscript ~/Output_plot.r LAM_plotdata.txt LAM_GWAS')




################# SNUH ######################
# prevalence : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5911525/ 10.9% 2013
library(kinship2)
fam <- read.table("/data/kare/allmerge/allmerge_withrela.fam",head=F,stringsAsFactor=F)
new.fam <- fam[-grep("^KNIH",fam$V1),]
total_ped <- with(new.fam,pedigree(id=V2,dadid=V3,momid=V4,sex=V5,famid=V1,missid='0'))
VV <- 2*as.matrix(kinship(total_ped))

# Only intercept
ind <- which(new.fam$V6!=-9)
dataset <- new.fam[ind,,drop=F]
dataset$V6 <- dataset$V6-1
V <- VV[ind,ind]
model <- V6 ~ 1
famid <- dataset$V1
PB.candi <- grep('PKS',dataset$V2)
PB.fam <- dataset$V1[PB.candi]
PB <- PB.candi[!duplicated(PB.fam)]
proband <- rep(0,nrow(dataset))
proband[PB] <- 1
out <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/1.SNUH/OnlyIC_181113.txt"

source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
write.table(data.frame('intercept','h2','n_iteration'),out,row.names=F,col.names=F)

res <- LTMH.asc(model=model,
				init_beta=1,
				init_h2=0.26,
				V=V,
				famid=famid,
				prev=0.109,
				dataset=dataset,
				n.cores=24,
				proband=proband,
				SAVE=T,
				out=out)

### CEST
# prevalence : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5911525/ 10.9% 2013
library(kinship2)
fam <- read.table("/data/kare/allmerge/allmerge_withrela.fam",head=F,stringsAsFactor=F)
new.fam <- fam[-grep("^KNIH",fam$V1),]
total_ped <- with(new.fam,pedigree(id=V2,dadid=V3,momid=V4,sex=V5,famid=V1,missid='0'))
VV <- 2*as.matrix(kinship(total_ped))

# Only intercept
ind <- which(new.fam$V6!=-9)
dataset <- new.fam[ind,,drop=F]
dataset$V6 <- dataset$V6-1
V <- VV[ind,ind]
model <- V6 ~ 1
famid <- dataset$V1
PB.candi <- grep('PKS',dataset$V2)
PB.fam <- dataset$V1[PB.candi]
PB <- PB.candi[!duplicated(PB.fam)]
proband <- rep(0,nrow(dataset))
proband[PB] <- 1
out <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/1.SNUH/CEST.txt"

init_beta=0.755458
prev <- 0.109
Y <- dataset[,'V6']
X <- model.matrix(model,dataset)				
n.cores=24
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
res <- Testing.h2.new(famid,V,init_beta,prev,X,Y,n.cores,proband)
write.table(res,out,row.names=F,quote=F)
				
				

##### With age
library(kinship2)
fam <- read.table("/data/kare/allmerge/allmerge_withrela.fam",head=F,stringsAsFactor=F)
new.fam <- fam[-grep("^KNIH",fam$V1),]
pheno <- read.csv("/home2/wjkim/project/kare+snu/snu/snu_pheno.csv",head=T,stringsAsFactor=F)
rela <- read.table("/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/1.SNUH/SNUH_rela_age.txt",head=T,stringsAsFactor=F)

case <- pheno[pheno$GWAS.ID!="",c('GWAS.ID','AGE')]
rela <- rela[,c('IID','AGE')]
colnames(case) <- c('IID','AGE')
age.info <- rbind(case,rela)

mer <- merge(new.fam,age.info,by.x='V2',by.y='IID',sort=F,all=F)
total_ped <- with(mer,pedigree(id=V2,dadid=V3,momid=V4,sex=V5,famid=V1,missid='0'))
VV <- 2*as.matrix(kinship(total_ped))

ind <- which(mer$V6!=-9 & !is.na(mer$AGE))
dataset <- mer[ind,,drop=F]
dataset$V6 <- dataset$V6-1
dataset$std_AGE <- (dataset$AGE-mean(dataset$AGE))/sd(dataset$AGE)
V <- VV[ind,ind]
model <- V6 ~ std_AGE-1
famid <- dataset$V1
PB.candi <- grep('PKS',dataset$V2)
PB.fam <- dataset$V1[PB.candi]
PB <- PB.candi[!duplicated(PB.fam)]
proband <- rep(0,nrow(dataset))
proband[PB] <- 1
out <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/1.SNUH/CEST_age.txt"

init_beta=1
prev <- 0.109
Y <- dataset[,'V6']
X <- model.matrix(model,dataset)				
n.cores=24
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
res <- Testing.h2.new(famid,V,init_beta,prev,X,Y,n.cores,proband)
write.table(res,out,row.names=F,quote=F)
				
### LTMH
source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
out <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/1.SNUH/LTMH_age.txt"
res <- LTMH.asc(model=model,
				init_beta=1,
				init_h2=0.26,
				V=V,
				famid=famid,
				prev=0.109,
				dataset=dataset,
				n.cores=24,
				proband=proband,
				SAVE=T,
				out=out)

## beta hat for unstandardized age				
setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/1.SNUH/")				
res <- read.table("LTMH_age.txt",head=F,stringsAsFactor=F)

fam <- read.table("/data/kare/allmerge/allmerge_withrela.fam",head=F,stringsAsFactor=F)
new.fam <- fam[-grep("^KNIH",fam$V1),]
pheno <- read.csv("/home2/wjkim/project/kare+snu/snu/snu_pheno.csv",head=T,stringsAsFactor=F)
rela <- read.table("/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/1.SNUH/SNUH_rela_age.txt",head=T,stringsAsFactor=F)

case <- pheno[pheno$GWAS.ID!="",c('GWAS.ID','AGE')]
rela <- rela[,c('IID','AGE')]
colnames(case) <- c('IID','AGE')
age.info <- rbind(case,rela)

mer <- merge(new.fam,age.info,by.x='V2',by.y='IID',sort=F,all=F)
ind <- which(mer$V6!=-9 & !is.na(mer$AGE))
dataset <- mer[ind,,drop=F]
beta.age <- res[nrow(res),1]/sd(dataset$AGE)



###### CRC ######
### 1.FDR
library(kinship2)

setwd("/home2/wjkim/project/CRC/20171020/1.model_based/2.mendel/0.data")
fin.dat <- read.table("reduced_FDR_mendel.in",head=F,sep=",",stringsAsFactor=F)
fin.dat[fin.dat[,3]=='',3] <- 0
fin.dat[fin.dat[,4]=='',4] <- 0
colnames(fin.dat)[c(1:5,14,17)] <- c('FID','IID','PID','MID','SEX','DummyAge1','CRC')
fin.dat$SEX <- ifelse(fin.dat$SEX=='M',1,2)

fin.dat$std_age <- (fin.dat$DummyAge1-mean(fin.dat$DummyAge1))/sd(fin.dat$DummyAge1)
fin.dat$std_sex <- (fin.dat$SEX-mean(fin.dat$SEX))/sd(fin.dat$SEX)

model <- CRC~std_age-1
total_ped <- with(fin.dat,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
V <- 2*as.matrix(kinship(total_ped))
famid <- as.character(fin.dat$FID)
prev=0.002326
max.iter=100;max.sub.iter=50
proband <- ifelse(fin.dat$V8=='proband',1,0)

source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
out <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/2.CRC/1.FDR/LTMH_FDR.txt"
write.table(data.frame('beta_stdage','h2','n.iter'),out,row.names=F,col.names=F,quote=F)
res <- LTMH.asc(model=model,
				init_beta=0.2,
				init_h2=0.12,
				V=V,
				famid=famid,
				prev=0.002326,
				dataset=fin.dat,
				n.cores=32,
				proband=proband,
				SAVE=T,
				out=out)
