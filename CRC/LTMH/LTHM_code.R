###### CRC ######
### 1.FDR
library(kinship2)

setwd("/home2/wjkim/project/CRC/20171020/1.model_based/2.mendel/0.data")
fin.dat <- read.table("reduced_FDR_mendel.in",head=F,sep=",",stringsAsFactor=F)
fin.dat[fin.dat[,3]=='',3] <- 0
fin.dat[fin.dat[,4]=='',4] <- 0
colnames(fin.dat)[c(1:5,14,17)] <- c('FID','IID','PID','MID','SEX','DummyAge1','CRC')
fin.dat$SEX <- ifelse(fin.dat$SEX=='M',0,1)

fin.dat$std_age <- (fin.dat$DummyAge1-mean(fin.dat$DummyAge1))/sd(fin.dat$DummyAge1)
fin.dat$std_sex <- (fin.dat$SEX-mean(fin.dat$SEX))/sd(fin.dat$SEX)

model <- CRC~std_age+std_sex-1
total_ped <- with(fin.dat,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
V <- 2*as.matrix(kinship(total_ped))q
famid <- as.character(fin.dat$FID)
prev=0.002326
max.iter=100;max.sub.iter=50
proband <- ifelse(fin.dat$V8=='proband',1,0)

source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
out <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/2.CRC/1.FDR/LTMH_FDR_agesex.txt"
write.table(data.frame('beta_stdage','beta_stdsex','h2','n.iter'),out,row.names=F,col.names=F,quote=F)
res <- LTMH.asc(model=model,
                init_beta=c(0.281342275,0.07401957413),
                init_h2=0.3512363539,
                V=V,
                famid=famid,
                prev=0.002326,
                dataset=fin.dat,
                n.cores=48,
                proband=proband,
                SAVE=T,
                out=out)




#### GCTA
library(gap)
library(kinship2)
setwd("/home2/wjkim/project/CRC/20171020/1.model_based/2.mendel/0.data")
fin.dat <- read.table("reduced_FDR_mendel.in",head=F,sep=",",stringsAsFactor=F)
fin.dat[fin.dat[,3]=='',3] <- 0
fin.dat[fin.dat[,4]=='',4] <- 0
colnames(fin.dat)[c(1:5,14,17)] <- c('FID','IID','PID','MID','SEX','DummyAge1','CRC')
fin.dat$SEX <- ifelse(fin.dat$SEX=='M',0,1)

fin.dat$std_age <- (fin.dat$DummyAge1-mean(fin.dat$DummyAge1))/sd(fin.dat$DummyAge1)
fin.dat$std_sex <- (fin.dat$SEX-mean(fin.dat$SEX))/sd(fin.dat$SEX)
dataset <- fin.dat

model <- CRC~std_age-1
total_ped <- with(fin.dat,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
V <- 2*as.matrix(kinship(total_ped))
grm <- V[upper.tri(V,diag=T)]
i='CRC'
prev=0.002326
WriteGRMBin(paste0('tmp_',i,'_grm'),grm,nrow(dataset)*(nrow(dataset)+1)/2,cbind(dataset[,'FID'],dataset[,'IID']))
write.table(dataset[,c('FID','IID','CRC')],paste0('tmp_',i,'_phen.txt'),row.names=F,col.names=F,quote=F)
write.table(dataset[,c('FID','IID','std_age')],paste0('tmp_',i,'_qcovar.txt'),row.names=F,col.names=F,quote=F)

system(paste0('gcta64 --grm tmp_',i,'_grm --reml --pheno tmp_',i,'_phen.txt --qcovar tmp_',i,'_qcovar.txt --prevalence ',prev,' --out tmp_',i,'_output'))
res <- system(paste0('grep V\\(G\\)\\/Vp_L tmp_',i,'_output.hsq'),intern=T)
res.I <- cbind(i,matrix(strsplit(res,'\t')[[1]],nrow=1)[,-1,drop=F])
colnames(res.I) <- c('Obs','h2','SE_h2')
print(res.I)
out <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/2.CRC/1.FDR/GCTA_FDR.txt"
write.table(res.I,out,col.names=F,row.names=F,quote=F,append=TRUE)

system(paste0('rm tmp_',i,'_*'))




#### CEST

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
X <- model.matrix(model,fin.dat)				
Y <- fin.dat[,'CRC']

source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
res <- Testing.h2.new(famid=famid,V=V,init_beta=0,prev=0.002326,X=X,Y=Y,n.cores=40,proband=proband)
out <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/2.CRC/1.FDR/CEST_FDR.txt"
write.table(res,out,row.names=F,quote=F)



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
