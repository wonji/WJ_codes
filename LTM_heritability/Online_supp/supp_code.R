LTMH.asc <- function(model,init_beta,init_h2,V,famid,prev,dataset,max.iter=100,max.sub.iter=50,n.cores=1,proband,SAVE=F,out=NULL)
  
LTMH <- function(model,init_beta,init_h2,V,famid,prev,data,max.iter=100,max.sub.iter=50,n.cores=1)

  
LTMH(model=model,init_beta=1,init_h2=init_h2,V=V,famid=famid,prev=prev,data=data,n.cores=n.cores)

CEST.h2 <- function(model,data,famid,V,init_beta=NULL,prev,n.cores=1,proband=NULL)
    
 
  
  
###### Random sample   
library(kinship2)

prev=0.1
h2=0.2
totalfam=500
n.cores=24

setwd("~/paper/heritability/ML_ver2/variousFam/")
fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)

set.seed(6)
famlist <- unique(fin.dat$FID)
target <- sample(famlist,totalfam)
dataset = fin.dat[fin.dat$FID%in%target,]
dataset <- dataset[,c('FID','IID','PID','MID','SEX','Y','snp')]
total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
V <- 2*as.matrix(kinship(total_ped))
write.table(dataset,'Online_supp/Random_sp.txt',row.names=F,quote=F)
write.table(V,'Online_supp/Random_kinship.txt',row.names=F,quote=F,col.names=F)

famid <- as.character(dataset$FID)
model <- Y~snp-1

res <- LTMH(model=model,init_h2=h2,V=V,famid=famid,prev=prev,data=dataset,n.cores=n.cores)
CEST.h2(model=model,data=dataset,famid=famid,V=V,init_beta=res$beta_std,prev=prev,n.cores=n.cores)


###### Ascertained sample
library(kinship2)

prev=0.1
h2=0.4
totalfam=500
n.cores=24

setwd("~/paper/heritability/ML_ver2/variousFam/")
fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
fin.dat <- fin.dat[fin.dat$FID%in%fin.dat$FID[fin.dat$ind==1 & fin.dat$Y==1],,drop=F]
PB <- 'ind'

set.seed(1)
famlist <- unique(fin.dat$FID)
target <- sample(famlist,totalfam)
dataset = fin.dat[fin.dat$FID%in%target,]
dataset <- dataset[,c('FID','IID','PID','MID','SEX','Y','snp','ind')]
proband <- dataset[,PB]
total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
V <- 2*as.matrix(kinship(total_ped))
famid <- as.character(dataset$FID)
model <- Y~snp-1
n.cores=24
write.table(dataset,'Online_supp/Ascertained_sp.txt',row.names=F,quote=F)
write.table(V,'Online_supp/Ascertained_kinship.txt',row.names=F,quote=F,col.names=F)


setwd("~/paper/heritability/ML_ver2/variousFam/Online_supp")
source("LTMH_Sourcecode.R")
dataset <- read.table("Random_sp.txt",head=T,stringsAsFactor=F)
V <- as.matrix(read.table("Random_kinship.txt",head=F))

LTMH(model=Y~snp-1,
     data=dataset,
     init_h2=0.2,
     V=V,
     famid=dataset$FID,
     prev=0.1,
     n.cores=24)

res <- LTMH.asc(model=model,init_h2=h2,V=V,famid=famid,prev=prev,data=dataset,n.cores=n.cores,proband=proband)
CEST.h2(model=model,data=dataset,famid=famid,V=V,init_beta=res$beta_std,prev=prev,n.cores=n.cores,proband=proband)

source("LTMH_Sourcecode.R")
dataset <- read.table("Ascertained_sp.txt",head=T,stringsAsFactor=F)
V <- as.matrix(read.table("Random_kinship.txt",head=F))

LTMH.asc(model=Y~snp-1,
         data=dataset,
         init_h2=0.4,
         V=V,
         famid=dataset$FID,
         prev=0.1,
         proband=dataset$ind,
         n.cores=24)

CEST.h2(model=Y~snp-1,
        data=dataset,
        init_beta=NULL,
        V=V,
        famid=dataset$FID,
        prev=0.1,
        proband=dataset$ind,
        n.cores=24)