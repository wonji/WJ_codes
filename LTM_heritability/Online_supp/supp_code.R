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
total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
V <- 2*as.matrix(kinship(total_ped))
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
proband <- dataset[,PB]
total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
V <- 2*as.matrix(kinship(total_ped))
famid <- as.character(dataset$FID)
model <- Y~snp-1
n.cores=24

res <- LTMH.asc(model=model,init_h2=h2,V=V,famid=famid,prev=prev,data=dataset,n.cores=n.cores,proband=proband)
CEST.h2(model=model,data=dataset,famid=famid,V=V,init_beta=res$beta_std,prev=prev,n.cores=n.cores,proband=proband)

