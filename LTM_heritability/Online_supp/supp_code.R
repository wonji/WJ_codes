LTMH.asc <- function(model,init_beta,init_h2,V,famid,prev,dataset,max.iter=100,max.sub.iter=50,n.cores=1,proband,SAVE=F,out=NULL)
  
LTMH <- function(model,init_beta,init_h2,V,famid,prev,data,max.iter=100,max.sub.iter=50,n.cores=1)

  
LTMH(model=model,init_beta=1,init_h2=init_h2,V=V,famid=famid,prev=prev,data=data,n.cores=n.cores)

CEST.h2 <- function(model,data,famid,V,init_beta=NULL,prev,n.cores=1,proband=NULL)
    
 
  
  
   
library(kinship2)

prev=0.1
h2=0.2
totalfam=500
n.cores=24

setwd("~/paper/heritability/ML_ver2/variousFam/")
fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)

set.seed(6)
famlist <- unique(fin.dat$FID)
target <- sample(famlist,totalfam)
dataset = fin.dat[fin.dat$FID%in%target,]
total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
V <- 2*as.matrix(kinship(total_ped))
famid <- as.character(dataset$FID)
model <- Y~snp-1
init_h2 <- 0.2
prev=0.1
n.cores=24

res <- LTMH(model=model,init_h2=init_h2,V=V,famid=famid,prev=prev,data=dataset,n.cores=n.cores)






fin.dat <- fin.dat[fin.dat$FID%in%fin.dat$FID[fin.dat$ind==1 & fin.dat$Y==1],,drop=F]
PB <- 'ind'


