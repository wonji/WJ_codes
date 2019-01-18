output <- LTMH(model=model,data=dataset,init_beta=NULL,init_h2=0.2,V=V,famid=famid,prev=prev,n.cores=20)

# random sample
prev=0.1
h2=0.2
i=1
seed=T
totalfam=500

library(kinship2)
source('~/LTMH_Sourcecode.R')
setwd("~/paper/heritability/ML_ver2/variousFam/")
fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
model <- Y~snp

print(i)
famlist <- unique(fin.dat$FID)
if(seed) set.seed(i)
target <- sample(famlist,totalfam)
dataset = fin.dat[fin.dat$FID%in%target,]
total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
V <- 2*as.matrix(kinship(total_ped))
famid <- as.character(dataset$FID)
n.cores=20

a <- LRT.beta(model,data=dataset,init_h2=h2,V=V,famid=famid,prev,test.beta='snp',max.iter=100,n.cores=20,proband=NULL)
  





prev=0.1
h2=0.2
i=1
seed=T
totalfam=500

library(kinship2)
source('~/LTMH_Sourcecode.R')
setwd("~/paper/heritability/ML_ver2/variousFam/")
fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
fin.dat <- fin.dat[fin.dat$FID%in%fin.dat$FID[fin.dat$ind==1 & fin.dat$Y==1],,drop=F]
proband <- 'ind'
model <- Y~snp

print(i)
famlist <- unique(fin.dat$FID)
if(seed) set.seed(i)
target <- sample(famlist,totalfam)
dataset = fin.dat[fin.dat$FID%in%target,]
total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
V <- 2*as.matrix(kinship(total_ped))
famid <- as.character(dataset$FID)
n.cores=20

a <- LRT.beta(model,data=dataset,init_h2=h2,V=V,famid=famid,prev,test.beta='snp',max.iter=100,n.cores=20,proband=dataset[,'ind'])
