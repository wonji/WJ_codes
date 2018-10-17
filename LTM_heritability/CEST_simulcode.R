############### CEST simulation code #################

###### Test h2 
#### Generate simulation data
# prev : 0.1,0.2 | h2 : 0.2,0.4
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/gen_simuldata.r')
for(prev in c(0.1,0.2)){
  for(h2 in c(0, 0.2, 0.4)){
    print(paste0('prevalence:',prev,', heritability:',h2))
    setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
    system(paste0("mkdir prev_",prev,"_h2_",h2))
    dataset = genNucFam(totalfam=10000,MAF=0.2,h2,ha2=0.05,prev,num_snp=1)
    write.table(dataset,paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),row.names=F,quote=F)
  }
}



source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
n.cores=32
for(prev in c(0.1,0.2)){
  for(h2 in c(0, 0.2, 0.4)){
    for(totalfam in c(500)){
      print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
      library(kinship2)
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
      fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
      
      # dataset : data.frame
      # dataset = fin.dat[fin.dat$FID%in%sample(unique(fin.dat$FID),500),]
      X <- as.matrix(dataset[,'std_snp',drop=F])
      Y <- as.matrix(dataset[,'Y',drop=F])
      init_beta <- matrix(2.5,1,1)
      famid <- dataset$FID
      total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
      V <- 2*as.matrix(kinship(total_ped))
      
      model <- Y~std_snp-1
      
      # families containing at least one affected member
      all.famlist <- unique(fin.dat$FID)
      num.aff <- sapply(all.famlist,function(i) sum(fin.dat$Y[fin.dat$FID==i]))
      fin.dat <- fin.dat[fin.dat$FID%in%all.famlist[num.aff>0],]
      
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/ascertained")
      write.table(data.frame('obs','beta_std_snp','h2','n_iteration'),paste0("prev_",prev,"_h2_",h2,"/Esth2_",totalfam,".txt"),col.names=F,row.names=F,quote=F)
      
      Esth2 <- sapply(1:300,getEsth2,fin.dat=fin.dat,init_beta=1,init_h2=h2,totalfam=totalfam,assumed_prev=prev,model=model,n.cores=n.cores)
    }
  }
}

