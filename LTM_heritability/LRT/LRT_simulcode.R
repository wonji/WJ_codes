############## Simulation code
###### generate directories & use CEST dataset
for(prev in c(0.1,0.2)){
  for(h2 in c(0, 0.2, 0.4)){
    print(paste0('prevalence:',prev,', heritability:',h2))
    setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/LRT/1.h2")
    system(paste0("mkdir prev_",prev,"_h2_",h2))
  }
}

###### size : LRT_h2_500_all.txt
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/LRT/LRT.R')
n.cores=20
for(h2 in c(0)){
  for(prev in c(0.1)){
    for(totalfam in c(500)){
      print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
      library(kinship2)
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
      fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
      init_beta=matrix(0,1,1)
      model <- Y~std_snp-1
      out <- paste0("prev_",prev,"_h2_",h2,"/LRT_h2_",totalfam,".txt")
      
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/LRT/1.h2")
      write.table(data.frame('MLE_beta','MLE_h2','Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)
      
      LRT_h2 <- sapply(1:1000,LRT.DoTest.h2,fin.dat=fin.dat,totalfam=totalfam,model=model,init_beta=init_beta,init_h2=h2,prev=prev,h2=h2,n.cores=n.cores,out=out)
    }
  }
}

source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/LRT/LRT.R')
n.cores=24
for(h2 in c(0)){
  for(prev in c(0.2)){
    for(totalfam in c(500)){
      print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
      library(kinship2)
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
      fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
      init_beta=matrix(0,1,1)
      model <- Y~std_snp-1
      out <- paste0("prev_",prev,"_h2_",h2,"/LRT_h2_",totalfam,".txt")
      
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/LRT/1.h2")
      write.table(data.frame('MLE_beta','MLE_h2','Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)
      
      LRT_h2 <- sapply(1:1000,LRT.DoTest.h2,fin.dat=fin.dat,totalfam=totalfam,model=model,init_beta=init_beta,init_h2=h2,prev=prev,h2=h2,n.cores=n.cores,out=out)
    }
  }
}
