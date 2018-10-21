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
	  init_beta=matrix(0.3,1,1)
      model <- Y~std_snp-1
      
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      write.table(data.frame('Chisq','Pvalue'),paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,".txt"),col.names=F,row.names=F,quote=F)
      
      CEST_h2 <- sapply(1:2000,DoTest.h2,fin.dat=fin.dat,totalfam=totalfam,model=model,init_beta=init_beta,prev=prev,h2=h2,n.cores=n.cores)
    }
  }
}


