##############################################################
################ Software Comparison : h2 ####################
##############################################################

############## GCTA ###############
#### Make directory
for(prev in c(0.05,0.1,0.2)){
  for(h2 in c(0.05,0.2, 0.4)){
    setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/otherSW/1.GCTA/")
    system(paste0("mkdir prev_",prev,"_h2_",h2))
  }
}


#### Do GCTA
source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
num_snp=1;n.cores=1
for(prev in c(0.05,0.1,0.2)){
  for(h2 in c(0.05,0.2,0.4)){
    for(totalfam in c(500)){
      print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
      library(kinship2)
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/")
      fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
      fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
      out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/otherSW/1.GCTA/prev_",prev,"_h2_",h2,"/Esth2_",totalfam,".txt")
      write.table(data.frame('Obs','h2','SE_h2'),out,col.names=F,row.names=F,quote=F)
      working_dir <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/otherSW/1.GCTA/prev_",prev,"_h2_",h2)
      Esth2 <- sapply(1:300,getEsth2.GCTA,fin.dat=fin.dat,totalfam=totalfam,prev=prev,res_var='Y',exp_var='std_snp',out=out,working_dir=working_dir,seed=T)
    }
  }
}


################### Ascertained
#### Convergence of heritability
## family size : 500
source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
n.cores=1
for(prev in c(0.05,0.1,0.2)){
  for(h2 in c(0.05,0.2,0.4)){
    for(totalfam in c(500)){
      print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
      library(kinship2)
      setwd("~/paper/heritability/ML_ver2/variousFam/")
      fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
      fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
      fin.dat <- fin.dat[fin.dat$FID%in%fin.dat$FID[fin.dat$ind==1 & fin.dat$Y==1],,drop=F]

      model <- Y~std_snp-1
      out <- paste0("~/paper/heritability/ML_ver2/variousFam/otherSW/1.GCTA/prev_",prev,"_h2_",h2,"/Esth2_",totalfam,"_asc.txt")
      write.table(data.frame('obs','beta_std_snp','h2','n_iteration'),out,col.names=F,row.names=F,quote=F)
      working_dir <- paste0("~/paper/heritability/ML_ver2/variousFam/otherSW/1.GCTA/prev_",prev,"_h2_",h2)
      Esth2 <- sapply(1:300,getEsth2.GCTA,fin.dat=fin.dat,totalfam=totalfam,prev=prev,res_var='Y',exp_var='std_snp',out=out,working_dir=working_dir,seed=T)
    }
  }
}


