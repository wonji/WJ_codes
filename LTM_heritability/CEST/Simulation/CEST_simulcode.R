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



###### Ver 1 : Information matrix latter part
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
n.cores=24
for(prev in c(0.1,0.2)){
  for(h2 in c(0, 0.2, 0.4)){
    for(totalfam in c(500)){
      print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
      library(kinship2)
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
      fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
	    # init_beta=matrix(c(1,0.3),2,1)
      # model <- Y~std_snp
	    # out <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,"_tau.txt")
      init_beta=matrix(0,1,1)
      model <- Y~std_snp-1
	    out <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,".txt")
      
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      write.table(data.frame('Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)
      
      CEST_h2 <- sapply(1:2000,DoTest.h2,fin.dat=fin.dat,totalfam=totalfam,model=model,init_beta=init_beta,prev=prev,h2=h2,n.cores=n.cores,out=out)
    }
  }
}

###### n9, Ver 2 : Remove Information matrix latter part
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
n.cores=24
for(prev in c(0.1,0.2)){
  for(h2 in c(0, 0.2, 0.4)){
    for(totalfam in c(500)){
      print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
      library(kinship2)
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
      fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
      # init_beta=matrix(c(1,0.3),2,1)
      # model <- Y~std_snp
      # out <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,"_tau.txt")
      init_beta=matrix(0,1,1)
      model <- Y~std_snp-1
      out <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,"_ver2.txt")
      
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      write.table(data.frame('Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)
      
      CEST_h2 <- sapply(1:2000,DoTest.h2,fin.dat=fin.dat,totalfam=totalfam,model=model,init_beta=init_beta,prev=prev,h2=h2,n.cores=n.cores,out=out,version=2)
    }
  }
}

###### n9, CEST_h2_500_all.txt : Remove Information matrix latter part & consider sign of score (S<0 -> 1/2, S>0 ~ chisq()/2 )
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
n.cores=24
for(h2 in c(0, 0.2, 0.4)){
  for(prev in c(0.1,0.2)){
    for(totalfam in c(500)){
      print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
      library(kinship2)
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
      fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
      init_beta=matrix(0,1,1)
      model <- Y~std_snp-1
      out <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,"_all.txt")
      
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      write.table(data.frame('Score','var_Score_ver1','var_Score_ver2','Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)
      
      CEST_h2 <- sapply(1:2000,DoTest.h2,fin.dat=fin.dat,totalfam=totalfam,model=model,init_beta=init_beta,prev=prev,h2=h2,n.cores=n.cores,out=out)
    }
  }
}

###### n2, CEST_h2_500_tau_all.txt : Remove Information matrix latter part & consider sign of score (S<0 -> 1/2, S>0 ~ chisq()/2 ) & tau
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
n.cores=20
for(h2 in c(0, 0.2, 0.4)){
  for(prev in c(0.1,0.2)){
    for(totalfam in c(500)){
      print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
      library(kinship2)
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
      fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
      init_beta=matrix(c(0,0.2),2,1)
      model <- Y~std_snp
      out <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,"_tau_all.txt")
      
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      write.table(data.frame('Score','var_Score_ver1','var_Score_ver2','Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)
      
      CEST_h2 <- sapply(1:2000,DoTest.h2,fin.dat=fin.dat,totalfam=totalfam,model=model,init_beta=init_beta,prev=prev,h2=h2,n.cores=n.cores,out=out)
    }
  }
}


###### n9, CEST_h2_700_all.txt & CEST_h2_1000_all.txt
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
n.cores=24
for(h2 in c(0, 0.2, 0.4)){
  for(prev in c(0.2,0.1)){
    for(totalfam in c(700,1000)){
      print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
      library(kinship2)
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
      fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
      init_beta=matrix(0,1,1)
      model <- Y~std_snp-1
      out <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,"_all.txt")
      
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      write.table(data.frame('Score','var_Score_ver1','var_Score_ver2','Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)
      
      CEST_h2 <- sapply(1:2000,DoTest.h2,fin.dat=fin.dat,totalfam=totalfam,model=model,init_beta=init_beta,prev=prev,h2=h2,n.cores=n.cores,out=out)
    }
  }
}



###### n9, CEST_h2_500_nobeta.txt : CEST with no beta parameters
## Data generation
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/gen_simuldata.r')
for(prev in c(0.1,0.2)){
  for(h2 in c(0, 0.2, 0.4)){
    num_snp <- 0
    print(paste0('prevalence:',prev,', heritability:',h2))
    setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2/1.NoBeta")
    system(paste0("mkdir prev_",prev,"_h2_",h2))
    dataset = genNucFam(totalfam=10000,MAF=0.2,h2=h2,ha2=0.05,prev=prev,num_snp=num_snp)
    write.table(dataset,paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),row.names=F,quote=F)
  }
}

## Simulation
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
n.cores=24
for(h2 in c(0, 0.2, 0.4)){
  for(prev in c(0.1,0.2)){
    for(totalfam in c(500)){
      print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
      library(kinship2)
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2/1.NoBeta")
      fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
  	  if(h2==0) {
  		  fin.dat$L <- rnorm(nrow(fin.dat))
  	    fin.dat$Y <- ifelse(fin.dat$L<qnorm(prev,lower.tail=F),0,1)
  	  }
	  
      model <- Y~1
	    init_beta <- NULL
      out <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,"_nobeta.txt")
      
      write.table(data.frame('Score','var_Score_ver1','var_Score_ver2','Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)
      
      CEST_h2 <- sapply(1:2000,DoTest.h2,fin.dat=fin.dat,totalfam=totalfam,model=model,init_beta=init_beta,prev=prev,h2=h2,n.cores=n.cores,out=out)
    }
  }
}


##################### Smaller prevalence (Channing server) ############################
###### Test h2 
#### Generate simulation data
# prev : 0.01, 0.05 | h2 : 0, 0.2,0.4
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/gen_simuldata.r')
for(prev in c(0.01,0.05)){
  for(h2 in c(0, 0.2, 0.4)){
    print(paste0('prevalence:',prev,', heritability:',h2))
    setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
    system(paste0("mkdir prev_",prev,"_h2_",h2))
    dataset = genNucFam(totalfam=10000,MAF=0.2,h2,ha2=0.05,prev,num_snp=1)
    write.table(dataset,paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),row.names=F,quote=F)
  }
}


###### n9, CEST_h2_500_all.txt : Remove Information matrix latter part & consider sign of score (S<0 -> 1/2, S>0 ~ chisq()/2 )
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
n.cores=32
if(h2==0) c.prev <- 0.05
if(h2!=0) c.prev <- c(0.01,0.05)

for(h2 in c(0, 0.2, 0.4)){
  for(prev in c.prev){
    for(totalfam in c(500)){
      print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
      library(kinship2)
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
      fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
      init_beta=matrix(0.2,1,1)
      model <- Y~std_snp-1
      out <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,"_all.txt")
      
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      write.table(data.frame('Score','var_Score_ver1','var_Score_ver2','Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)
      
      CEST_h2 <- sapply(1:2000,DoTest.h2,fin.dat=fin.dat,totalfam=totalfam,model=model,init_beta=init_beta,prev=prev,h2=h2,n.cores=n.cores,out=out)
    }
  }
}




###### Test beta
#### Generate simulation data
# prev : 0.1,0.2 | h2 : 0.2,0.4 | size and power
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/gen_simuldata.r')
for(prev in c(0.1,0.2)){
  for(h2 in c(0.2, 0.4)){
    for(num_snp in c(0,1)){
      print(paste0('prevalence:',prev,', heritability:',h2,ifelse(num_snp==0,", size",", power")))
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/2.beta")
      system(paste0("mkdir prev_",prev,"_h2_",h2,"_",ifelse(num_snp==0,"size","power")))
      dataset = genNucFam(totalfam=10000,MAF=0.2,h2,ha2=0.005,prev,num_snp=num_snp)
      write.table(dataset,paste0("prev_",prev,"_h2_",h2,"_",ifelse(num_snp==0,"size","power"),"/dataset.txt"),row.names=F,quote=F)
    }
  }
}


##### version 2 (size : n14, power : n8)
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
n.cores=32
for(type in c('size')){
  for(h2 in c(0.2, 0.4)){
    for(prev in c(0.1,0.2)){
      totalfam <- 500
      print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam,', Test type: ',type))
      library(kinship2)
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/2.beta")
      fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"_",type,"/dataset.txt"),head=T,stringsAsFactor=F)
      fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
      model <- Y~std_snp
      
      test.beta <- 'std_snp'
      init_h2 <- h2
      init_beta=matrix(1,1,1) # initial value for beta which is not for the H0.
      out <- paste0("prev_",prev,"_h2_",h2,"_",type,"/CEST_beta_",totalfam,"_tau.txt")

      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/2.beta")
      write.table(data.frame('Score','var_Score_ver1','var_Score_ver2','Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)
      
      CEST_beta <- sapply(1:500,DoTest.beta,fin.dat=fin.dat,totalfam=totalfam,model=model,init_beta=init_beta,test.beta='std_snp',prev=prev,init_h2=h2,n.cores=n.cores,out=out)
    }
  }
}

source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
n.cores=24
for(type in c('power')){
  for(h2 in c(0.2, 0.4)){
    for(prev in c(0.1,0.2)){
      totalfam <- 500
      print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam,', Test type: ',type))
      library(kinship2)
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/2.beta")
      fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"_",type,"/dataset.txt"),head=T,stringsAsFactor=F)
      fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
      model <- Y~std_snp
      
      test.beta <- 'std_snp'
      init_h2 <- h2
      init_beta=matrix(1,1,1) # initial value for beta which is not for the H0.
      out <- paste0("prev_",prev,"_h2_",h2,"_",type,"/CEST_beta_",totalfam,"_tau.txt")

      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/2.beta")
      write.table(data.frame('Score','var_Score_ver1','var_Score_ver2','Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)
      
      CEST_beta <- sapply(1:500,DoTest.beta,fin.dat=fin.dat,totalfam=totalfam,model=model,init_beta=init_beta,test.beta='std_snp',prev=prev,init_h2=h2,n.cores=n.cores,out=out)
    }
  }
}


##### Summarize #####


setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
for(prev in c(0.1,0.2)){
  for(h2 in c(0,0.2,0.4)){
    for(totalfam in c(500)){
      IN <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,".txt")
      dat <- read.table(IN,head=T,stringsAsFactor=F)
      dat[,'Pvalue'] <- pchisq(dat[,'Chisq'],df=1,lower.tail=F)/2
      write.table(dat,IN,row.names=F,quote=F)
      
      IN <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,"_ver2.txt")
      dat <- read.table(IN,head=T,stringsAsFactor=F)
      dat[,'Pvalue'] <- pchisq(dat[,'Chisq'],df=1,lower.tail=F)/2
      write.table(dat,IN,row.names=F,quote=F)
      
    }
  }
}

setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
for(prev in c(0.1,0.2)){
  for(h2 in c(0,0.2,0.4)){
    for(totalfam in c(500)){
      IN <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,".txt")
      OUT <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_prev_",prev,"_h2_",h2,"_",totalfam)
      get.Summary(IN,OUT,prev,h2,totalfam,alpha=c(0.005,0.01,0.05,0.1),type=ifelse(h2==0,'size','power'),dist='mixture')
      system(paste0("scp ",OUT,"_power.txt wjkim@ng:~"))
      if(h2==0) system(paste0("scp ",OUT,"_QQ.png wjkim@ng:~"))
    }
  }
}

setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
for(prev in c(0.1,0.2)){
  for(h2 in c(0,0.2,0.4)){
    for(totalfam in c(500)){
      IN <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,"_ver2.txt")
      OUT <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_prev_",prev,"_h2_",h2,"_",totalfam,"_ver2")
      get.Summary(IN,OUT,prev,h2,totalfam,alpha=c(0.005,0.01,0.05,0.1),type=ifelse(h2==0,'size','power'),dist='mixture')
      system(paste0("scp ",OUT,"_power.txt wjkim@ng:~"))
      if(h2==0) system(paste0("scp ",OUT,"_QQ.png wjkim@ng:~"))
    }
  }
}

# Beta
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/2.beta")
for(prev in c(0.1,0.2)){
  for(h2 in c(0.2,0.4)){
    for(type in c('size')){
      for(totalfam in c(500)){
        IN <- paste0("prev_",prev,"_h2_",h2,"_",type,"/CEST_beta_",totalfam,"_tau.txt")
        OUT <- paste0("prev_",prev,"_h2_",h2,"_",type,"/CEST_beta_prev_",prev,"_h2_",h2,"_",totalfam,"_tau_",type)
        get.Summary(IN,OUT,prev,h2,totalfam,alpha=c(0.005,0.01,0.05,0.1),type=type,dist='chisq')
        system(paste0("scp ",OUT,"_power.txt wjkim@ng:~"))
        if(type=='size') system(paste0("scp ",OUT,"_QQ.png wjkim@ng:~"))
      }
      
    }
  }
}




##### get max logL
getlogL <- function(model,V,famid,prev,dataset,n.cores=1){
  ## This function depends on packages; tmvtnorm, parallele, mvtnorm
  ## model : Y~X 
  ## V : genetic relationship matrix (theoritically, kinship coefficient matrix)
  ## famid : family id
  ## prev : disease prevalence
  ## dataset : data.frame containing variables specified at formula
  
  if(!require(tmvtnorm))	{
    install.packages("tmvtnorm")
    library(tmvtnorm)
  }	
  
  if(!require(parallel))	{
    install.packages("parallel")
    library(parallel)
  }
  
  if(!require(mvtnorm))	{
    install.packages("mvtnorm")
    library(mvtnorm)
  }
  
  
  Y <- dataset[,as.character(model)[2]]
  
  logL <- function(h2){
    getab <- function(i){
      # i : family id
      fam <- which(famid==i)
      nn <- length(fam)			# number of family members for family i
      Yi <- matrix(Y[fam],nrow=nn)		# disease status for family i
      Vi <- V[fam,fam]			# GRM matric for family i
      Si <- h2*Vi+(1-h2)*diag(nn)		# Sigma for family i
      
      # calculating A and B
      k <- which(Yi==1) # aff
      kk <- which(Yi==0) # un
      t <- qnorm(prev,0,1,lower.tail=F)
      a <- rep(-Inf,nn); a[k] <- t
      b <- rep(Inf,nn); b[kk] <- t
      
      return(list(fam=fam,nn=nn,Vi=Vi,Si=Si,a=a,b=b))
    }
    resab <- lapply(unique(famid),getab)
    
    O <- Reduce('+',lapply(1:length(unique(famid)),function(i) with(resab[[i]],log(pmvnorm(lower=a,upper=b,mean=rep(0,nn),sig=Si)[[1]]))))
    return(c(h2,O))
  }
  n.point <- seq(-0.5,0.5,by=0.01)
  Y.point <- do.call(rbind,mclapply(n.point,logL,mc.cores=n.cores))
  colnames(Y.point) <- c("n.point","Y.point")
  return(Y.point)
}

for(prev in c(0.1,0.2)){
  h2=0
  totalfam=1000
  n.cores=24
  
  library(kinship2)
  setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
  fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
  
  model <- Y~snp
  dataset <- fin.dat[fin.dat$FID%in%sample(unique(fin.dat$FID),totalfam),,drop=F]
  famid <- dataset$FID
  total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
  V <- 2*as.matrix(kinship(total_ped))
  logL <- getlogL(model,V,famid,prev,dataset,n.cores)
  png(paste0("logL_prev_",prev,".png"))
  plot(logL[,1],logL[,2])
  dev.off()
}





###########################################################
############### New simul with seed #######################
###########################################################

###### Scenario 1, n14, CEST_h2_500_all_181116.txt
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
n.cores=32

for(h2 in c(0.2, 0.4)){
  for(prev in c(0.05,0.1,0.2)){
    for(totalfam in c(500)){
      print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
      library(kinship2)
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
      fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
      init_beta=matrix(0.2,1,1)
      model <- Y~std_snp-1
      out <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,"_all_181116.txt")
      
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
      if(paste0("CEST_h2_",totalfam,"_all_181116.txt") %in% system(paste0("ls prev_",prev,"_h2_",h2),intern=T)){
        start <- as.numeric(strsplit(system(paste0("wc -l ",out),intern=T)," ")[[1]][1])
      } else {
        start <- 1
        setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
        write.table(data.frame('Score','var_Score_ver1','var_Score_ver2','Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)
      }
      CEST_h2 <- sapply(start:2000,DoTest.h2,fin.dat=fin.dat,totalfam=totalfam,model=model,init_beta=init_beta,prev=prev,h2=h2,n.cores=n.cores,out=out,seed=T)
    }
  }
}


###### Scenario 2, CEST_h2_500_asc_181120.txt
## data generation
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/gen_simuldata.r')
for(prev in c(0.1,0.2)){
  for(h2 in c(0, 0.2, 0.4)){
    print(paste0('prevalence:',prev,', heritability:',h2))
    setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
    system(paste0("mkdir prev_",prev,"_h2_",h2))
    dataset = genNucFam(totalfam=50000,MAF=0.2,h2,ha2=0.05,prev,num_snp=1,n.cores=24)
    write.table(dataset,paste0("prev_",prev,"_h2_",h2,"/dataset_asc.txt"),row.names=F,quote=F)
  }
}



source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
do.CEST.S2 <- function(PREV,H2,N.fam,n.cores){
	for(h2 in H2){
	  for(prev in PREV){
		for(totalfam in N.fam){
		  print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
		  library(kinship2)
		  setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
		  fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset_asc.txt"),head=T,stringsAsFactor=F)
		  fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
		  fin.dat <- fin.dat[fin.dat$FID%in%fin.dat$FID[fin.dat$ind==1 & fin.dat$Y==1],,drop=F]
		  
		  init_beta=matrix(0.2,1,1)
		  model <- Y~std_snp-1
		  out <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,"_asc_181120.txt")
		  
		  setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
		  if(paste0("CEST_h2_",totalfam,"_asc_181120.txt") %in% system(paste0("ls prev_",prev,"_h2_",h2),intern=T)){
			start <- as.numeric(strsplit(system(paste0("wc -l ",out),intern=T)," ")[[1]][1])
		  } else {
			start <- 1
			setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
			write.table(data.frame('Score','var_Score_ver1','var_Score_ver2','Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)
		  }
		  CEST_h2 <- sapply(start:2000,DoTest.h2,fin.dat=fin.dat,totalfam=totalfam,model=model,init_beta=init_beta,prev=prev,h2=h2,n.cores=n.cores,out=out,seed=T,Proband="ind")
		}
	  }
	}
}

# n4, 24 cores
do.CEST.S2(0.05,0,500,24)

# n5, 24 cores
do.CEST.S2(0.05,0.2,500,24)

# n6, 24 cores
do.CEST.S2(0.05,0.4,500,24)

# n9, 24 cores
do.CEST.S2(0.1,0,500,24)

# n10, 24 cores
do.CEST.S2(0.1,0.2,500,24)

# n11, 24 cores
do.CEST.S2(0.1,0.4,500,24)

# n14, 32 cores
do.CEST.S2(0.2,0.2,500,32)

# n15, 32 cores
do.CEST.S2(0.2,0.4,500,32)

# n15, 32 cores
do.CEST.S2(0.2,0,500,32)

