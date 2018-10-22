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

###### Ver 2 : Remove Information matrix latter part
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
      dataset = genNucFam(totalfam=10000,MAF=0.2,h2,ha2=0.05,prev,num_snp=1)
      write.table(dataset,paste0("prev_",prev,"_h2_",h2,"_",ifelse(num_snp==0,"size","power"),"/dataset.txt"),row.names=F,quote=F)
    }
  }
}



source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/CEST.R')
n.cores=24
for(prev in c(0.1,0.2)){
  for(h2 in c(0.2, 0.4)){
    for(type in c('size','power')){
      totalfam <- 500
      print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam,' Test type: ',type))
      library(kinship2)
      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/2.beta")
      fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"_",type,"/dataset.txt"),head=T,stringsAsFactor=F)
      fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
      model <- Y~std_snp-1
      
      test.beta <- 'std_snp'
      init_h2 <- h2
      init_beta=matrix(1,1,1) # initial value for beta which is not for the H0.
      out <- paste0("prev_",prev,"_h2_",h2,"_",type,"/CEST_beta_",totalfam,"_tau.txt")

      setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/2.beta")
      write.table(data.frame('Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)
      
      CEST_h2 <- sapply(1:2000,DoTest.beta,fin.dat=fin.dat,totalfam=totalfam,model=model,init_beta=init_beta,prev=prev,h2=h2,n.cores=n.cores,out=out)
    }
  }
}


##### Summarize #####

get.Summary <- function(IN,OUT,prev,h2,totalfam,alpha){
  ## Table
  dat <- read.table(IN,head=T)
  pval <- dat[,2]
  new.pval <- pchisq(dat[,1],df=1,lower.tail=F)
  res <- lapply(alpha,function(i) data.frame(pvalue=sum(pval<i)/length(pval),new_pvalue=sum(new.pval<i)/length(new.pval)))
  res.I <- do.call(rbind,res)
  rownames(res.I) <- alpha
  colnames(res.I) <- c("Version 1","Version 2")
  write.table(res.I,paste0(OUT,"_power.txt"),quote=F)

  ## QQ plots
  p.val <- cbind(pval,new.pval)
  png(paste0(OUT,"_QQ.png"), width=2000, height=1000, units = "px", bg = "white", res = 200)
  par(mfrow=c(1,2))

  title <- c("Version 1","Version 2")
  for(i in 1:2){
    y <- -log(p.val[,i],10)

    # QQ plot
    o.y <- sort(y[is.finite(y)],decreasing=T)
    xx<- (1:length(o.y))/length(o.y)
    x <- -log(xx,10)
    ifelse(max(o.y[is.finite(o.y)])>=x[1],YY<-o.y,YY<-x)
    
    plot(YY, YY, main=title[i-1],type = "n", xlab = expression(paste("Expected ",-log[10],"(p-value)")), ylab = expression(paste("Observed ",-log[10],"(p-value)")))
    N=length(o.y)
    c95 <- rep(0,N)
    c05 <- rep(0,N)
    for(i in 1:N){
      c95[i] <- qbeta(0.95,i,N-i+1)
      c05[i] <- qbeta(0.05,i,N-i+1)
    }
    abline(h = c(0:max(max(x),max(o.y[is.finite(o.y)]))), v =c(0:max(max(x),max(o.y[is.finite(o.y)]))), col = "darkgray", lty=3)
    polygon(c(x,sort(x)),c(-log(c95,10),sort(-log(c05,10))),col=c("gray"),border=NA)
    abline(a=0,b=1,col="black", lty=1)
    points(x, o.y, cex = .5, col = "dark red")
    
  }
  dev.off()
  return(res.I)
}

setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
for(prev in c(0.1,0.2)){
  for(h2 in c(0, 0.2, 0.4)){
    for(totalfam in c(500)){
      IN <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,".txt")
      OUT <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_prev_",prev,"_h2_",h2,"_",totalfam)
      get.Summary(IN,OUT,prev,h2,totalfam,alpha=c(0.005,0.01,0.05,0.1))
      system(paste0("scp ",OUT,"_* wjkim@ng:~"))
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
  n.cores=32
  
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

