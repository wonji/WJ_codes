###################################################################
########### log-likelihood plot when score is negative #############
###################################################################

Testing.h2 <- function(famid,V,init_beta,prev,X,Y,n.cores=1){
  library(parallel)
  library(tmvtnorm)
  
  ### X : N x p dim
  ### Y : N x 1 dim
  ### init_beta : P x 1 dim
  
  ## Theshold
  t <- qnorm(prev,0,1,lower.tail=F)
  
  getAB <- function(ij){
    Yij <- Y[ij,1]
    Xij <- X[ij,,drop=F]
    a <- ifelse(Yij==0,-Inf,t)
    b <- ifelse(Yij==0,t,Inf)
    
    para <- mtmvnorm(mean=as.vector(Xij%*%beta.old),sigma=1,lower=a,upper=b,doComputeVariance=TRUE)
    Bij <- matrix(para$tmean,ncol=1)
    Aij <- para$tvar+Bij%*%t(Bij)
    
    O1ij <- t(Xij)%*%Xij
    O2ij <- t(Xij)%*%Bij
    
    return(list(ij=ij,Yij=Yij,Xij=Xij,Aij=Aij,Bij=Bij,O1ij=O1ij,O2ij=O2ij))
  }
  
  beta.old <- init_beta
  n.iter = 0
  while(1){
    n.iter = n.iter+1
    resAB <- mclapply(1:nrow(Y),getAB,mc.cores=n.cores)
    O1 <- Reduce('+',lapply(1:nrow(Y),function(ij) resAB[[ij]]$O1ij))
    O2 <- Reduce('+',lapply(1:nrow(Y),function(ij) resAB[[ij]]$O2ij))
    beta.new <- solve(O1)%*%O2
    
    epsilon <- sqrt(sum((beta.old-beta.new)^2))
    print(data.frame(beta.new,epsilon,n.iter))
    if(epsilon<1e-5){
      beta.old <- beta.new
      break
    } else {
      beta.old <- beta.new
    }
  }
  
  #### Now, we have MLE for beta under the H0: h2=0.
  resAB <- mclapply(1:nrow(Y),getAB,mc.cores=n.cores)
  finalAB <- do.call(rbind,lapply(1:nrow(Y),function(ij) as.matrix(data.frame(obs=resAB[[ij]]$ij,Aij=resAB[[ij]]$Aij,Bij=resAB[[ij]]$Bij))))
  finalAB <- finalAB[order(finalAB[,'obs']),]
  
  ## Score under the H0: h2=0
  FAMID <- unique(famid)
  tr <- function(m) sum(diag(m))
  
  get.Score <- function(fid){
    fam.idx <- which(famid==fid)
    Vi <- V[fam.idx,fam.idx,drop=F]
    Xi <- X[fam.idx,,drop=F]
    B0i <- finalAB[fam.idx,'Bij',drop=F]
    A0i <- B0i%*%t(B0i)
    diag(A0i) <- finalAB[fam.idx,'Aij']
    
    S0i_1 <- t(Xi)%*%B0i - t(Xi)%*%Xi%*%beta.new
    S0i_2 <- -tr(Vi-diag(length(fam.idx)))/2 + tr((Vi-diag(length(fam.idx)))%*%A0i)/2-t(Xi%*%beta.new)%*%(Vi-diag(length(fam.idx)))%*%(B0i-Xi%*%beta.new/2)
    S0i <- rbind(S0i_1,S0i_2)
    return(S0i)
  }
  
  Scores <- lapply(FAMID,get.Score)
  S <- Reduce('+',Scores)
  h2.idx <- nrow(S)
  S.h2 <- S[h2.idx,,drop=F]
  
  # version 1 variance
  M <- Reduce('+',lapply(1:length(FAMID),function(i) Scores[[i]]%*%t(Scores[[i]]))) - S%*%t(S)/length(FAMID)
  varS.h2.1 <- M[h2.idx,h2.idx,drop=F]-M[h2.idx,-h2.idx,drop=F]%*%solve(M[-h2.idx,-h2.idx,drop=F])%*%M[-h2.idx,h2.idx,drop=F]
  
  # version 1 variance
  M <- Reduce('+',lapply(1:length(FAMID),function(i) Scores[[i]]%*%t(Scores[[i]])))
  varS.h2.2 <- M[h2.idx,h2.idx,drop=F]-M[h2.idx,-h2.idx,drop=F]%*%solve(M[-h2.idx,-h2.idx,drop=F])%*%M[-h2.idx,h2.idx,drop=F]
  
  if(S.h2[1,1]<=0){
    T.h2 <- 0
    pval <- 1/2
  } else{
    T.h2 <- S.h2%*%solve(varS.h2.1)%*%S.h2
    pval <- pchisq(T.h2,df=1,lower.tail=F)/2
  }
  fin.res <- data.frame(Score=S.h2,Beta=beta.new,var_Score_ver1=varS.h2.1,var_Score_ver2=varS.h2.2,Chisq=T.h2,Pvalue=pval)
  return(fin.res)
}


getlogL <- function(model,V,famid,prev,dataset,n.cores=1,Beta){
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
  X <- model.matrix(model,data=dataset)
  
  logL <- function(h2){
    getab <- function(i){
      # i : family id
      fam <- which(famid==i)
      nn <- length(fam)			# number of family members for family i
      Yi <- matrix(Y[fam],nrow=nn)		# disease status for family i
      Xi <- X[fam,,drop=F]
      Vi <- V[fam,fam]			# GRM matric for family i
      Si <- h2*Vi+(1-h2)*diag(nn)		# Sigma for family i
      
      # calculating A and B
      k <- which(Yi==1) # aff
      kk <- which(Yi==0) # un
      t <- qnorm(prev,0,1,lower.tail=F)
      a <- rep(-Inf,nn); a[k] <- t
      b <- rep(Inf,nn); b[kk] <- t
      
      return(list(fam=fam,nn=nn,Vi=Vi,Si=Si,a=a,b=b,Xi=Xi))
    }
    resab <- lapply(unique(famid),getab)
    
    O <- Reduce('+',lapply(1:length(unique(famid)),function(i) with(resab[[i]],log(pmvnorm(lower=a,upper=b,mean=as.vector(Xi%*%Beta),sig=Si)[[1]]))))
    return(c(h2,O))
  }
  n.point <- seq(-0.5,1.5,by=0.01)
  Y.point <- do.call(rbind,mclapply(n.point,logL,mc.cores=n.cores))
  colnames(Y.point) <- c("n.point","Y.point")
  return(Y.point)
}



library(kinship2)
prev=0.1
h2=0
n.cores=24
totalfam=500
setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
model <- Y~std_snp-1

ii=0
while(1){
  
  # Sampling test data
  dataset <- fin.dat[fin.dat$FID%in%sample(unique(fin.dat$FID),totalfam),,drop=F]
  Y <- dataset[,as.character(model)[2],drop=F]
  X <- model.matrix(model,dataset)
  init_beta <- matrix(0.2,1,1)
  famid <- dataset$FID
  total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
  V <- 2*as.matrix(kinship(total_ped))
  
  # Testing
  res <- Testing.h2(famid,V,init_beta,prev,X,Y,n.cores)
  
  if(res[1,1]<=0){
    ii=ii+1
    print(paste0(ii,'th log-Likelihood'))
    # Loglikelihood
    logLik <- getlogL(model,V,famid,prev,dataset,n.cores,Beta=matrix(res[,'Beta']))
    png(paste0("logL_prev_",prev,"_negaScore_",ii,".png"))
    plot(logLik[,1],logLik[,2])
    dev.off()
    if(ii==10) break
  }

}



####### Distribution of Score

prev=0.2
h2=0
totalfam=500
alpha=c(0.01,0.05,0.1)

setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
out <- paste0("prev_",prev,"_h2_",h2,"/CEST_h2_",totalfam,"_all.txt")
dat <- read.table(out,head=T)
Stat.norm <- dat[,'Score']/sqrt(dat[,'var_Score_ver1'])
png(paste0("prev_",prev,"_h2_",h2,"/Score_dist_prev",prev,".png"))
h <- hist(Stat.norm,breaks=seq(-4.25,3.25,0.5))
xfit <- seq(min(Stat.norm),max(Stat.norm),length=40)
yfit <- dnorm(xfit)
yfit <- yfit*diff(h$mids[1:2])*length(Stat.norm) 
lines(xfit,yfit,col='blue',lwd=2)
dev.off()


get.power <- function(alpha){
   Tstat <- dat[,'Score']/sqrt(dat[,'var_Score_ver2'])
   pval <- pnorm(Tstat,lower.tail=F)
   res <- lapply(alpha,function(i) data.frame(pvalue=sum(pval<i)/length(pval)))
   res.I <- do.call(rbind,res)
   rownames(res.I) <- alpha
   colnames(res.I) <- "Proportion"
   return(res.I)  
}

