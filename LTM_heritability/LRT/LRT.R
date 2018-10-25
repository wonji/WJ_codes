########################################################
####### Likelihood Ratio Test (LRT) for H0:h2=0 ########
########################################################

## Function to estimate beta under the H0:h2=0
LRT.getMLE.H0 <- function(famid,V,init_beta=NULL,prev,X=NULL,Y,n.cores=1){
  library(parallel)
  library(tmvtnorm)
  
  ### X : N x p dim
  ### Y : N x 1 dim
  ### init_beta : P x 1 dim
  
  ## Theshold
  t <- qnorm(prev,0,1,lower.tail=F)
  
  getAB <- function(ij){
    Yij <- Y[ij,1]
    if(!is.null(X)) Xij <- X[ij,,drop=F]
    a <- ifelse(Yij==0,-Inf,t)
    b <- ifelse(Yij==0,t,Inf)
    
    if(is.null(X)){
      para <- mtmvnorm(mean=rep(0,length(Yij)),sigma=1,lower=a,upper=b,doComputeVariance=TRUE)
    } else {
      para <- mtmvnorm(mean=as.vector(Xij%*%beta.old),sigma=1,lower=a,upper=b,doComputeVariance=TRUE)
    }
    Bij <- matrix(para$tmean,ncol=1)
    Aij <- para$tvar+Bij%*%t(Bij)
    
    if(is.null(X)){
      return(list(ij=ij,Yij=Yij,Aij=Aij,Bij=Bij))
    } else {
      O1ij <- t(Xij)%*%Xij
      O2ij <- t(Xij)%*%Bij
      return(list(ij=ij,Yij=Yij,Xij=Xij,Aij=Aij,Bij=Bij,O1ij=O1ij,O2ij=O2ij))
    }
  }
  
  if(!is.null(X)){
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
  }
  
  return(beta.new)
}

## Function to estimate h2 & beta under the (H0 U H1) == Just use estimation function!
LTMH <- function(model,init_beta,init_h2,V,famid,prev,data,max.iter=100,max.sub.iter=50,n.cores=1){
  ## This function depends on a package; tmvtnorm
  ## model : Y~X 
  ## init_beta : initial value for beta (p+1 vector)
  ## init_h2 :  initial value for heritability
  ## V : genetic relationship matrix (theoritically, kinship coefficient matrix)
  ## famid : family id
  ## prev : disease prevalence
  ## data : data.frame containing variables specified at formula
  
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
  
  # calculate A and B
  getAB <- function(i,h2.old,beta.old,famid){
    ## i : family id
    fam <- which(famid==i)
    nn <- length(fam)			# number of family members for family i
    Yi <- matrix(Y[fam],nrow=nn)		# disease status for family i
    Xi <- matrix(X[fam,],nrow=nn)		# design matric for family i
    Vi <- V[fam,fam]			# GRM matric for family i
    Si <- h2.old*Vi+(1-h2.old)*diag(nn)	# Sigma for family i
    
    # calculating A and B
    k <- which(Yi==1) # aff
    kk <- which(Yi==0) # unaff
    t <- qnorm(prev,0,1,lower.tail=F)
    a <- rep(-Inf,nn); a[k] <- t
    b <- rep(Inf,nn); b[kk] <- t
    para <- mtmvnorm(mean=as.vector(Xi%*%beta.old),sigma=Si,lower=a,upper=b,doComputeVariance=TRUE)
    Bi <- matrix(para$tmean,ncol=1)
    Ai <- para$tvar+Bi%*%t(Bi)
    
    return(list(fam=fam,a=a,b=b,Yi=Yi,Xi=Xi,Vi=Vi,nn=nn,Ai=Ai,Bi=Bi))
  }

  getAB_ij <- function(ij){
    Yij <- Y[ij]
    if(!is.null(X)) Xij <- X[ij,,drop=F]
    t <- qnorm(prev,0,1,lower.tail=F)
	a <- ifelse(Yij==0,-Inf,t)
    b <- ifelse(Yij==0,t,Inf)
    
    if(is.null(X)){
      para <- mtmvnorm(mean=rep(0,length(Yij)),sigma=1,lower=a,upper=b,doComputeVariance=TRUE)
    } else {
      para <- mtmvnorm(mean=as.vector(Xij%*%beta.old),sigma=1,lower=a,upper=b,doComputeVariance=TRUE)
    }
    Bij <- matrix(para$tmean,ncol=1)
    Aij <- para$tvar+Bij%*%t(Bij)
    
    if(is.null(X)){
      return(list(ij=ij,Yij=Yij,Aij=Aij,Bij=Bij))
    } else {
      O1ij <- t(Xij)%*%Xij
      O2ij <- t(Xij)%*%Bij
      return(list(ij=ij,Yij=Yij,Xij=Xij,Aij=Aij,Bij=Bij,O1ij=O1ij,O2ij=O2ij))
    }
  }
    
  getELE <- function(i,resAB,h2.old){
    # resAB : result of getAB
    res <- resAB[[i]]
    Si <- with(res,h2.old*Vi+(1-h2.old)*diag(nn))	# Sigma for family i
    invSi <- solve(Si)
    Ci <- with(res,-invSi%*%(Vi-diag(nn))%*%invSi)
    Hi <- with(res,-2*invSi%*%(Vi-diag(nn))%*%Ci)
    O1 <- with(res,t(Xi)%*%invSi%*%Xi)
    O2 <- with(res,t(Xi)%*%Ci%*%Xi)
    O3 <- with(res,t(Xi)%*%invSi%*%Bi)
    O4 <- with(res,t(Xi)%*%Ci%*%Bi)
    return(list(Si=Si,invSi=invSi,Ci=Ci,Hi=Hi,O1=O1,O2=O2,O3=O3,O4=O4))
  }
  
  getELE.0 <- function(i,resAB){
    # resAB : result of getAB
    res <- resAB[[i]]
    Si <- invSi <- with(res,diag(nn))	# Sigma & invSigma for family i
    Ci <- with(res,-(Vi-diag(nn)))
    Hi <- with(res,-2*(Vi-diag(nn))%*%Ci)
    O1 <- with(res,t(Xi)%*%Xi)
    O2 <- with(res,t(Xi)%*%Ci%*%Xi)
    O3 <- with(res,t(Xi)%*%Bi)
    O4 <- with(res,t(Xi)%*%Ci%*%Bi)
    return(list(Si=Si,invSi=invSi,Ci=Ci,Hi=Hi,O1=O1,O2=O2,O3=O3,O4=O4))
  }
  
  getELE.subiter <- function(i,resAB,h2.old){
    # resAB : result of getAB
    res <- resAB[[i]]
    Si <- with(res,h2.old*Vi+(1-h2.old)*diag(nn))	# Sigma for family i
    invSi <- solve(Si)
    O1 <- with(res,t(Xi)%*%invSi%*%Xi)
    O3 <- with(res,t(Xi)%*%invSi%*%Bi)
    return(list(Si=Si,invSi=invSi,O1=O1,O3=O3))
  }
  
  tr <- function(m) sum(diag(m))
  getff <- function(i){
    f <- with(resAB[[i]],with(resELE[[i]],-tr(invSi%*%(Vi-diag(nn)))/2-tr(Ci%*%Ai)/2+t(Xi%*%beta.hat)%*%Ci%*%(Bi-Xi%*%beta.hat/2))) 
    f_prime <- with(resAB[[i]],with(resELE[[i]],-tr(Ci%*%(Vi-diag(nn)))/2-tr(Hi%*%Ai)/2+t(Xi%*%F)%*%Ci%*%(Bi-Xi%*%beta.hat/2)+t(Xi%*%beta.hat)%*%Hi%*%(Bi-Xi%*%beta.hat/2)-t(Xi%*%beta.hat)%*%Ci%*%Xi%*%F/2))
    return(c(f,f_prime))
  }
  
  getf.01 <- function(i,h2){
    if(h2==0) f <- with(resAB[[i]],with(resELE[[i]],-tr(Vi-diag(nn))/2-tr(Ci%*%Ai)/2+t(Xi%*%beta.hat)%*%Ci%*%(Bi-Xi%*%beta.hat/2))) 
    if(h2==1) f <- with(resAB[[i]],with(resELE[[i]],-tr(invSi%*%(Vi-diag(nn)))/2-tr(Ci%*%Ai)/2+t(Xi%*%beta.hat)%*%Ci%*%(Bi-Xi%*%beta.hat/2))) 
    return(f)
  }
  
  Y <- data[,as.character(model)[2]]
  X <- model.matrix(model,data)
  beta.old <- init_beta
  h2.old <- init_h2
  
  n.iter = 0
  epsilon = 1
  while(1){
    if(n.iter==max.iter) {
      return(paste("not converge after ",max.iter," iterations",sep=""))
    }
    n.iter=n.iter+1
    
    ## Calculate A and B
    resAB <- mclapply(unique(famid),getAB,beta.old=beta.old,h2.old=h2.old,famid=famid,mc.cores=n.cores)
    
    #### Check constraint
    ## h2=0
    # Calculate elements
    resELE <- lapply(1:length(unique(famid)),getELE.0,resAB=resAB)
    
    O1 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O1))
    O2 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O2))
    O3 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O3))
    O4 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O4))
    
    D <- solve(O1)
    E <- -D%*%O2%*%D
    F <- E%*%O3+D%*%O4
    beta.hat <- D%*%O3
    
    # f and f_prime
    f <- sum(sapply(1:length(unique(famid)),getf.01,h2=0),na.rm=T)
    if(f<=0) {
      h2.new <- 0
      epsilon <- sqrt(sum((beta.old-beta.hat)^2))
      print(data.frame(h2.new,epsilon,n.iter))
	  
      if(epsilon<1e-5){
        beta.new <- beta.hat
        return(list(beta=beta.new,h2=h2.new,n_iter=n.iter,resAB=resAB))
      } else {
        h2.old <- h2.new
        beta.old <- beta.hat
		
		while(1){
		  resAB <- mclapply(1:length(Y),getAB_ij,mc.cores=n.cores)
		  O1 <- Reduce('+',lapply(1:length(Y),function(ij) resAB[[ij]]$O1ij))
		  O2 <- Reduce('+',lapply(1:length(Y),function(ij) resAB[[ij]]$O2ij))
		  beta.new <- solve(O1)%*%O2
		  epsilon <- sqrt(sum((beta.old-beta.new)^2))
		  print(data.frame(h2.new,epsilon,n.iter))
		  
		  beta.old <- beta.new
		  if(epsilon<1e-5) break
		}
      }			
    } else {
      ## h2=1
      # Calculate elements
      resELE <- lapply(1:length(unique(famid)),getELE,resAB=resAB,h2.old=1)
      
      O1 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O1))
      O2 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O2))
      O3 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O3))
      O4 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O4))
      
      D <- solve(O1)
      E <- -D%*%O2%*%D
      F <- E%*%O3+D%*%O4
      beta.hat <- D%*%O3
      
      f <- sum(sapply(1:length(unique(famid)),getf.01,h2=1),na.rm=T)
      if(f>=0){
        h2.new <- 1
        epsilon <- sqrt(sum((beta.old-beta.hat)^2))
        
        if(epsilon<1e-5){
          beta.new <- beta.hat
          print(data.frame(h2.new,epsilon,n.iter))
          return(list(beta=beta.new,h2=h2.new,n_iter=n.iter,resAB=resAB))
        } else {
          h2.old <- h2.new
          beta.old <- beta.hat
		  print(data.frame(h2.new,epsilon,n.iter))
        }			
      } else {
        
        # Sub-iterations
        n.sub.iter=0
        e=1
        while(1){
          if(n.sub.iter==max.sub.iter) {
            return(paste("not converge at ",n.iter," iterations",sep=""))
          }
          n.sub.iter=n.sub.iter+1
          
          # Calculate elements
          resELE <- lapply(1:length(unique(famid)),getELE,resAB=resAB,h2.old=h2.old)
          
          O1 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O1))
          O2 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O2))
          O3 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O3))
          O4 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O4))
          
          D <- solve(O1)
          E <- -D%*%O2%*%D
          F <- E%*%O3+D%*%O4
          beta.hat <- D%*%O3
          
          # f and f_prime
          ff <- rowSums(sapply(1:length(unique(famid)),getff))
          f <- ff[1]; f_prime <- ff[2]
          
          # update h2
          h2.new <- h2.old -f/f_prime
          e <- abs(h2.new-h2.old)
          
          if(e<1e-5){
            h2.old <- h2.new
            resELE <- lapply(1:length(unique(famid)),getELE.subiter,resAB=resAB,h2.old=h2.old)
            O1 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O1))
            D <- solve(O1)
            O3 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O3))
            beta.new <- D%*%O3
            print(data.frame(h2.new,epsilon,n.iter))
            break
          } else {
            h2.old <- h2.new
          }
        }
        
        epsilon <- sqrt(sum((beta.old-beta.new)^2))
        if(epsilon<1e-5){
          return(list(beta=beta.new,h2=h2.new,n_iter=n.iter,resAB=resAB))
        } else {
          beta.old <- beta.new
        }
      }
    }
  }
}			


LRT.DoTest.h2 <- function(obs,fin.dat,totalfam,model,init_beta=NULL,init_h2,prev,h2,n.cores=1,out){
  print(obs) 
  
  # Sampling test data
  dataset <- fin.dat[fin.dat$FID%in%sample(unique(fin.dat$FID),totalfam),,drop=F]
  Y <- dataset[,as.character(model)[2],drop=F]
  if(!is.null(init_beta)) {
    X <- model.matrix(model,dataset)
  } else {
    X <- NULL
  }
  famid <- dataset$FID
  total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
  V <- 2*as.matrix(kinship(total_ped))
  
  ### Get MLE
  # Under H0
  MLE.H0 <- LRT.getMLE.H0(famid=famid,V=V,init_beta=init_beta,prev=prev,X=X,Y=Y,n.cores=n.cores)
  
  # Under (H0 U H1)
  MLE.H01 <- LTMH(model=model,init_beta=init_beta,init_h2=init_h2,V=V,famid=famid,prev=prev,data=dataset,max.iter=100,max.sub.iter=50,n.cores=n.cores)
  
  ### Testing
  thres <- qnorm(prev,0,1,lower.tail=F)
  a <- ifelse(Y==0,-Inf,thres)
  b <- ifelse(Y==0,thres,Inf)
  
  # Under H0
  Lij <- log(pnorm(b,mean=X%*%MLE.H0,sd=1,lower.tail=T) - pnorm(a,mean=X%*%MLE.H0,sd=1,lower.tail=T))
  logL_H0 <- sum(Lij)
  
  # Under H0 U H1
  MLE.H01.beta <- MLE.H01$beta
  MLE.H01.h2 <- MLE.H01$h2
  
  get.faminfo <- function(ii){
	resAB <- MLE.H01$resAB[[ii]]
	Si <- MLE.H01.h2*resAB$Vi+(1- MLE.H01.h2)*diag(resAB$nn)
	return(list(Xi=resAB$Xi,Si=Si,a=resAB$a,b=resAB$b))
  }
  resFAM <- lapply(1:length(unique(famid)),get.faminfo)
  logL_H01 <- Reduce('+',lapply(1:length(unique(famid)),function(i) with(resFAM[[i]],sum(log((pmvnorm(lower=a,upper=b,mean=as.vector(Xi%*%MLE.H01.beta),sig=Si)[[1]]))))))
  
  Chisq <- -2*(logL_H0 - logL_H01)
  Pvalue <- ifelse(MLE.H01$h2==0,1/2,pchisq(Chisq,df=1,lower.tail=F)/2)
  Testing.res <- data.frame(MLE_beta=MLE.H01.beta,MLE_h2=MLE.H01.h2,Chisq,Pvalue)
  write.table(Testing.res,out,col.names=F,row.names=F,quote=F,append=TRUE)
}
