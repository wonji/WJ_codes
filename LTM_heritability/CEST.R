#############################################################################################
############### Conditional Expectation Score Test (CEST) Oct 17, 2018 ######################
#############################################################################################

########## H0 : h2=0
Testing.h2 <- function(famid,V,init_beta,prev,X,Y,learning.rate=10,n.cores=1){
  library(tmvtnorm)
  library(parallel)
  
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
	  resAB <- lapply(1:length(Y),getAB)
	  O1 <- Reduce('+',lapply(1:length(Y),function(ij) resAB[[ij]]$O1ij))
	  O2 <- Reduce('+',lapply(1:length(Y),function(ij) resAB[[ij]]$O2ij))
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
  resAB <- lapply(1:length(Y),getAB)
  finalAB <- do.call(rbind,lapply(1:length(Y),function(ij) as.matrix(data.frame(Aij=resAB[[ij]]$Aij,Bij=resAB[[ij]]$Bij))))

  ## Score under the H0: h2=0
  FAMID <- unique(famid)
  tr <- function(m) sum(diag(m))
  
  get.Score <- function(fid){
    fam.idx <- which(famid==fid)
    Vi <- V[fam.idx,fam.idx,drop=F]
    Xi <- X[fam.idx,,drop=F]
    A0i <- finalAB[fam.idx,1,drop=F]
    B0i <- finalAB[fam.idx,2,drop=F]
    
    S0i_1 <- t(Xi)%*%B0i - t(Xi)%*%Xi%*%beta.new
    S0i_2 <- -tr(Vi-diag(length(fam.idx)))/2 + tr((Vi-diag(length(fam.idx)))%*%A0i)/2-t(Xi%*%beta.new)%*%(Vi-diag(length(fam.idx)))%*%(B0i-Xi%*%beta.new/2)
    S0i <- rbind(S0i_1,S0i_2)
    return(S0i)
  }
  
  Score <- Reduce('+',lapply(FAMID,get.Score))
  SS <- Score%*%t(Score)
  W <- solve(SS)
  
  
  
}


