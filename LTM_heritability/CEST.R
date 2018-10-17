#############################################################################################
############### Conditional Expectation Score Test (CEST) Oct 17, 2018 ######################
#############################################################################################

########## H0 : h2=0
Testing.h2 <- function(famid,V,init_beta,prev,X,Y,learning.rate=10,n.cores=1){

  ### X : N x p dim
  ### Y : N x 1 dim
  ### init_beta : P x 1 dim

  ## Theshold
  t <- qnorm(prev,0,1,lower.tail=F)
  
  beta.old <- init_beta
  n.iter = 0
  while(1){
    n.iter = n.iter+1
    ## NR algorithm to get estimate of beta under the h2=0
    phi <- dnorm(t-X%*%beta.old)
    PHI <- pnorm(X%*%beta.old)
    B0 <- X%*%beta.old+ifelse(Y==0,-phi/PHI,phi/(1-PHI))
  
    f <- t(X)%*%B0-t(X)%*%X%*%beta.old
    B0_prime <- do.call(rbind,lapply(1:length(Y),function(ij) return(ifelse(Y[ij,1]==0,
                                                                   X[ij,,drop=F]-phi[ij,1]*((t-X[ij,,drop=F]%*%beta.old)/PHI[ij,1]+phi[ij,1]/(PHI[ij,1])^2)*X[ij,,drop=F],
                                                                   X[ij,,drop=F]+phi[ij,1]*((t-X[ij,,drop=F]%*%beta.old)/(1-PHI[ij,1])-phi[ij,1]/(1-PHI[ij,1])^2)*X[ij,1]))))
    f_prime <- t(X)%*%B0_prime - t(X)%*%X
    beta.new <- beta.old - learning.rate*solve(f_prime)%*%f
    
    epsilon <- sqrt(sum((beta.old-beta.new)^2))
    print(data.frame(beta.new,epsilon,n.iter))
    if(epsilon<1e-4){
      break
    } else {
      beta.old <- beta.new
    }
  }

  phi <- dnorm(t-X%*%beta.new)
  PHI <- pnorm(X%*%beta.new)
  B0 <- X%*%beta.new+ifelse(Y==0,-phi/PHI,phi/(1-PHI))
  var0 <- 1+ifelse(Y==0,-(t-X%*%beta.new)*phi/PHI-(phi/PHI)^2,(t-X%*%beta.new)*phi/(1-PHI)-(phi/(1-PHI))^2)
  A0 <- var0+B0^2
    
  ## Score under the H0: h2=0
  FAMID <- unique(famid)
  tr <- function(m) sum(diag(m))
  
  get.Score <- function(fid){
    fam.idx <- which(famid==fid)
    Vi <- V[fam.idx,fam.idx,drop=F]
    Xi <- X[fam.idx,,drop=F]
    A0i <- A0[fam.idx,,drop=F]
    B0i <- B0[fam.idx,,drop=F]
    
    S0i_1 <- t(Xi)%*%B0i - t(Xi)%*%Xi%*%beta.new
    S0i_2 <- -tr(Vi-diag(length(fam.idx)))/2 + tr((Vi-diag(length(fam.idx)))%*%A0i)/2-t(Xi%*%beta.new)%*%(Vi-diag(length(fam.idx)))%*%(B0i-Xi%*%beta.new/2)
    S0i <- rbind(S0i_1,S0i_2)
    return(S0i)
  }
  
  Score <- Reduce('+',lapply(FAMID,get.Score))
  SS <- Score%*%t(Score)
  W <- solve(SS)
  
  
  
}


