#######################################################
####################### LTMH ##########################
#######################################################
########### MLE for h2, beta ###########

LTMH <- function(model,data,init_beta=NULL,init_h2,V,famid,prev,max.iter=100,max.sub.iter=50,n.cores=1){

  if(!require(tmvtnorm))	{
    install.packages("tmvtnorm")
    library(tmvtnorm)
  }	
  
  if(!require(parallel))	{
    install.packages("parallel")
    library(parallel)
  }
  
  # calculate A and B
  getAB <- function(i,h2.old,beta.old,famid){
    ## i : family id
    fam <- which(famid==i)
    nn <- length(fam)			# number of family members for family i
    Yi <- matrix(Y[fam],nrow=nn)		# disease status for family i
    if(!is.null(X)) {
      Xi <- matrix(X[fam,],nrow=nn)		# design matric for family i
    } else {
      Xi <- NULL
    }
    Vi <- V[fam,fam]			# GRM matric for family i
    Si <- h2.old*Vi+(1-h2.old)*diag(nn)	# Sigma for family i
    rownames(Si) <- colnames(Si)
    
    # calculating A and B
    k <- which(Yi==1) # aff
    kk <- which(Yi==0) # unaff
    t <- qnorm(prev,0,1,lower.tail=F)
    a <- rep(-Inf,nn); a[k] <- t
    b <- rep(Inf,nn); b[kk] <- t
    if(!is.null(X)){
      para <- mtmvnorm(mean=as.vector(Xi%*%beta.old),sigma=Si,lower=a,upper=b,doComputeVariance=TRUE)
    } else {
      para <- mtmvnorm(mean=rep(0,nn),sigma=Si,lower=a,upper=b,doComputeVariance=TRUE)
    }
    
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
      para <- mtmvnorm(mean=0,sigma=1,lower=a,upper=b,doComputeVariance=TRUE)
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
    
    if(!is.null(res$Xi)){
      O1 <- with(res,t(Xi)%*%invSi%*%Xi)
      O2 <- with(res,t(Xi)%*%Ci%*%Xi)
      O3 <- with(res,t(Xi)%*%invSi%*%Bi)
      O4 <- with(res,t(Xi)%*%Ci%*%Bi)
      return(list(Si=Si,invSi=invSi,Ci=Ci,Hi=Hi,O1=O1,O2=O2,O3=O3,O4=O4))
    } else {
      return(list(Si=Si,invSi=invSi,Ci=Ci,Hi=Hi))
    }
  }
  
  getELE.0 <- function(i,resAB){
    # resAB : result of getAB
    res <- resAB[[i]]
    Si <- invSi <- with(res,diag(nn))	# Sigma & invSigma for family i
    Ci <- with(res,-(Vi-diag(nn)))
    Hi <- with(res,-2*(Vi-diag(nn))%*%Ci)
    
    if(!is.null(res$Xi)){
      O1 <- with(res,t(Xi)%*%Xi)
      O2 <- with(res,t(Xi)%*%Ci%*%Xi)
      O3 <- with(res,t(Xi)%*%Bi)
      O4 <- with(res,t(Xi)%*%Ci%*%Bi)
      return(list(Si=Si,invSi=invSi,Ci=Ci,Hi=Hi,O1=O1,O2=O2,O3=O3,O4=O4))
    } else {
      return(list(Si=Si,invSi=invSi,Ci=Ci,Hi=Hi))
    }
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
    if(!is.null(X)){
      f <- with(resAB[[i]],with(resELE[[i]],-tr(invSi%*%(Vi-diag(nn)))/2-tr(Ci%*%Ai)/2+t(Xi%*%beta.hat)%*%Ci%*%(Bi-Xi%*%beta.hat/2))) 
      f_prime <- with(resAB[[i]],with(resELE[[i]],-tr(Ci%*%(Vi-diag(nn)))/2-tr(Hi%*%Ai)/2+t(Xi%*%F)%*%Ci%*%(Bi-Xi%*%beta.hat/2)+t(Xi%*%beta.hat)%*%Hi%*%(Bi-Xi%*%beta.hat/2)-t(Xi%*%beta.hat)%*%Ci%*%Xi%*%F/2))
    } else {
      f <- with(resAB[[i]],with(resELE[[i]],-tr(invSi%*%(Vi-diag(nn)))/2-tr(Ci%*%Ai)/2)) 
      f_prime <- with(resAB[[i]],with(resELE[[i]],-tr(Ci%*%(Vi-diag(nn)))/2-tr(Hi%*%Ai)/2))
    }
    return(c(f,f_prime))
  }
  
  getf.01 <- function(i,h2){
    if(!is.null(X)){
      if(h2==0) f <- with(resAB[[i]],with(resELE[[i]],-tr(Vi-diag(nn))/2-tr(Ci%*%Ai)/2+t(Xi%*%beta.hat)%*%Ci%*%(Bi-Xi%*%beta.hat/2))) 
      if(h2==1) f <- with(resAB[[i]],with(resELE[[i]],-tr(invSi%*%(Vi-diag(nn)))/2-tr(Ci%*%Ai)/2+t(Xi%*%beta.hat)%*%Ci%*%(Bi-Xi%*%beta.hat/2))) 
    } else {
      if(h2==0) f <- with(resAB[[i]],with(resELE[[i]],-tr(Vi-diag(nn))/2-tr(Ci%*%Ai)/2)) 
      if(h2==1) f <- with(resAB[[i]],with(resELE[[i]],-tr(invSi%*%(Vi-diag(nn)))/2-tr(Ci%*%Ai)/2))
    }
    
    return(f)
  }
  
  Y <- data[,as.character(model)[2]]
  h2.old <- h2.OLD <- init_h2
  
  if(as.character(model)[3]=='-1'){
    X <- NULL
  } else {
    X <- model.matrix(model,data)
    notint <- which(colnames(X)!='(Intercept)')
    if(length(notint)!=0) {
      XX <- X
      for(jj in notint) X[,jj] <- (XX[,jj]-mean(XX[,jj],na.rm=T))/sd(XX[,jj],na.rm=T)
    }
    if(as.character(model)[3]!='-1' & is.null(init_beta))
      init_beta <- glm(Y~X-1,family=binomial(link=probit))$coef
    beta.old <- matrix(init_beta,ncol=1)
  }
  
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
    
    if(!is.null(X)){
      O1 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O1))
      O2 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O2))
      O3 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O3))
      O4 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O4))
      
      D <- solve(O1)
      E <- -D%*%O2%*%D
      F <- E%*%O3+D%*%O4
      beta.hat <- D%*%O3
    }
    
    # f and f_prime
    f <- sum(sapply(1:length(unique(famid)),getf.01,h2=0),na.rm=T)
    if(f<=0) {
      h2.new <- 0
      if(is.null(X)){
        epsilon <- abs(h2.old-h2.new)
        print(data.frame(h2.new,epsilon,n.iter))
        if(epsilon<1e-5){
          logL <- Reduce('+',lapply(1:length(unique(famid)),function(i) with(resAB[[i]],sum(log((pmvnorm(lower=a,upper=b,mean=rep(0,nn),sig=h2.new*Vi+(1-h2.new)*diag(nn))[[1]]))))))
          return(list(h2=h2.new,n_iter=n.iter))
        } else {
          h2.old <- h2.new
        }
      } else {
        epsilon <- sqrt(sum((beta.old-beta.hat)^2))
        print(data.frame(h2.new,epsilon,n.iter))
        
        if(epsilon<1e-5){
          beta.new <- beta.hat
          beta.std <- data.frame(t(beta.new)); colnames(beta.std) <- colnames(X)
          if(ncol(X)!=1 & '(Intercept)'%in%colnames(X)){
            Mean <- colMeans(XX)
            SD <- apply(XX,2,sd)
            for(gg in 2:ncol(XX)){
              beta.new[gg,1] <- beta.new[gg,1]/SD[gg]
              beta.new[1,1] <- beta.new[1,1]-Mean[gg]/SD[gg]
            }
            beta.unstd <- data.frame(t(beta.new)); colnames(beta.unstd) <- colnames(X)
          } else if(!'(Intercept)'%in%colnames(X)){
            Mean <- colMeans(XX)
            SD <- apply(XX,2,sd)
            itc <- 0
            for(gg in 1:ncol(XX)){
              beta.new[gg,1] <- beta.new[gg,1]/SD[gg]
              itc <- itc - Mean[gg]/SD[gg]
            }
            beta.unstd <- data.frame(t(beta.new)); colnames(beta.unstd) <- colnames(X)
            beta.unstd <- data.frame(itc,beta.unstd)
            colnames(beta.unstd)[1] <- '(Intercept)'
          } else {
            beta.unstd <- beta.std
          }
          logL <- Reduce('+',lapply(1:length(unique(famid)),function(i) with(resAB[[i]],sum(log((pmvnorm(lower=a,upper=b,mean=as.vector(Xi%*%matrix(as.matrix(beta.std))),sig=h2.new*Vi+(1-h2.new)*diag(nn))[[1]]))))))
          return(list(beta_std=beta.std,beta_unstd=beta.unstd,h2=h2.new,n_iter=n.iter,logL=logL))
        } else {
          h2.old <- h2.new
          beta.old <- beta.hat
          
          while(1){
            n.iter=n.iter+1
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
      }
    } else {
      ## h2=1
      # Calculate elements
      resELE <- lapply(1:length(unique(famid)),getELE,resAB=resAB,h2.old=1)
      
      if(!is.null(X)){
        O1 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O1))
        O2 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O2))
        O3 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O3))
        O4 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O4))
        
        D <- solve(O1)
        E <- -D%*%O2%*%D
        F <- E%*%O3+D%*%O4
        beta.hat <- D%*%O3
      }
      
      f <- sum(sapply(1:length(unique(famid)),getf.01,h2=1),na.rm=T)
      if(f>=0){
        h2.new <- 1
        if(is.null(X)){
          epsilon <- abs(h2.old-h2.new)
          print(data.frame(h2.new,epsilon,n.iter))
          if(epsilon<1e-5){
            logL <- Reduce('+',lapply(1:length(unique(famid)),function(i) with(resAB[[i]],sum(log((pmvnorm(lower=a,upper=b,mean=rep(0,nn),sig=h2.new*Vi+(1-h2.new)*diag(nn))[[1]]))))))
            return(list(h2=h2.new,n_iter=n.iter,logL=logL))
          } else {
            h2.old <- h2.new
          }
        } else {
          epsilon <- sqrt(sum((beta.old-beta.hat)^2))
          print(data.frame(h2.new,epsilon,n.iter))
          
          if(epsilon<1e-5){
            beta.new <- beta.hat
            beta.std <- data.frame(t(beta.new)); colnames(beta.std) <- colnames(X)
            if(ncol(X)!=1 & '(Intercept)'%in%colnames(X)){
              Mean <- colMeans(XX)
              SD <- apply(XX,2,sd)
              for(gg in 2:ncol(XX)){
                beta.new[gg,1] <- beta.new[gg,1]/SD[gg]
                beta.new[1,1] <- beta.new[1,1]-Mean[gg]/SD[gg]
              }
              beta.unstd <- data.frame(t(beta.new)); colnames(beta.unstd) <- colnames(X)
            } else if(!'(Intercept)'%in%colnames(X)){
              Mean <- colMeans(XX)
              SD <- apply(XX,2,sd)
              itc <- 0
              for(gg in 1:ncol(XX)){
                beta.new[gg,1] <- beta.new[gg,1]/SD[gg]
                itc <- itc - Mean[gg]/SD[gg]
              }
              beta.unstd <- data.frame(t(beta.new)); colnames(beta.unstd) <- colnames(X)
              beta.unstd <- data.frame(itc,beta.unstd)
              colnames(beta.unstd)[1] <- '(Intercept)'
            } else {
              beta.unstd <- beta.std
            }
            logL <- Reduce('+',lapply(1:length(unique(famid)),function(i) with(resAB[[i]],sum(log((pmvnorm(lower=a,upper=b,mean=as.vector(Xi%*%matrix(as.matrix(beta.std))),sig=h2.new*Vi+(1-h2.new)*diag(nn))[[1]]))))))
            return(list(beta_std=beta.std,beta_unstd=beta.unstd,h2=h2.new,n_iter=n.iter,logL=logL)) 
          } else {
            h2.old <- h2.new
            beta.old <- beta.hat
          }			
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
          
          if(!is.null(X)){
            O1 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O1))
            O2 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O2))
            O3 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O3))
            O4 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O4))
            
            D <- solve(O1)
            E <- -D%*%O2%*%D
            F <- E%*%O3+D%*%O4
            beta.hat <- D%*%O3
          }        
          # f and f_prime
          ff <- rowSums(sapply(1:length(unique(famid)),getff))
          f <- ff[1]; f_prime <- ff[2]
          
          # update h2
          h2.new <- h2.old -f/f_prime
          e <- abs(h2.new-h2.old)
          
          if(e<1e-5){
            if(!is.null(X)){
              h2.old <- h2.new
              resELE <- lapply(1:length(unique(famid)),getELE.subiter,resAB=resAB,h2.old=h2.old)
              O1 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O1))
              D <- solve(O1)
              O3 <- Reduce('+',lapply(1:length(unique(famid)),function(i) resELE[[i]]$O3))
              beta.new <- D%*%O3
            }
            break
          } else {
            h2.old <- h2.new
          }
        }
        
        if(!is.null(X)){
          epsilon <- sqrt(sum((beta.old-beta.new)^2))
          print(data.frame(h2.new,epsilon,n.iter))
          if(epsilon<1e-5){
            beta.std <- data.frame(t(beta.new)); colnames(beta.std) <- colnames(X)
            if(ncol(X)!=1 & '(Intercept)'%in%colnames(X)){
              Mean <- colMeans(XX)
              SD <- apply(XX,2,sd)
              for(gg in 2:ncol(XX)){
                beta.new[gg,1] <- beta.new[gg,1]/SD[gg]
                beta.new[1,1] <- beta.new[1,1]-Mean[gg]/SD[gg]
              }
              beta.unstd <- data.frame(t(beta.new)); colnames(beta.unstd) <- colnames(X)
            } else if(!'(Intercept)'%in%colnames(X)){
              Mean <- colMeans(XX)
              SD <- apply(XX,2,sd)
              itc <- 0
              for(gg in 1:ncol(XX)){
                beta.new[gg,1] <- beta.new[gg,1]/SD[gg]
                itc <- itc - Mean[gg]/SD[gg]
              }
              beta.unstd <- data.frame(t(beta.new)); colnames(beta.unstd) <- colnames(X)
              beta.unstd <- data.frame(itc,beta.unstd)
              colnames(beta.unstd)[1] <- '(Intercept)'
            } else {
              beta.unstd <- beta.std
            }
            logL <- Reduce('+',lapply(1:length(unique(famid)),function(i) with(resAB[[i]],sum(log((pmvnorm(lower=a,upper=b,mean=as.vector(Xi%*%matrix(as.matrix(beta.std))),sig=h2.new*Vi+(1-h2.new)*diag(nn))[[1]]))))))
            return(list(beta_std=beta.std,beta_unstd=beta.unstd,h2=h2.new,n_iter=n.iter,logL=logL))
          } else {
            beta.old <- beta.new
          }
        } else {
          epsilon <- abs(h2.OLD-h2.new)
          print(data.frame(h2.new,epsilon,n.iter))
          if(epsilon<1e-5){
            logL <- Reduce('+',lapply(1:length(unique(famid)),function(i) with(resAB[[i]],sum(log((pmvnorm(lower=a,upper=b,mean=rep(0,nn),sig=h2.new*Vi+(1-h2.new)*diag(nn))[[1]]))))))
            return(list(h2=h2.new,n_iter=n.iter,logL=logL))
          } else {
            h2.old <- h2.OLD <- h2.new
          }
        }
      }
    }
  }
}			

########################## Ascertainment adjust LTMH

LTMH.asc <- function(model,data,init_beta=NULL,init_h2,V,famid,prev,max.iter=100,n.cores=1,proband,SAVE=F,out=NULL){

  if(as.character(model)[3]=='-1'){
    res <- LTMH(model=new.model,init_beta=new.init_beta,init_h2=init_h2,V=V,famid=famid,prev=prev,data=data,n.cores=n.cores)
    return(res)
  }
  
  if(!require(tmvtnorm))	{
    install.packages("tmvtnorm")
    library(tmvtnorm)
  }	
  
  if(!require(parallel))	{
    install.packages("parallel")
    library(parallel)
  }
  
  tr <- function(m) sum(diag(m))
  
  # calculate A and B
	# calculate A and B
	getAB <- function(i,h2.old,beta.old,famid,proband){
		## i : family id
		fam <- which(famid==i)
		PB <- proband[fam]
		nn <- length(fam)			# number of family members for family i
		Yi <- matrix(Y[fam],nrow=nn)		# disease status for family i
		Xi <- matrix(X[fam,],nrow=nn)		# design matric for family i
		Vi <- V[fam,fam]			# GRM matric for family i
		Si <- h2.old*Vi+(1-h2.old)*diag(nn)	# Sigma for family i
		rownames(Si) <- colnames(Si)

		# calculating A and B
		k <- which(Yi==1) # aff
		kk <- which(Yi==0) # unaff
		t <- qnorm(prev,0,1,lower.tail=F)
		a <- rep(-Inf,nn); a[k] <- t
		b <- rep(Inf,nn); b[kk] <- t
		if(h2.old!=0){
		  para <- mtmvnorm(mean=as.vector(Xi%*%beta.old),sigma=Si,lower=a,upper=b,doComputeVariance=TRUE)
		  Bi <- matrix(para$tmean,ncol=1)
		  Ai <- para$tvar+Bi%*%t(Bi)
		} else {
		  getAB_ij <- function(jj) {
		    para <- mtmvnorm(mean=as.vector(Xi[jj,,drop=F]%*%beta.old),sigma=1,lower=a[jj],upper=b[jj],doComputeVariance = TRUE)
		    Bi <- matrix(para$tmean,ncol=1)
		    Ai <- para$tvar+Bi%*%t(Bi)
		    return(c(Ai,Bi))
		  }
		  res <- lapply(1:nn,getAB_ij)
		  res <- do.call(rbind,res)
		  Bi <- res[,2,drop=F]
		  Ai <- Bi%*%t(Bi)
		  diag(Ai) <- res[,1]
		}

		return(list(fam=fam,a=a,b=b,Yi=Yi,Xi=Xi,Vi=Vi,nn=nn,Ai=Ai,Bi=Bi,PB=PB,thres=t))
	}
	
	getELE <- function(i,resAB,h2.old,beta.old){
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
		
		dBeta <- O3 - O1%*%beta.old
		dh2 <- with(res,-tr(invSi%*%(Vi-diag(nn)))/2 - tr(Ci%*%Ai)/2 + t(Xi%*%beta.old)%*%Ci%*%(Bi-Xi%*%beta.old/2))
		
		d2Beta <- -O1
		d2h2Beta <- O4 - O2%*%beta.old
		d2h2 <- with(res,-tr(Ci%*%(Vi-diag(nn)))/2 - tr(Hi%*%Ai)/2 + t(Xi%*%beta.old)%*%Hi%*%(Bi-Xi%*%beta.old/2))
		
		return(list(Si=Si,invSi=invSi,Ci=Ci,Hi=Hi,O1=O1,O2=O2,O3=O3,O4=O4,dBeta=dBeta,dh2=dh2,d2Beta=d2Beta,d2h2Beta=d2h2Beta,d2h2=d2h2))
	}

	getELE.0 <- function(i,resAB,beta.old){
		# resAB : result of getAB
		res <- resAB[[i]]
		Si <- invSi <- with(res,diag(nn))	# Sigma for family i
		Ci <- with(res,-(Vi-diag(nn)))
		Hi <- with(res,-2*(Vi-diag(nn))%*%Ci)
		O1 <- with(res,t(Xi)%*%Xi)
		O2 <- with(res,t(Xi)%*%Ci%*%Xi)
		O3 <- with(res,t(Xi)%*%Bi)
		O4 <- with(res,t(Xi)%*%Ci%*%Bi)
		
		dBeta <- O3 - O1%*%beta.old
		dh2 <- with(res,-tr(Vi-diag(nn))/2 - tr(Ci%*%Ai)/2 + t(Xi%*%beta.old)%*%Ci%*%(Bi-Xi%*%beta.old/2))
		
		d2Beta <- -O1
		d2h2Beta <- O4 - O2%*%beta.old
		d2h2 <- with(res,-tr(Ci%*%(Vi-diag(nn)))/2 - tr(Hi%*%Ai)/2 + t(Xi%*%beta.old)%*%Hi%*%(Bi-Xi%*%beta.old/2))
		
		return(list(Si=Si,invSi=invSi,Ci=Ci,Hi=Hi,O1=O1,O2=O2,O3=O3,O4=O4,dBeta=dBeta,dh2=dh2,d2Beta=d2Beta,d2h2Beta=d2h2Beta,d2h2=d2h2))
	}
	
	getELE.PB <- function(i,resAB,beta.old){
		# resAB : result of getAB
		res <- resAB[[i]]

		Yip <- with(res,Yi[PB==1])
		Xip <- with(res,Xi[PB==1,drop=F])
		Mui <- with(res,1-pnorm(thres-Xip%*%beta.old,lower.tail=T))
		Alphai <- log(Mui/(1-Mui))
		dBeta <- with(res,(Yip-Mui)/(Mui*(1-Mui))*dnorm(thres-Xip%*%beta.old)*t(Xip))
		d2Beta <- with(res,dnorm(thres-Xip%*%beta.old)/(Mui*(1-Mui))*(-dnorm(thres-Xip%*%beta.old)+((Yip-Mui)*(2*Mui-1)*dnorm(thres-Xip%*%beta.old))/(Mui*(1-Mui))+(Yip-Mui)*(thres-Xip%*%beta.old))*Xip%*%t(Xip))
		
		return(list(dBeta=dBeta,d2Beta=d2Beta))
	}
  

  Y <- data[,as.character(model)[2]]
  if(as.character(model)[3]=='-1'){
    X <- NULL
  } else {
    X <- model.matrix(model,data)
    notint <- which(colnames(X)!='(Intercept)')
    if(length(notint)!=0) {
      XX <- X
      for(jj in notint) X[,jj] <- (XX[,jj]-mean(XX[,jj],na.rm=T))/sd(XX[,jj],na.rm=T)
    }
    if(as.character(model)[3]!='-1' & is.null(init_beta))
      init_beta <- glm(Y~X-1,family=binomial(link=probit))$coef
    beta.old <- beta.OLD <- matrix(init_beta,ncol=1)
  }
  h2.old <- h2.OLD <- init_h2
  theta.old <- theta.OLD <- rbind(beta.old,h2.old)
  
  ## Estimate h2
  n.iter = 0
  epsilon = 1
  beta.OLD <- beta.old
  while(1){
    if(n.iter==max.iter) {
      return(paste("not converge after ",max.iter," iterations",sep=""))
    }
    print(data.frame(h2=theta.OLD[nrow(theta.OLD),],epsilon=epsilon,n.iter=n.iter))
    n.iter=n.iter+1
    
    ## Calculate A and B
    resAB <- mclapply(unique(famid),getAB,beta.old=beta.old,h2.old=h2.old,famid=famid,proband=proband,mc.cores=n.cores)
    
    ## Calculate first and second derivative
    while(1){
      
      resELE <- mclapply(1:length(unique(famid)),getELE,resAB=resAB,beta.old=beta.old,h2.old=h2.old,mc.cores=n.cores)
      resELE.PB <- mclapply(1:length(unique(famid)),getELE.PB,resAB=resAB,beta.old=beta.old,mc.cores=n.cores)
      
      ## f & f_prime
      f_Q <- Reduce('+',mclapply(1:length(unique(famid)),function(i) with(resELE[[i]],c(dBeta,dh2)),mc.cores=n.cores))
      J_Q <- Reduce('+',mclapply(1:length(unique(famid)),function(i) with(resELE[[i]],rbind(cbind(d2Beta,d2h2Beta),cbind(t(d2h2Beta),d2h2))),mc.cores=n.cores))
      
      # proband
      f.PB <- Reduce('+',mclapply(1:length(unique(famid)),function(i) with(resELE.PB[[i]],dBeta),mc.cores=n.cores))
      J.PB <- Reduce('+',mclapply(1:length(unique(famid)),function(i) with(resELE.PB[[i]],d2Beta),mc.cores=n.cores))
      
      f <- f_Q
      J <- J_Q
      f[-length(f)] <- f_Q[-length(f_Q)]-f.PB
      J[-nrow(J),-nrow(J)] <- J_Q[-nrow(J_Q),-nrow(J_Q)] - J.PB
      
      theta.new <- theta.old - solve(J)%*%f
      
      e <- sqrt(sum((theta.old-theta.new)^2))
      theta.old <- theta.new
      beta.old <- theta.old[-nrow(theta.old),,drop=F]
      h2.old <- theta.old[nrow(theta.old),]
      
      if(e<1e-5) break
    }
    
    if(h2.old<=0){
      h2.old <- 0
      while(1){
        resAB <- mclapply(unique(famid),getAB,beta.old=beta.old,h2.old=0,famid=famid,proband=proband,mc.cores=n.cores)
        
        while(1){
          resELE <- mclapply(1:length(unique(famid)),getELE.0,resAB=resAB,beta.old=beta.OLD,mc.cores=n.cores)
          resELE.PB <- mclapply(1:length(unique(famid)),getELE.PB,resAB=resAB,beta.old=beta.OLD,mc.cores=n.cores)
          
          ## f & f_prime
          f_Q <- Reduce('+',mclapply(1:length(unique(famid)),function(i) with(resELE[[i]],dBeta),mc.cores=n.cores))
          J_Q <- Reduce('+',mclapply(1:length(unique(famid)),function(i) with(resELE[[i]],d2Beta),mc.cores=n.cores))
          
          # proband
          f.PB <- Reduce('+',mclapply(1:length(unique(famid)),function(i) with(resELE.PB[[i]],dBeta),mc.cores=n.cores))
          J.PB <- Reduce('+',mclapply(1:length(unique(famid)),function(i) with(resELE.PB[[i]],d2Beta),mc.cores=n.cores))
          
          f <- f_Q - f.PB
          J <- J_Q - J.PB
          
          beta.NEW <- beta.OLD - solve(J)%*%f
          epsilon <- sqrt(sum((beta.NEW-beta.OLD)^2))
          beta.OLD <- beta.NEW
          if(epsilon < 1e-5) break
        }
        
        epsilon <- sqrt(sum((beta.OLD - beta.old)^2))
        beta.old <- beta.OLD
        if(epsilon < 1e-5) break
      }
      theta.new <- rbind(beta.old,h2.old)
    }
    
    if(h2.old>=1){
      h2.old <- 1
      while(1){
        resAB <- mclapply(unique(famid),getAB,beta.old=beta.old,h2.old=1,famid=famid,proband=proband,mc.cores=n.cores)
        
        while(1){
          resELE <- mclapply(1:length(unique(famid)),getELE,resAB=resAB,beta.old=beta.OLD,h2.old=1,mc.cores=n.cores)
          resELE.PB <- mclapply(1:length(unique(famid)),getELE.PB,resAB=resAB,beta.old=beta.OLD,mc.cores=n.cores)
          
          ## f & f_prime
          f_Q <- Reduce('+',mclapply(1:length(unique(famid)),function(i) with(resELE[[i]],dBeta),mc.cores=n.cores))
          J_Q <- Reduce('+',mclapply(1:length(unique(famid)),function(i) with(resELE[[i]],d2Beta),mc.cores=n.cores))
          
          # proband
          f.PB <- Reduce('+',mclapply(1:length(unique(famid)),function(i) with(resELE.PB[[i]],dBeta),mc.cores=n.cores))
          J.PB <- Reduce('+',mclapply(1:length(unique(famid)),function(i) with(resELE.PB[[i]],d2Beta),mc.cores=n.cores))
          
          f <- f_Q - f.PB
          J <- J_Q - J.PB
          
          beta.NEW <- beta.OLD - solve(J)%*%f
          epsilon <- sqrt(sum((beta.NEW-beta.OLD)^2))
          beta.OLD <- beta.NEW
          if(epsilon < 1e-5) break
        }
        
        epsilon <- sqrt(sum((beta.OLD - beta.old)^2))
        beta.old <- beta.OLD
        if(epsilon < 1e-5) break
      }
      theta.new <- rbind(beta.old,h2.old)
    }
    epsilon <- sqrt(sum((theta.OLD-theta.new)^2))
    theta.OLD <- theta.new
    if(SAVE){
      output <- list(beta=beta.old,h2=h2.old,n_iter=n.iter)
      write.table(t(as.matrix(c(unlist(output),n_iter=n.iter))),out,row.names=F,col.names=F,append=T)
    }
    if(epsilon<1e-5) {
	  beta.new <- beta.std <- theta.new[-nrow(theta.new),,drop=F]
	  colnames(beta.std) <- colnames(X)
	  h2.new <- theta.new[nrow(theta.new),]
	  names(h2.new) <- NULL
      if(ncol(X)!=1 & '(Intercept)'%in%colnames(X)){
        Mean <- colMeans(XX)
        SD <- apply(XX,2,sd)
        for(gg in 2:ncol(XX)){
          beta.new[gg,1] <- beta.new[gg,1]/SD[gg]
          beta.new[1,1] <- beta.new[1,1]-Mean[gg]/SD[gg]
        }
        beta.unstd <- data.frame(t(beta.new))
		colnames(beta.unstd) <- colnames(X)
      } else if(!'(Intercept)'%in%colnames(X)){
        Mean <- colMeans(XX)
        SD <- apply(XX,2,sd)
        itc <- 0
        for(gg in 1:ncol(XX)){
          beta.new[gg,1] <- beta.new[gg,1]/SD[gg]
          itc <- itc - Mean[gg]/SD[gg]
        }
        beta.unstd <- data.frame(t(beta.new))
		colnames(beta.unstd) <- colnames(X)
        beta.unstd <- data.frame(itc,beta.unstd)
        colnames(beta.unstd)[1] <- '(Intercept)'
      } else {
        beta.unstd <- beta.std
      }
	  output <- LTMH(model=model,data=dataset,init_beta=NULL,init_h2=0.2,V=V,famid=famid,prev=prev,n.cores=20)
    return(list(beta_std=beta.std,beta_unstd=beta.unstd,h2=h2.new,n_iter=n.iter,logL=logL))
    break
    }
  }
}





#######################################################
####################### CEST ##########################
#######################################################


########## H0 : h2=0
CEST.h2 <- function(model,data,famid,V,init_beta=NULL,prev,n.cores=1,proband=NULL){
  library(parallel)
  library(tmvtnorm)
  
  ## Theshold
  thres <- qnorm(prev,0,1,lower.tail=F)
  
  getAB <- function(ij,Y,XX,thres,beta.old){
    Yij <- Y[ij]
    if(!is.null(X)) Xij <- XX[ij,,drop=F]
    a <- ifelse(Yij==0,-Inf,thres)
    b <- ifelse(Yij==0,thres,Inf)
    
    if(is.null(X)){
      para <- mtmvnorm(mean=0,sigma=1,lower=a,upper=b,doComputeVariance=TRUE)
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
  
  getELE <- function(ij,Y,XX,beta.old,thres){
    Yij <- Y[ij]
    Xij <- XX[ij,,drop=F]
    Muij <- 1-pnorm(thres-Xij%*%beta.old,lower.tail=T)
    dBeta <- (Yij-Muij)/(Muij*(1-Muij))*dnorm(thres-Xij%*%beta.old)*t(Xij)
    d2Beta <- dnorm(thres-Xij%*%beta.old)/(Muij*(1-Muij))*(-dnorm(thres-Xij%*%beta.old)+((Yij-Muij)*(2*Muij-1)*dnorm(thres-Xij%*%beta.old))/(Muij*(1-Muij))+(Yij-Muij)*(thres-Xij%*%beta.old))*Xij%*%t(Xij)
    
    return(list(dBeta=dBeta,d2Beta=d2Beta))
  }

  Y <- data[,as.character(model)[2]]
  if(as.character(model)[3]=='-1'){
    X <- NULL
  } else {
    X <- model.matrix(model,data)
    notint <- which(colnames(X)!='(Intercept)')
    if(length(notint)!=0) {
      XX <- X
      for(jj in notint) X[,jj] <- (XX[,jj]-mean(XX[,jj],na.rm=T))/sd(XX[,jj],na.rm=T)
    }
    if(as.character(model)[3]!='-1' & is.null(init_beta))
      init_beta <- glm(Y~X-1,family=binomial(link=probit))$coef
    beta.old <- matrix(init_beta,ncol=1)
  }

  if(!is.null(X)){
    if(is.null(proband)){
      Y.actual <- Y
      X.actual <- X
    } else {
      Y.actual <- Y[proband==0]
      X.actual <- X[proband==0,,drop=F]
    }
    beta.old <- matrix(init_beta,ncol=1)
    n.iter = 0
    epsilon <- 1
    while(1){
      print(data.frame(beta.old,epsilon,n.iter))
      n.iter <- n.iter+1
      resELE <- mclapply(1:length(Y.actual),getELE,beta.old=beta.old,Y=Y.actual,XX=X.actual,thres=thres,mc.cores=n.cores)
      
      # proband
      f <- Reduce('+',mclapply(1:length(Y.actual),function(ij) with(resELE[[ij]],dBeta),mc.cores=n.cores))
      J <- Reduce('+',mclapply(1:length(Y.actual),function(ij) with(resELE[[ij]],d2Beta),mc.cores=n.cores))
      
      beta.new <- beta.old - solve(J)%*%f
      epsilon <- sqrt(sum((beta.new-beta.old)^2))
      beta.old <- beta.new
      if(epsilon < 1e-5) break
    }
  }
  
  #### Now, we have MLE for beta under the H0: h2=0.
  resAB <- mclapply(1:length(Y),getAB,Y=Y,XX=X,thres=thres,beta.old=beta.new,mc.cores=n.cores)
  finalAB <- do.call(rbind,mclapply(1:length(Y),function(ij) as.matrix(data.frame(obs=resAB[[ij]]$ij,Aij=resAB[[ij]]$Aij,Bij=resAB[[ij]]$Bij)),mc.cores=n.cores))
  finalAB <- finalAB[order(finalAB[,'obs']),]
  
  ## Score under the H0: h2=0
  FAMID <- unique(famid)
  tr <- function(m) sum(diag(m))
  
  get.Score <- function(fid){
    fam.idx <- which(famid==fid)
    Vi <- V[fam.idx,fam.idx,drop=F]
    if(!is.null(X)) Xi <- X[fam.idx,,drop=F]
    B0i <- finalAB[fam.idx,'Bij',drop=F]
    A0i <- B0i%*%t(B0i)
    diag(A0i) <- finalAB[fam.idx,'Aij']
    
    if(is.null(X)){
      S0i <- matrix(-tr(Vi-diag(length(fam.idx)))/2 + tr((Vi-diag(length(fam.idx)))%*%A0i)/2,1,1)
    } else {
      S0i_1 <- t(Xi)%*%B0i - t(Xi)%*%Xi%*%beta.new
      S0i_2 <- -tr(Vi-diag(length(fam.idx)))/2 + tr((Vi-diag(length(fam.idx)))%*%A0i)/2-t(Xi%*%beta.new)%*%(Vi-diag(length(fam.idx)))%*%(B0i-Xi%*%beta.new/2)
      S0i <- rbind(S0i_1,S0i_2)
    }
    return(S0i)
  }
  
  Scores <- lapply(FAMID,get.Score)
  S <- Reduce('+',Scores)
  h2.idx <- nrow(S)
  if(!is.null(proband)){
    Y.PB <- Y[proband==1]
    X.PB <- X[proband==1,,drop=F]
    resELE <- mclapply(1:length(Y.PB),getELE,beta.old=beta.new,Y=Y.PB,XX=X.PB,thres=thres,mc.cores=n.cores)
    S.PB <- Reduce('+',mclapply(1:length(Y.PB),function(ij) with(resELE[[ij]],dBeta),mc.cores=n.cores))
    S[-h2.idx,] <- S[-h2.idx,,drop=F] - S.PB
  }
  S.h2 <- S[h2.idx,,drop=F]
  
  # version 1 variance
  M <- Reduce('+',lapply(1:length(FAMID),function(i) Scores[[i]]%*%t(Scores[[i]]))) - S%*%t(S)/length(FAMID)
  varS.h2.1 <- ifelse(is.null(X),M,M[h2.idx,h2.idx,drop=F]-M[h2.idx,-h2.idx,drop=F]%*%solve(M[-h2.idx,-h2.idx,drop=F])%*%M[-h2.idx,h2.idx,drop=F])
  
  # version 2 variance
  M <- Reduce('+',lapply(1:length(FAMID),function(i) Scores[[i]]%*%t(Scores[[i]])))
  varS.h2.2 <- ifelse(is.null(X),M,M[h2.idx,h2.idx,drop=F]-M[h2.idx,-h2.idx,drop=F]%*%solve(M[-h2.idx,-h2.idx,drop=F])%*%M[-h2.idx,h2.idx,drop=F])
  
  if(S.h2[1,1]<=0){
    T.h2 <- 0
    pval <- 1/2
  } else{
    T.h2 <- t(S.h2)%*%solve(varS.h2.1)%*%S.h2
    pval <- pchisq(T.h2,df=1,lower.tail=F)/2
  }
  fin.res <- data.frame(Score=S.h2,var_Score_ver1=varS.h2.1,Chisq=T.h2,DF=1,Pvalue=pval)
  colnames(fin.res) <- c('Score','var_Score','Chisq','DF','Pvalue')
  return(fin.res)
}



########## H0 : beta=0
CEST.beta <- function(model,data,init_h2,V,famid,prev,test.beta,max.iter=100,n.cores=1,proband=NULL){

  if(!require(tmvtnorm))	{
    install.packages("tmvtnorm")
    library(tmvtnorm)
  }	
  
  if(!require(parallel))	{
    install.packages("parallel")
    library(parallel)
  }
  
  # calculate A and B
  getAB <- function(i,h2.old,beta.old,famid){
    ## i : family id
    fam <- which(famid==i)
    nn <- length(fam)			# number of family members for family i
    Yi <- Y[fam,,drop=F]		# disease status for family i
    Xi <- X[fam,,drop=F]		# design matric for family i
    Vi <- V[fam,fam,drop=F]			# GRM matric for family i
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
    Yij <- Y[ij,1]
    if(!is.null(X)) Xij <- X[ij,,drop=F]
    t <- qnorm(prev,0,1,lower.tail=F)
    a <- ifelse(Yij==0,-Inf,t)
    b <- ifelse(Yij==0,t,Inf)
    
    if(is.null(X)){
      para <- mtmvnorm(mean=0,sigma=1,lower=a,upper=b,doComputeVariance=TRUE)
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
  
  model.comp <- as.character(model)
  model.comp[3] <- gsub(test.beta,'',model.comp[3])
  if(model.comp[3]=='') model.comp[3] <- 1
  new.model <- as.formula(paste0(model.comp[2],'~',model.comp[3]))
  
  if(model.comp[3]=='-1'){
    new.init_beta <- NULL
  } else {
    new.init_beta <- glm(new.model,data=data,family=binomial(link=probit))$coef
  }
    
  if(is.null(proband)){
    res <- LTMH(model=new.model,init_beta=new.init_beta,init_h2=init_h2,V=V,famid=famid,prev=prev,data=data,n.cores=n.cores)
  } else {
    res <- LTMH.asc(model=new.model,init_beta=new.init_beta,init_h2=init_h2,V=V,famid=famid,prev=prev,data=data,n.cores=n.cores,proband=proband)
  }

  #### Now, we have MLE for h2 under the H0: beta=0.
  FAMID <- unique(famid)
  X <- model.matrix(model,data=data)
  test.beta.idx <- which(colnames(X)%in%test.beta)
  final.beta <- matrix(0,ncol(X),1)
  if(!is.null(res$beta)) final.beta[-test.beta.idx,1] <- res$beta[,1]
  
  resAB <- mclapply(FAMID,getAB,h2.old=h2.new,beta.old=final.beta,famid=famid,mc.cores=n.cores)
  
  ## Score under the H0: beta=0
  get.Score <- function(ii){
    
    fam.idx <- resAB[[ii]]$fam
    Vi <- resAB[[ii]]$Vi
    Si <- h2.new*Vi+(1-h2.new)*diag(length(fam.idx))
    invSi <- solve(Si)
    
    Ci <- -invSi%*%(Vi-diag(length(fam.idx)))%*%invSi
    Xi <- matrix(X[fam.idx,,drop=F],nrow=length(fam.idx))
    B0i <- resAB[[ii]]$Bi
    A0i <- resAB[[ii]]$Ai
    
    S0i_1 <- t(Xi)%*%invSi%*%B0i - t(Xi)%*%invSi%*%Xi%*%final.beta
    S0i_2 <- -tr(invSi%*%(Vi-diag(length(fam.idx))))/2 - tr(Ci%*%A0i)/2 + t(Xi%*%final.beta)%*%Ci%*%(B0i-Xi%*%final.beta/2)
    S0i <- rbind(S0i_1,S0i_2)
    return(S0i)
  }
  
  Scores <- lapply(1:length(FAMID),get.Score)
  S <- Reduce('+',Scores)
  
  # version 1 variance
  M <- Reduce('+',lapply(1:length(FAMID),function(i) Scores[[i]]%*%t(Scores[[i]]))) - S%*%t(S)/length(FAMID)
  varS.beta.1 <- NA
  try(varS.beta.1 <- M[test.beta.idx,test.beta.idx,drop=F]-M[test.beta.idx,-test.beta.idx,drop=F]%*%solve(M[-test.beta.idx,-test.beta.idx,drop=F])%*%M[-test.beta.idx,test.beta.idx,drop=F],silent=TRUE)
  
  # version 2 variance
  M <- Reduce('+',lapply(1:length(FAMID),function(i) Scores[[i]]%*%t(Scores[[i]])))
  varS.beta.2 <- NA
  try(varS.beta.2 <- M[test.beta.idx,test.beta.idx,drop=F]-M[test.beta.idx,-test.beta.idx,drop=F]%*%solve(M[-test.beta.idx,-test.beta.idx,drop=F])%*%M[-test.beta.idx,test.beta.idx,drop=F],silent=TRUE)
  
  S.beta <- S[test.beta.idx,,drop=F]
  if(!is.na(varS.beta.1)){
    T.beta <- t(S.beta)%*%solve(varS.beta.1)%*%S.beta
    pval <- pchisq(T.beta,df=length(test.beta.idx),lower.tail=F)
  } else if(!is.na(varS.beta.2)){
    T.beta <- t(S.beta)%*%solve(varS.beta.2)%*%S.beta
    pval <- pchisq(T.beta,df=length(test.beta.idx),lower.tail=F)
  } else {
    T.beta <- pval <- NA
  }
  
  return(data.frame(Score=S.beta,var_Score=ifelse(!is.na(varS.beta.1),varS.beta.1,varS.beta.2),Chisq=T.beta,DF=length(test.beta),length(test.beta.idx),Pvalue=pval))
}			

########## H0 : beta=0
LRT.beta <- function(model,data,init_h2,V,famid,prev,test.beta,max.iter=100,n.cores=1,proband=NULL){

  model.comp <- as.character(model)
  model.comp[3] <- gsub(test.beta,'',model.comp[3])
  model.comp[3] <- gsub('\\+\\s+\\+','\\+',model.comp[3])
  if(model.comp[3]=='') model.comp[3] <- 1
  
  new.model <- as.formula(paste0(model.comp[2],'~',model.comp[3]))
  
  # Under H0
  if(model.comp[3]=='-1'){
    new.init_beta <- NULL
  } else {
    new.init_beta <- glm(new.model,data=data,family=binomial(link=probit))$coef
  }
  if(is.null(proband)){
    MLE.H0 <- LTMH(model=new.model,init_beta=new.init_beta,init_h2=init_h2,V=V,famid=famid,prev=prev,data=data,n.cores=n.cores)
  } else {
    MLE.H0 <- LTMH.asc(model=new.model,init_beta=new.init_beta,init_h2=init_h2,V=V,famid=famid,prev=prev,data=data,n.cores=n.cores,proband=proband)
  }
  
  # Under H0 U H1
  if(as.character(model)[3]=='-1'){
    init_beta <- NULL
  } else {
    init_beta <- glm(model,data=data,family=binomial(link=probit))$coef
  }
  if(is.null(proband)){
    MLE.H01 <- LTMH(model=model,init_beta=init_beta,init_h2=init_h2,V=V,famid=famid,prev=prev,data=data,n.cores=n.cores)
  } else {
    MLE.H01 <- LTMH.asc(model=model,init_beta=init_beta,init_h2=init_h2,V=V,famid=famid,prev=prev,data=data,n.cores=n.cores,proband=proband)
  }
  
  Chisq <- -2*(MLE.H0$logL - MLE.H01$logL)
  Pvalue <- pchisq(Chisq,df=length(test.beta),lower.tail=F)
  Testing.res <- list(MLE_beta_std=MLE.H01$beta_std,MLE_beta_unstd=MLE.H01$beta_unstd,MLE_h2=MLE.H01$h2,Test.beta=test.beta,Chisq=Chisq,DF=length(test.beta),Pvalue=Pvalue)
  return(Testing.res)
}
