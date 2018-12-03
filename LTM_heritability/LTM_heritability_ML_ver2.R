########### MLE for h2, beta ###########

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

	Y <- data[,as.character(model)[2]]
	X <- model.matrix(model,data)
	beta.old <- matrix(init_beta,nrow=length(init_beta))
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
				return(list(beta=beta.new,h2=h2.new,n_iter=n.iter))
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
				print(data.frame(h2.new,epsilon,n.iter))
				
				if(epsilon<1e-5){
				  beta.new <- beta.hat
					return(list(beta=beta.new,h2=h2.new,n_iter=n.iter))
				} else {
					h2.old <- h2.new
					beta.old <- beta.hat
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
						break
					} else {
						h2.old <- h2.new
					}
				}

				epsilon <- sqrt(sum((beta.old-beta.new)^2))
				print(data.frame(h2.new,epsilon,n.iter))
				if(epsilon<1e-5){
					return(list(beta=beta.new,h2=h2.new,n_iter=n.iter))
				} else {
					beta.old <- beta.new
				}
			}
		}
	}
}			

########################## Ascertainment adjust LTMH

LTMH.asc <- function(model,init_beta,init_h2,V,famid,prev,dataset,max.iter=100,max.sub.iter=50,n.cores=1,proband,SAVE=F,out=NULL){
	## This function depends on a package; tmvtnorm, parallel, mvtnorm
	## model : Y~X 
	## init_beta : initial value for beta (p+1 vector)
	## init_h2 :  initial value for heritability
	## V : genetic relationship matrix (theoritically, kinship coefficient matrix)
	## famid : family id
	## prev : disease prevalence
	## dataset : data.frame containing variables specified at formula
	## proband : proband indicator (1 = proband, 0 = non-proband)

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
	
	tr <- function(m) sum(diag(m))

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

	Y <- dataset[,as.character(model)[2]]
	X <- model.matrix(model,dataset)
	beta.old <- beta.OLD <- matrix(init_beta,nrow=length(init_beta))
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
	    return(list(beta=beta.old,h2=h2.old,n_iter=n.iter))
	    break
	  }
	}
}




##### get max logL
getlogL <- function(model,V,famid,prev,dataset,max.iter=100,max.sub.iter=50,n.cores=1){
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
	n.point <- seq(0,1,by=0.01)
	Y.point <- do.call(rbind,mclapply(n.point,logL,mc.cores=n.cores))
	colnames(Y.point) <- c("n.point","Y.point")
	return(Y.point)
}


getmaxlogL <- function(fin.dat,totalfam,assumed_prev,model,n.cores=1){
	print(paste0("number of families : ",totalfam))
	famlist <- unique(fin.dat$FID)
	target <- sample(famlist,totalfam)
	dataset = fin.dat[fin.dat$FID%in%target,]
	total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
	V <- 2*as.matrix(kinship(total_ped))
	famid <- as.character(dataset$FID)
	output <- getlogL(model,V,famid,assumed_prev,dataset,max.iter=100,max.sub.iter=50,n.cores)
	maxlogL <- with(data.frame(output),n.point[Y.point==max(Y.point)])
	return(data.frame(numfam=totalfam,maxlogL=maxlogL))
}




getEsth2 <- function(i,fin.dat,init_beta,init_h2,totalfam,assumed_prev,model,n.cores=1,seed=F){
	print(i)
	famlist <- unique(fin.dat$FID)
	if(seed) set.seed(i)
	target <- sample(famlist,totalfam)
	dataset = fin.dat[fin.dat$FID%in%target,]
	total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
	V <- 2*as.matrix(kinship(total_ped))
	famid <- as.character(dataset$FID)
	output <- LTMH(model=model,init_beta=init_beta,init_h2=init_h2,V=V,famid=famid,prev=assumed_prev,data=dataset,n.cores=n.cores)
	write.table(t(as.matrix(c(obs=i,unlist(output)))),paste0("prev_",assumed_prev,"_h2_",init_h2,"/Esth2_",totalfam,".txt"),col.names=F,row.names=F,quote=F,append=TRUE)
	return(c(i,output))
}


getEsth2_ver2 <- function(i,fin.dat,init_beta,init_h2,totalfam,assumed_prev,model,n.cores=1,out_path){
	print(i)
	famlist <- unique(fin.dat$FID)
	target <- sample(famlist,totalfam)
	dataset = fin.dat[fin.dat$FID%in%target,]
	total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
	V <- 2*as.matrix(kinship(total_ped))
	famid <- as.character(dataset$FID)
	output <- LTMH(model=model,init_beta=init_beta,init_h2=init_h2,V=V,famid=famid,prev=assumed_prev,data=dataset,n.cores=n.cores)
	fin.res <- data.frame(obs=i,unlist(output))
	write.table(fin.res,outpath,col.names=F,row.names=F,quote=F,append=TRUE)
	return(c(i,output))
}


getEsth2.out <- function(i,fin.dat,init_beta,init_h2,totalfam,assumed_prev,model,n.cores=1,out,seed=F){
  print(i)
  famlist <- unique(fin.dat$FID)
  if(seed) set.seed(i)
  target <- sample(famlist,totalfam)
  dataset = fin.dat[fin.dat$FID%in%target,]
  total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
  V <- 2*as.matrix(kinship(total_ped))
  famid <- as.character(dataset$FID)
  output <- LTMH(model=model,init_beta=init_beta,init_h2=init_h2,V=V,famid=famid,prev=assumed_prev,data=dataset,n.cores=n.cores)
  write.table(t(as.matrix(c(obs=i,unlist(output)))),out,col.names=F,row.names=F,quote=F,append=TRUE)
  return(c(i,output))
}

getEsth2.GCTA <- function(i,fin.dat,totalfam,prev,res_var,exp_var,out,working_dir,seed=F){
  library(gap)
  library(kinship2)
  print(i)
  setwd(working_dir)
  famlist <- unique(fin.dat$FID)
  if(seed) set.seed(i)
  target <- sample(famlist,totalfam)
  dataset = fin.dat[fin.dat$FID%in%target,]
  total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
  V <- 2*as.matrix(kinship(total_ped))
  grm <- V[upper.tri(V,diag=T)]
  WriteGRMBin(paste0('tmp_',i,'_grm'),grm,nrow(dataset)*(nrow(dataset)+1)/2,cbind(dataset[,'FID'],dataset[,'IID']))
  write.table(dataset[,c('FID','IID',res_var)],paste0('tmp_',i,'_phen.txt'),row.names=F,col.names=F,quote=F)
  write.table(dataset[,c('FID','IID',exp_var)],paste0('tmp_',i,'_qcovar.txt'),row.names=F,col.names=F,quote=F)
  
  system(paste0('gcta64 --grm tmp_',i,'_grm --reml --pheno tmp_',i,'_phen.txt --qcovar tmp_',i,'_qcovar.txt --prevalence ',prev,' --out tmp_',i,'_output'))
  res <- system(paste0('grep V\\(G\\)\\/Vp_L tmp_',i,'_output.hsq'),intern=T)
  res.I <- cbind(i,matrix(strsplit(res,'\t')[[1]],nrow=1)[,-1,drop=F])
  colnames(res.I) <- c('Obs','h2','SE_h2')
  print(res.I)
  write.table(res.I,out,col.names=F,row.names=F,quote=F,append=TRUE)
  
  system(paste0('rm tmp_',i,'_*'))
}

getLRT.GCTA <- function(i,fin.dat,totalfam,prev,res_var,exp_var,out,working_dir,seed=F){
  library(gap)
  library(kinship2)
  print(i)
  setwd(working_dir)
  famlist <- unique(fin.dat$FID)
  if(seed) set.seed(i)
  target <- sample(famlist,totalfam)
  dataset = fin.dat[fin.dat$FID%in%target,]
  total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
  V <- 2*as.matrix(kinship(total_ped))
  grm <- V[upper.tri(V,diag=T)]
  WriteGRMBin(paste0('LRT_',i,'_grm'),grm,nrow(dataset)*(nrow(dataset)+1)/2,cbind(dataset[,'FID'],dataset[,'IID']))
  write.table(dataset[,c('FID','IID',res_var)],paste0('LRT_',i,'_phen.txt'),row.names=F,col.names=F,quote=F)
  write.table(dataset[,c('FID','IID',exp_var)],paste0('LRT_',i,'_qcovar.txt'),row.names=F,col.names=F,quote=F)
  
  system(paste0('gcta64 --grm LRT_',i,'_grm --reml --pheno LRT_',i,'_phen.txt --qcovar LRT_',i,'_qcovar.txt --prevalence ',prev,' --out LRT_',i,'_output'))
  res <- system(paste0('grep -E "LRT|df|Pval" LRT_',i,'_output.hsq'),intern=T)
  res <- matrix(gsub('.+\t','',res),nrow=1)
  res.I <- cbind(i,res)
  colnames(res.I) <- c('Obs','LRT','DF','Pvalue')
  print(res.I)
  write.table(res.I,out,col.names=F,row.names=F,quote=F,append=TRUE)
  
  system(paste0('rm LRT_',i,'_*'))
}


getEsth2.asc <- function(i,fin.dat,init_beta,init_h2,totalfam,assumed_prev,model,n.cores=1,out,seed=F,PB){
  print(i)
  famlist <- unique(fin.dat$FID)
  if(seed) set.seed(i)
  target <- sample(famlist,totalfam)
  dataset = fin.dat[fin.dat$FID%in%target,]
  proband <- dataset[,PB]
  total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
  V <- 2*as.matrix(kinship(total_ped))
  famid <- as.character(dataset$FID)
  output <- LTMH.asc(model=model,init_beta=init_beta,init_h2=init_h2,V=V,famid=famid,prev=assumed_prev,dataset=dataset,n.cores=n.cores,proband=proband)
  write.table(t(as.matrix(c(obs=i,unlist(output)))),out,col.names=F,row.names=F,quote=F,append=TRUE)
  return(c(i,output))
}
