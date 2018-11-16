#############################################################################################
############### Conditional Expectation Score Test (CEST) Oct 17, 2018 ######################
#############################################################################################

########## H0 : h2=0
Testing.h2 <- function(famid,V,init_beta=NULL,prev,X=NULL,Y,n.cores=1){
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
  fin.res <- data.frame(Score=S.h2,var_Score_ver1=varS.h2.1,var_Score_ver2=varS.h2.2,Chisq=T.h2,Pvalue=pval)
  return(fin.res)
}


DoTest.h2 <- function(obs,fin.dat,totalfam,model,init_beta=NULL,prev,h2,n.cores=1,out){
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

	# Testing
	Testing.res <- Testing.h2(famid,V,init_beta,prev,X,Y,n.cores)
	write.table(Testing.res,out,col.names=F,row.names=F,quote=F,append=TRUE)
}



########## H0 : beta=0
Testing.beta <- function(init_beta,init_h2,V,famid,prev,Y,X,full.X,test.beta.idx,max.iter=100,max.sub.iter=50,n.cores=1){
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
			
			if(epsilon<1e-5){
				beta.new <- beta.hat
				print(data.frame(h2.new,epsilon,n.iter))
				break
			} else {
				h2.old <- h2.new
				beta.old <- beta.hat
				
				while(1){
				  n.iter=n.iter+1
				  resAB <- mclapply(1:nrow(Y),getAB_ij,mc.cores=n.cores)
				  O1 <- Reduce('+',lapply(1:nrow(Y),function(ij) resAB[[ij]]$O1ij))
				  O2 <- Reduce('+',lapply(1:nrow(Y),function(ij) resAB[[ij]]$O2ij))
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
					break
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
					break
				} else {
					beta.old <- beta.new
				}
			}
		}
	}
	
	#### Now, we have MLE for h2 under the H0: beta=0.
	FAMID <- unique(famid)
	X <- full.X
	final.beta <- matrix(0,ncol(full.X),1)
	final.beta[-test.beta.idx,1] <- beta.new[,1]
	
	resAB <- mclapply(FAMID,getAB,h2.old=h2.new,beta.old=final.beta,famid=famid,mc.cores=n.cores)

	## Score under the H0: beta=0
	get.Score <- function(ii){
	  
	  fam.idx <- resAB[[ii]]$fam
	  Vi <- resAB[[ii]]$Vi
	  Si <- h2.new*Vi+(1-h2.new)*diag(length(fam.idx))
	  invSi <- solve(Si)
	  
	  Ci <- -invSi%*%(Vi-diag(length(fam.idx)))%*%invSi
	  Xi <- matrix(full.X[fam.idx,,drop=F],nrow=length(fam.idx))
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

	return(data.frame(Score=S.beta,var_Score_ver1=varS.beta.1,var_Score_ver2=varS.beta.2,Chisq=T.beta,Pvalue=pval))
}			



DoTest.beta <- function(obs,fin.dat,totalfam,model,init_beta,test.beta,prev,init_h2,n.cores=1,out,seed=F){
	print(obs) 

	# Sampling test data
  if(seed) set.seed(i)
  dataset <- fin.dat[fin.dat$FID%in%sample(unique(fin.dat$FID),totalfam),,drop=F]
	Y <- dataset[,as.character(model)[2],drop=F]
	full.X <- model.matrix(model,dataset)
	test.beta.idx <- which(colnames(full.X)%in%test.beta)
	X <- full.X[,-test.beta.idx,drop=F]
	
  famid <- dataset$FID
  total_ped <- with(dataset,pedigree(id=IID,dadid=PID,momid=MID,sex=SEX,famid=FID,missid='0'))
  V <- 2*as.matrix(kinship(total_ped))

	# Testing
	Testing.res <- Testing.beta(init_beta,init_h2,V,famid,prev,Y,X,full.X,test.beta.idx,max.iter=100,max.sub.iter=50,n.cores)
	write.table(Testing.res,out,col.names=F,row.names=F,quote=F,append=TRUE)
}


##### Summarize : QQ plot & power table
get.Summary <- function(IN,OUT,prev,h2,totalfam,alpha,type=c("size","power"),dist=c("mixture","chisq")){
  ## Table
  dat <- read.table(IN,head=T)
  # pval <- pchisq(dat[,'Chisq'],df=1,lower.tail=F)/2
  pval <- dat[,'Pvalue']
  res <- lapply(alpha,function(i) data.frame(pvalue=sum(pval<i)/length(pval)))
  res.I <- do.call(rbind,res)
  rownames(res.I) <- alpha
  colnames(res.I) <- "Proportion"
  write.table(res.I,paste0(OUT,"_power.txt"),quote=F)
  
  if(type=="size"){
    ## QQ plots
    p.val <- pval
    png(paste0(OUT,"_QQ.png"), width=1000, height=1000, units = "px", bg = "white", res = 200)

    y <- -log(p.val,10)
    o.y <- sort(y[is.finite(y)],decreasing=T)
    if(dist=="mixture"){
      xx<- (1:length(o.y))/length(o.y)/2
    } else {
      xx<- (1:length(o.y))/length(o.y)
    }
    x <- -log(xx,10)
    ifelse(max(o.y[is.finite(o.y)])>=x[1],YY<-o.y,YY<-x)
    
    plot(YY, YY, main="Empirical Size",type = "n", xlab = expression(paste("Expected ",-log[10],"(p-value)")), ylab = expression(paste("Observed ",-log[10],"(p-value)")))
    N=length(o.y)
    c95 <- rep(0,N)
    c05 <- rep(0,N)
    for(i in 1:N){
      c95[i] <- ifelse(dist=="mixture",qbeta(0.95,i,N-i+1)/2,qbeta(0.95,i,N-i+1))
      c05[i] <- ifelse(dist=="mixture",qbeta(0.05,i,N-i+1)/2,qbeta(0.05,i,N-i+1))
    }
    abline(h = c(0:max(max(x),max(o.y[is.finite(o.y)]))), v =c(0:max(max(x),max(o.y[is.finite(o.y)]))), col = "darkgray", lty=3)
    polygon(c(x,sort(x)),c(-log(c95,10),sort(-log(c05,10))),col=c("gray"),border=NA)
    abline(a=0,b=1,col="black", lty=1)
    points(x, o.y, cex = .5, col = "dark red")

    dev.off()
  }
  return(res.I)
}


#####################################
### CEST for beta 
#####################################

DoCEST.beta.WGS <- function(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx=NULL,out){
  # chr : chromosome
  # model : formula except SNP
  # FID : Variable name for family ID in phenotype file
  # FID : Variable name for individual ID in phenotype file
  # working_dir : working directory containing raw data splitted into chromosomes.
  # phe : phenotype file. It should be same order of individuals with that for the raw file.
  # VV : Relationship matrix (matrix)
  # V.name : individual names correspond to the each row (vector)
  # prev : prevalence (scalar)
  # init_beta : initial value for beta except SNP. It should be p x 1 matrix.
  # init_h2 : initial value for h2. (scalar)
  # n.cores : the number of cores which are used for parallel analysis.
  # valid.idx : row index to be included for the analysis
  # out : output path
  
  setwd(working_dir)
  raw <- read.table(paste0("chr",chr,"/chr",chr,".raw"),head=T,stringsAsFactor=F)
  if(is.null(valid.idx)) valid.idx <- 1:nrow(phe)
  new.model <- formula(paste0(c(as.character(model)[c(2,1,3)],' + stdSNP'),collapse=' '))
  
  DoCEST.beta.SNP <- function(ii,phe,valid.idx,new.model,FID,IID,VV,V.name,prev,init_beta,init_h2,n.cores,out){
    print(ii)
    dataset <- cbind(phe,SNP=raw[,ii+6])[valid.idx,]
    dataset <- dataset[!is.na(dataset$SNP),]
    dataset$stdSNP <- (dataset$SNP-mean(dataset$SNP))/sd(dataset$SNP)
    test.beta <- 'stdSNP'
    
    Y <- dataset[,as.character(new.model)[2],drop=F]
    full.X <- model.matrix(new.model,dataset)
    test.beta.idx <- which(colnames(full.X)%in%test.beta)
    X <- full.X[,-test.beta.idx,drop=F]
    famid <- dataset[,FID]
    ord.idx <- match(dataset[,IID],V.name)
    V <- VV[ord.idx,ord.idx]
    colnames(V) <- NULL
    
    # Testing
    Testing.res <- Testing.beta(init_beta,init_h2,V,famid,prev,Y,X,full.X,test.beta.idx,max.iter=100,max.sub.iter=50,n.cores)
    res <- cbind(CHR=chr,SNP=gsub('_[ACGT]','',colnames(raw)[ii+6]),Testing.res)
    write.table(res,out,col.names=F,row.names=F,quote=F,append=TRUE)
  }
  
  write.table(data.frame('CHR','SNP','Score','var_Score_ver1','var_Score_ver2','Chisq','Pvalue'),out,col.names=F,row.names=F,quote=F)
  res <- lapply(1:(ncol(raw)-6),DoCEST.beta.SNP,phe=phe,valid.idx=valid.idx,new.model=new.model,FID=FID,IID=IID,VV=VV,V.name=V.name,prev=prev,init_beta=init_beta,init_h2=init_h2,n.cores=n.cores,out=out)
  return(NULL)
}


############# Temp function for already runned chr
tmp.DoCEST.beta.WGS <- function(chr,model,FID,IID,working_dir,phe,VV,V.name,prev,init_beta,init_h2,n.cores,valid.idx=NULL,out){
  # chr : chromosome
  # model : formula except SNP
  # FID : Variable name for family ID in phenotype file
  # FID : Variable name for individual ID in phenotype file
  # working_dir : working directory containing raw data splitted into chromosomes.
  # phe : phenotype file. It should be same order of individuals with that for the raw file.
  # VV : Relationship matrix (matrix)
  # V.name : individual names correspond to the each row (vector)
  # prev : prevalence (scalar)
  # init_beta : initial value for beta except SNP. It should be p x 1 matrix.
  # init_h2 : initial value for h2. (scalar)
  # n.cores : the number of cores which are used for parallel analysis.
  # valid.idx : row index to be included for the analysis
  # out : output path
  
  setwd(working_dir)
  raw <- read.table(paste0("chr",chr,"/chr",chr,".raw"),head=T,stringsAsFactor=F)
  if(is.null(valid.idx)) valid.idx <- 1:nrow(phe)
  new.model <- formula(paste0(c(as.character(model)[c(2,1,3)],' + stdSNP'),collapse=' '))
  
  DoCEST.beta.SNP <- function(ii,phe,valid.idx,new.model,FID,IID,VV,V.name,prev,init_beta,init_h2,n.cores,out){
    print(ii)
    dataset <- cbind(phe,SNP=raw[,ii+6])[valid.idx,]
    dataset <- dataset[!is.na(dataset$SNP),]
    dataset$stdSNP <- (dataset$SNP-mean(dataset$SNP))/sd(dataset$SNP)
    test.beta <- 'stdSNP'
    
    Y <- dataset[,as.character(new.model)[2],drop=F]
    full.X <- model.matrix(new.model,dataset)
    test.beta.idx <- which(colnames(full.X)%in%test.beta)
    X <- full.X[,-test.beta.idx,drop=F]
    famid <- dataset[,FID]
    ord.idx <- match(dataset[,IID],V.name)
    V <- VV[ord.idx,ord.idx]
    colnames(V) <- NULL
    
    # Testing
    Testing.res <- Testing.beta(init_beta,init_h2,V,famid,prev,Y,X,full.X,test.beta.idx,max.iter=100,max.sub.iter=50,n.cores)
    res <- cbind(CHR=chr,SNP=gsub('_[ACGT]','',colnames(raw)[ii+6]),Testing.res)
    write.table(res,out,col.names=F,row.names=F,quote=F,append=TRUE)
  }
  
  existing.res <- read.table(out,head=T)
  res <- lapply((nrow(existing.res)+1):(ncol(raw)-6),DoCEST.beta.SNP,phe=phe,valid.idx=valid.idx,new.model=new.model,FID=FID,IID=IID,VV=VV,V.name=V.name,prev=prev,init_beta=init_beta,init_h2=init_h2,n.cores=n.cores,out=out)
  return(NULL)
}




#######################################################
################## Ascertained ver ####################
#######################################################

########## H0 : h2=0
Testing.h2.new <- function(famid,V,init_beta=NULL,prev,X=NULL,Y,n.cores=1,proband=NULL){
  library(parallel)
  library(tmvtnorm)
  
  ### X : N x p dim
  ### Y : N x 1 dim
  ### init_beta : P x 1 dim

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
  fin.res <- data.frame(Score=S.h2,var_Score_ver1=varS.h2.1,var_Score_ver2=varS.h2.2,Chisq=T.h2,Pvalue=pval)
  colnames(fin.res) <- c('Score','var_Score_ver1','var_Score_ver2','Chisq','Pvalue')
  return(fin.res)
}


DoTest.h2 <- function(obs,fin.dat,totalfam,model,init_beta=NULL,prev,h2,n.cores=1,out,seed=F){
	print(obs) 

	# Sampling test data
  if(seed) set.seed(obs)
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

	# Testing
	Testing.res <- Testing.h2(famid,V,init_beta,prev,X,Y,n.cores)
	write.table(Testing.res,out,col.names=F,row.names=F,quote=F,append=TRUE)
}



########## H0 : beta=0
Testing.beta <- function(init_beta,init_h2,V,famid,prev,Y,X,full.X,test.beta.idx,max.iter=100,max.sub.iter=50,n.cores=1){
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
			
			if(epsilon<1e-5){
				beta.new <- beta.hat
				print(data.frame(h2.new,epsilon,n.iter))
				break
			} else {
				h2.old <- h2.new
				beta.old <- beta.hat
				
				while(1){
				  n.iter=n.iter+1
				  resAB <- mclapply(1:nrow(Y),getAB_ij,mc.cores=n.cores)
				  O1 <- Reduce('+',lapply(1:nrow(Y),function(ij) resAB[[ij]]$O1ij))
				  O2 <- Reduce('+',lapply(1:nrow(Y),function(ij) resAB[[ij]]$O2ij))
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
					break
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
					break
				} else {
					beta.old <- beta.new
				}
			}
		}
	}
	
	#### Now, we have MLE for h2 under the H0: beta=0.
	FAMID <- unique(famid)
	X <- full.X
	final.beta <- matrix(0,ncol(full.X),1)
	final.beta[-test.beta.idx,1] <- beta.new[,1]
	
	resAB <- mclapply(FAMID,getAB,h2.old=h2.new,beta.old=final.beta,famid=famid,mc.cores=n.cores)

	## Score under the H0: beta=0
	get.Score <- function(ii){
	  
	  fam.idx <- resAB[[ii]]$fam
	  Vi <- resAB[[ii]]$Vi
	  Si <- h2.new*Vi+(1-h2.new)*diag(length(fam.idx))
	  invSi <- solve(Si)
	  
	  Ci <- -invSi%*%(Vi-diag(length(fam.idx)))%*%invSi
	  Xi <- matrix(full.X[fam.idx,,drop=F],nrow=length(fam.idx))
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

	return(data.frame(Score=S.beta,var_Score_ver1=varS.beta.1,var_Score_ver2=varS.beta.2,Chisq=T.beta,Pvalue=pval))
}			



