#### data.generate
### Nuclear family
genNucFam <- function(totalfam,MAF,h2,ha2,prev,num_snp){
	library(MASS)
	## nucelar family ##
	numfam <- sample(3:6,totalfam,replace=T,prob=c(0.2,0.3,0.3,0.2))
	famid 	<- fatherid <- motherid <- famsize <- indiv <- iid <- pid <- mid <- sex <- c()
	genFam <- function(i,totalfam,numfam){
		famid <- rep(paste("fam",i,sep=""),numfam[i])
		IID <- paste("fam",i,"_",1:numfam[i],sep="")
		temp1 <- temp2 <- ind <- PID <- MID <- rep(0,numfam[i])
		PID[3:numfam[i]] <- IID[1]
		MID[3:numfam[i]] <- IID[2]
		SEX <- c(1,2,rep(1,numfam[i]-2))
		temp1[1] <- temp2[2] <-1
		famsize	<- rep(numfam[i],numfam[i])
		ind[3] 	<- 1
		res <- data.frame(famid,father=temp1,mother=temp2,famsize=famsize,individual=ind,FID=famid,IID,PID,MID,SEX,ind)
		return(res)
	}
	fam_res <- lapply(1:length(numfam),genFam,totalfam=totalfam,numfam=numfam)
	fam_res_I <- do.call(rbind,fam_res)
	dat <- fam_res_I[,1:5]
	fam <- fam_res_I[,6:11]

	## Design X, beta
	Gen_offs <- function(j){						## function for family size j
		nfam_j <- sum(numfam==j)						## It is better to calculate once
		Geno  <- matrix(0,j,nfam_j)					## It is better to define the required matrix size
		Geno[1:2,] <- matrix(rbinom(nfam_j*2,2,MAF),nrow=2)
		for(jj in 1:(j-2)){
			o1 <- which(Geno[1,]/2 > runif(nfam_j))
			o2 <- which(Geno[2,]/2 > runif(nfam_j))
			P1 <- P2 <- rep(0, nfam_j)
			P1[o1] <- P2[o2] <- 1
			off <- P1+P2
			Geno[(jj+2),] <- off
		}
		attr(Geno,'dim') <- NULL
		return(Geno)							## attr(Geno,'dim')<-NULL is better in this case.
	}

	DesignX <- function(i){
		X <- rep(0,nrow(dat))						## It is better to define the matrix clearly.
		X[which(dat$famsize==3)] <- Gen_offs(3)
		X[which(dat$famsize==4)] <- Gen_offs(4)
		X[which(dat$famsize==5)] <- Gen_offs(5)
		X[which(dat$famsize==6)] <- Gen_offs(6)
		return(X)
	}

	geno <- as.matrix(DesignX(1),ncol=1)

	## Generating Liability ##
	Gen_L <- function(i){
		phe <- phe_1 <- rep(0,nrow(dat))
		poly <- error <- c()
		for(j in 3:6){
			oo <- which(dat$famsize==j)
			V <- matrix(1/2,j,j)
			V[1,2] <- V[2,1] <- 0
			diag(V) <- 1

			pol <- mvrnorm(n=sum(numfam==j),mu=rep(0,j),Sigma=h2*V)
			err <- mvrnorm(n=sum(numfam==j),mu=rep(0,j),Sigma=(1-h2)*diag(j))
			pol <- t(pol); err <- t(err)
			attr(pol,'dim') <- NULL; attr(err,'dim') <- NULL

			poly[oo] <- pol; error[oo] <- err
		}
		liab <- geno[,i]*beta + poly + error
		liab_1 <- poly + error
		phe[which(liab>thres)] <- phe_1[which(liab_1>thres)] <- 1
		res <- matrix(c(liab,phe,liab_1,phe_1),ncol=4)
		res <- rbind(seq((4*i-3),4*i),res)
		return(res)
	}
	beta <- sqrt((ha2)/(2*MAF*(1-MAF)*(1-ha2)))
	thres <- qnorm(prev,2*beta*MAF,sqrt(2*beta^2*MAF*(1-MAF)+1),lower.tail=F)
	outII_Gen_L <- Gen_L(1)

	L <- outII_Gen_L[-1,order(outII_Gen_L[1,])]
	if(num_snp==0){
	  L_pow 	<- L[,3]
	  pheno_pow 	<- L[,4]
	} else {
	  L_pow 	<- L[,seq(1,num_snp*4,by=4)]
	  pheno_pow 	<- L[,seq(2,num_snp*4,by=4)]
	  
  }
	fin.dat <- data.frame(fam,Y=pheno_pow,L=L_pow,snp=geno)
	return(fin.dat)
}

