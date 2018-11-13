kare_f <- function(i){
	rela <- kare[i,which(!is.na(kare[i,]))]

	num3 <- sum(rela==3)
	num4 <- sum(rela==4)
	ind_fam <- kare_fam[i,]

	if(!num3==0 & !num4==0){
		rela_fam <- matrix(0,4+num3,6)
		rela_fam[,5] <- 2	#set to female
		rela_fam[c(1,3),5] <- 1	#set to male for father
		rela_fam[,6] <- -9	#phenotype is set to missing
		rela_fam[,1] <- kare_fam[i,1]	#famid
		rela_fam[,2] <- paste(kare_fam[i,1],seq(nrow(rela_fam)),sep="_")	#individual id

		rela_fam[3,3] <- rela_fam[1,2]	#grand parent id
		rela_fam[3,4] <- rela_fam[2,2]	#grand parent id
		rela_fam[5:nrow(rela_fam),3] <- rela_fam[3,2]	#parent id
		rela_fam[5:nrow(rela_fam),4] <- rela_fam[4,2]	#parent id
		rela_fam[c(1:num4,5:nrow(rela_fam)),6] <- 2	#set to being affected
		
		if(any(rela==1)) rela_fam[3,6] <- 2	#set to being affected
		if(any(rela==2)) rela_fam[4,6] <- 2	#set to being affected
	}else if(!num3==0){
		rela_fam <- matrix(0,2+num3,6)
		rela_fam[,5] <- 2	#set to female
		rela_fam[1,5] <- 1	#set to male for father
		rela_fam[,6] <- -9	#phenotype is set to missing
		rela_fam[,1] <- kare_fam[i,1]	#famid
		rela_fam[,2] <- paste(kare_fam[i,1],seq(nrow(rela_fam)),sep="_")	#individual id

		rela_fam[3:nrow(rela_fam),3] <- rela_fam[1,2]	#parent id
		rela_fam[3:nrow(rela_fam),4] <- rela_fam[2,2]	#parent id
		rela_fam[c(3:nrow(rela_fam)),6] <- 2	#set to being affected
		
		if(any(rela==1)) rela_fam[1,6] <- 2	#set to being affected
		if(any(rela==2)) rela_fam[2,6] <- 2	#set to being affected
	}else if(!num4==0){
		rela_fam <- matrix(0,4,6)
		rela_fam[,5] <- 2
		rela_fam[c(1,3),5] <- 1
		rela_fam[,6] <- -9
		rela_fam[,1] <- kare_fam[i,1]
		rela_fam[,2] <- paste(kare_fam[i,1],seq(nrow(rela_fam)),sep="_")
	
		rela_fam[3,3] <- rela_fam[1,2]	#grand parent id
		rela_fam[3,4] <- rela_fam[2,2]	#grand parent id
		rela_fam[c(1:num4),6] <- 2

		if(any(rela==1)) rela_fam[3,6] <- 2	#set to being affected
		if(any(rela==2)) rela_fam[4,6] <- 2	#set to being affected
	}else{
		rela_fam <- matrix(0,2,6)
 		rela_fam[,5] <- 2
		rela_fam[1,5] <- 1
		rela_fam[,6] <- -9
		rela_fam[,1] <- kare_fam[i,1]
		rela_fam[,2] <- paste(kare_fam[i,1],seq(nrow(rela_fam)),sep="_")

		if(any(rela==1)) rela_fam[1,6] <- 2	#set to being affected
		if(any(rela==2)) rela_fam[2,6] <- 2	#set to being affected
	}
	kare_fam[i,] <- ind_fam
	return(rela_fam)
}


snu_f <- function(j) {
	res[[j]] <- which(snu[,1]==famid[j,1])
	fam_pheno <- snu[res[[j]],]
	rela_snu <- fam_pheno[fam_pheno$GWAS.ID=="",]
	rela_snu <- rela_snu[order(rela_snu$RELATIONSHIP),]

	rela <- rela_snu$RELATIONSHIP
	num4 <- sum(rela==4)
	num6 <- sum(rela==6)
	num7 <- sum(rela==7)
	
	if(num6==1){
		rela_fam <- matrix(0,4+num4+num6+num7,9)
		
		#sex
		rela_fam[,5] <- 2
		rela_fam[c(1,3),5] <- 1
		rela_fam[5,5] <- rela_snu[rela_snu$RELATIONSHIP==6,"SEX"]
		if(num4 > 0) rela_fam[(1:num4)+5,5] <- rela_snu[rela_snu$RELATIONSHIP==4,"SEX"]
		if(num7 > 0) rela_fam[(1:num7)+5+num4,5] <- rela_snu[rela_snu$RELATIONSHIP==7,"SEX"]
		
		# age
		rela_fam[,9] <- NA
		rela_fam[5,9] <- rela_snu[rela_snu$RELATIONSHIP==6,"AGE"]
		if(num4 > 0) rela_fam[(1:num4)+5,9] <- rela_snu[rela_snu$RELATIONSHIP==4,"AGE"]
		if(num7 > 0) rela_fam[(1:num7)+5+num4,9] <- rela_snu[rela_snu$RELATIONSHIP==7,"AGE"]

		#phenotype
		rela_fam[,6] <- -9
		if(any(rela==1)){
			grand <- rela_snu[rela_snu$RELATIONSHIP==1,]
			for(ii in 1:nrow(grand)){
				if(grand[ii,"SEX"]==1){
					rela_fam[1,6] <- 1
					rela_fam[1,6] <- grand[ii,"T2DM_IGT"]
				} else {
					rela_fam[2,6] <- 1
					rela_fam[2,6] <- grand[ii,"T2DM_IGT"]
				}
			}
		}
		if(any(rela==2)){
			parent <- rela_snu[rela_snu$RELATIONSHIP==2,]
			for(ii in 1:nrow(parent)){
				if(parent[ii,"SEX"]==1){
					rela_fam[3,6] <- 1
					rela_fam[3,6] <- parent[ii,"T2DM_IGT"]
				} else {
					rela_fam[4,6] <- 1
					rela_fam[4,6] <- parent[ii,"T2DM_IGT"]
				}
			}
		}
		rela_fam[5,6] <- rela_snu[rela_snu$RELATIONSHIP==6,"T2DM_IGT"]
		if(num4 > 0) rela_fam[(1:num4)+5,6] <- rela_snu[rela_snu$RELATIONSHIP==4,"T2DM_IGT"]
		if(num7 > 0) rela_fam[(1:num7)+5+num4,6] <- rela_snu[rela_snu$RELATIONSHIP==7,"T2DM_IGT"]

		#ID
		rela_fam[,1] <- fam_pheno[1,1]	#FID
		rela_fam[,2] <- paste(fam_pheno[1,1],1:nrow(rela_fam),sep="_")
		rela_fam[3,3] <- rela_fam[1,2]
		rela_fam[3,4] <- rela_fam[2,2]
		if(num4 > 0) {
			rela_fam[(1:num4)+5,3] <- rela_fam[3,2]
			rela_fam[(1:num4)+5,4] <- rela_fam[4,2]
		}
		if(num7 > 0) {
			if(fam_pheno[fam_pheno$RELATIONSHIP==5,"SEX"]==1){
				rela_fam[(1:num7)+5+num4,3] <- fam_pheno[fam_pheno$RELATIONSHIP==5,"GWAS.ID"]
				rela_fam[(1:num7)+5+num4,4] <- rela_fam[5,2]
			} else {
				rela_fam[(1:num7)+5+num4,3] <- rela_fam[5,2]
				rela_fam[(1:num7)+5+num4,4] <- fam_pheno[fam_pheno$RELATIONSHIP==5,"GWAS.ID"]
			}
		}

		#relationship
		rela_fam[1:2,7] <- 1
		rela_fam[3:4,7] <- 2
		rela_fam[5,7] <- 6
		if(num4 > 0) rela_fam[(1:num4)+5,7] <- 4
		if(num7 > 0) rela_fam[(1:num7)+5+num4,7] <- 7

		return(rela_fam)
	} else {
		rela_fam <- matrix(0,5+num4+num6+num7,9)

		#sex
		rela_fam[,5] <- 2
		rela_fam[c(1,3),5] <- 1
		if(num4 > 0) rela_fam[(1:num4)+5,5] <- rela_snu[rela_snu$RELATIONSHIP==4,"SEX"]
		if(fam_pheno[fam_pheno$RELATIONSHIP==5,"SEX"]==2) rela_fam[5,5] <- 1
		if(num7 > 0) rela_fam[(1:num7)+5+num4,5] <- rela_snu[rela_snu$RELATIONSHIP==7,"SEX"]

		# age
		rela_fam[,9] <- NA
		if(num4 > 0) rela_fam[(1:num4)+5,9] <- rela_snu[rela_snu$RELATIONSHIP==4,"AGE"]
		if(num7 > 0) rela_fam[(1:num7)+5+num4,9] <- rela_snu[rela_snu$RELATIONSHIP==7,"AGE"]
		
		#phenotype
		rela_fam[,6] <- -9
		if(any(rela==1)){
			grand <- rela_snu[rela_snu$RELATIONSHIP==1,]
			for(ii in 1:nrow(grand)){
				if(grand[ii,"SEX"]==1){
					rela_fam[1,6] <- 1
					rela_fam[1,6] <- grand[ii,"T2DM_IGT"]
				} else {
					rela_fam[2,6] <- 1
					rela_fam[2,6] <- grand[ii,"T2DM_IGT"]
				}
			}
		}
		if(any(rela==2)){
			parent <- rela_snu[rela_snu$RELATIONSHIP==2,]
			for(ii in 1:nrow(parent)){
				if(parent[ii,"SEX"]==1){
					rela_fam[3,6] <- 1
					rela_fam[3,6] <- parent[ii,"T2DM_IGT"]
				} else {
					rela_fam[4,6] <- 1
					rela_fam[4,6] <- parent[ii,"T2DM_IGT"]
				}
			}
		}
		if(num4 > 0) rela_fam[(1:num4)+5,6] <- rela_snu[rela_snu$RELATIONSHIP==4,"T2DM_IGT"]
		if(num7 > 0) rela_fam[(1:num7)+5+num4,6] <- rela_snu[rela_snu$RELATIONSHIP==7,"T2DM_IGT"]


		#ID
		rela_fam[,1] <- fam_pheno[1,1]	#FID
		rela_fam[,2] <- paste(fam_pheno[1,1],1:nrow(rela_fam),sep="_")
		rela_fam[3,3] <- rela_fam[1,2]
		rela_fam[3,4] <- rela_fam[2,2]
		if(num4 > 0) {
			rela_fam[(1:num4)+5,3] <- rela_fam[3,2]
			rela_fam[(1:num4)+5,4] <- rela_fam[4,2]
		}
		if(num7 > 0) {
			if(fam_pheno[fam_pheno$RELATIONSHIP==5,"SEX"]==1){
				rela_fam[(1:num7)+5+num4,3] <- fam_pheno[fam_pheno$RELATIONSHIP==5,"GWAS.ID"]
				rela_fam[(1:num7)+5+num4,4] <- rela_fam[5,2]
			} else {
				rela_fam[(1:num7)+5+num4,3] <- rela_fam[5,2]
				rela_fam[(1:num7)+5+num4,4] <- fam_pheno[fam_pheno$RELATIONSHIP==5,"GWAS.ID"]
			}
		}

		#relationship
		rela_fam[1:2,7] <- 1
		rela_fam[3:4,7] <- 2
		rela_fam[5,7] <- 6
		if(num4 > 0) rela_fam[(1:num4)+5,7] <- 4
		if(num7 > 0) rela_fam[(1:num7)+5+num4,7] <- 7

		return(rela_fam)
	}
}
