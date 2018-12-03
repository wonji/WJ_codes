#### Generate simulation data
# prev : 0.05,0.1,0.2 | h2 : 0.05,0.2,0.4
source('/home2/wjkim/paper/heritability/ML_ver2/variousFam/gen_simuldata.r')
for(prev in c(0.05,0.1,0.2)){
	for(h2 in c(0.05,0.2, 0.4)){
		print(paste0('prevalence:',prev,', heritability:',h2))
		setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/")
		system(paste0("mkdir prev_",prev,"_h2_",h2))
		dataset = genNucFam(totalfam=10000,MAF=0.2,h2,ha2=0.05,prev,num_snp=1)
		write.table(dataset,paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),row.names=F,quote=F)
	}
}




#### Convergence of heritability
## family size : 500
source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
num_snp=1;n.cores=20
for(prev in c(0.1,0.2)){
	for(h2 in c(0.2, 0.4)){
		for(totalfam in c(500)){
			print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
			library(kinship2)
			setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/")
			fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
			fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
			
			
			
			model <- Y~std_snp-1
			write.table(data.frame('obs','beta_std_snp','h2','n_iteration'),paste0("prev_",prev,"_h2_",h2,"/Esth2_",totalfam,".txt"),col.names=F,row.names=F,quote=F)
			Esth2 <- sapply(1:300,getEsth2,fin.dat=fin.dat,init_beta=1,init_h2=h2,totalfam=totalfam,assumed_prev=prev,model=model,n.cores=n.cores)
		}
	}
}


#### When prevalence is misspecified
## scenario 1 : assumped prev = true prev+0.1, 
source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
num_snp=1;n.cores=20
for(prev in c(0.1,0.2)){
	for(h2 in c(0.2, 0.4)){
		for(totalfam in c(500)){
			assumed_prev=prev+0.1
			print(paste0('true prevalence:',prev,', assumed prevalence:',assumed_prev,' heritability:',h2,', number of families:',totalfam))
			library(kinship2)

			setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam")
			system(paste0("mkdir misspecified/prev_",assumed_prev,"_h2_",h2))
			fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
			fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
			model <- Y~std_snp-1
			setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/misspecified")
			write.table(data.frame('obs','beta_std_snp','h2','n_iteration'),paste0("prev_",assumed_prev,"_h2_",h2,"/Esth2_",totalfam,".txt"),col.names=F,row.names=F,quote=F)
			Esth2 <- sapply(1:300,getEsth2,fin.dat=fin.dat,init_beta=1,init_h2=h2,totalfam=totalfam,assumed_prev=assumed_prev,model=model,n.cores=n.cores)
		}
	}
}



## robust version
source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
num_snp=1;n.cores=20
for(prev in c(0.1,0.2)){
	for(h2 in c(0.2, 0.4)){
		for(totalfam in c(500)){
			assumed_prev=prev+0.1
			print(paste0('true prevalence:',prev,', assumed prevalence:',assumed_prev,' heritability:',h2,', number of families:',totalfam))
			library(kinship2)

			setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam")
			system(paste0("mkdir robust_prev_intercept/prev_",assumed_prev,"_h2_",h2))
			fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
			fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
			model <- Y~std_snp
			setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/robust_prev_intercept")
			write.table(data.frame('obs','tau','beta_std_snp','h2','n_iteration'),paste0("prev_",assumed_prev,"_h2_",h2,"/Esth2_",totalfam,".txt"),col.names=F,row.names=F,quote=F)
			Esth2 <- sapply(1:300,getEsth2,fin.dat=fin.dat,init_beta=c(0,1),init_h2=h2,totalfam=totalfam,assumed_prev=assumed_prev,model=model,n.cores=n.cores)
		}
	}
}





## robust version
source('~/paper/heritability/ML_ver2/LTM_heritability_ML_robust.R')
num_snp=1;n.cores=20
for(prev in c(0.1,0.2)){
	for(h2 in c(0.2, 0.4)){
		for(totalfam in c(500)){
			assumed_prev=prev+0.1
			print(paste0('true prevalence:',prev,', assumed prevalence:',assumed_prev,' heritability:',h2,', number of families:',totalfam))
			library(kinship2)

			setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam")
			system(paste0("mkdir robust_prev/prev_",assumed_prev,"_h2_",h2))
			fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
			fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
			model <- Y~std_snp-1
			setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/robust_prev")
			write.table(data.frame('obs','beta_std_snp','h2','c','n_iteration'),paste0("prev_",assumed_prev,"_h2_",h2,"/Esth2_",totalfam,".txt"),col.names=F,row.names=F,quote=F)
			Esth2 <- sapply(1:300,getEsth2,fin.dat=fin.dat,init_beta=1,init_h2=h2,init_c=0,totalfam=totalfam,assumed_prev=assumed_prev,model=model,n.cores=n.cores)
		}
	}
}








setwd("C:\\Users\\rewki\\Documents\\paper\\my_paper\\heritability\\simulation\\well_specified\\prev_0.1_h2_0.2")
aa <- read.table("Esth2_500.txt",head=F)
hist(aa$V3,xlim=c(0,1),main='Prevalence:0.1 & Heritability:0.2',xlab='h2')

h2 = 0.2; prev=0.1
Esth2 <- aa
			x.point <- c(1:length(Esth2[,3]))
			Y.point <- cumsum(Esth2[,3])/x.point
			png(paste0("Esth2_",totalfam,".png"), width=1700, height=1000, units = "px", bg = "white", res = 200)
			plot(x.point,Y.point,xlab="Cumulative number of families",ylab="Sample mean of estmated heritabilities",main="Convergence of estimate",ylim=c(0,1),type="l")
			legend("topright",legend=paste0("True h2 : ",h2," & prev : ",prev))
			abline(h=h2,lty=2,col="dark gray")
			dev.off()


setwd("C:\\Users\\rewki\\Documents\\paper\\my_paper\\heritability\\simulation\\well_specified\\prev_0.1_h2_0.4")
aa <- read.table("Esth2_500.txt",head=F)
hist(aa$V3,xlim=c(0,1),main='Prevalence:0.1 & Heritability:0.4',xlab='h2')

h2 = 0.4; prev=0.1
Esth2 <- aa
			x.point <- c(1:length(Esth2[,3]))
			Y.point <- cumsum(Esth2[,3])/x.point
			png(paste0("Esth2_",totalfam,".png"), width=1700, height=1000, units = "px", bg = "white", res = 200)
			plot(x.point,Y.point,xlab="Cumulative number of families",ylab="Sample mean of estmated heritabilities",main="Convergence of estimate",ylim=c(0,1),type="l")
			legend("topright",legend=paste0("True h2 : ",h2," & prev : ",prev))
			abline(h=h2,lty=2,col="dark gray")
			dev.off()







source('~/paper/heritability/ML/LTM_heritability_ML.R')



#######################################################
############### Ascertained samples ###################
#######################################################

#### Ascertained samples - n14
## family size : 500
source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
num_snp=1;n.cores=32
for(prev in c(0.1,0.2)){
	for(h2 in c(0.2, 0.4)){
		for(totalfam in c(500)){
			print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
			library(kinship2)
			setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/")
			fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
			fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
			model <- Y~std_snp-1

			# families containing at least one affected member
			all.famlist <- unique(fin.dat$FID)
			num.aff <- sapply(all.famlist,function(i) sum(fin.dat$Y[fin.dat$FID==i]))
			fin.dat <- fin.dat[fin.dat$FID%in%all.famlist[num.aff>0],]
			out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/ascertained/prev_",prev,"_h2_",h2,"/Esth2_",totalfam,".txt")
			write.table(data.frame('obs','beta_std_snp','h2','n_iteration'),out,col.names=F,row.names=F,quote=F)

			Esth2 <- sapply(1:300,getEsth2.out,fin.dat=fin.dat,init_beta=1,init_h2=h2,totalfam=totalfam,assumed_prev=prev,model=model,n.cores=n.cores,out=out)
		}
	}
}

source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
num_snp=1;n.cores=32
for(prev in c(0.1,0.2)){
	for(h2 in c(0.2, 0.4)){
		for(totalfam in c(500)){
			print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
			library(kinship2)
			setwd("/home2/wjkim/paper/heritability/ML_ver2/variousFam/")
			fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
			fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
			model <- Y~std_snp-1

			# families containing at least one affected member
			all.famlist <- unique(fin.dat$FID)
			num.aff <- sapply(all.famlist,function(i) sum(fin.dat$Y[fin.dat$FID==i]))
			fin.dat <- fin.dat[fin.dat$FID%in%all.famlist[num.aff>0],]
			out <- paste0("/home2/wjkim/paper/heritability/ML_ver2/variousFam/ascertained/prev_",prev,"_h2_",h2,"/Esth2_",totalfam,"_181029.txt")
			write.table(data.frame('obs','beta_std_snp','h2','n_iteration'),out,col.names=F,row.names=F,quote=F)

			Esth2 <- sapply(1:300,getEsth2.out,fin.dat=fin.dat,init_beta=1,init_h2=h2,totalfam=totalfam,assumed_prev=prev,model=model,n.cores=n.cores,out=out)
		}
	}
}



##########################################################
######### New Estimating 
##########################################################

###### Data generation
source('~/paper/heritability/ML_ver2/variousFam/gen_simuldata.r')
for(prev in c(0.1,0.2)){
	for(h2 in c(0.2,0.4)){
		print(paste0('prevalence:',prev,', heritability:',h2))
		setwd("~/paper/heritability/ML_ver2/variousFam/")
		system(paste0("mkdir prev_",prev,"_h2_",h2))
		dataset = genNucFam(totalfam=10000,MAF=0.2,h2,ha2=0.05,prev,num_snp=1)
		write.table(dataset,paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),row.names=F,quote=F)
	}
}

source('~/paper/heritability/ML_ver2/variousFam/gen_simuldata.r')
for(prev in c(0.05)){
	for(h2 in c(0.05,0.2,0.4)){
		print(paste0('prevalence:',prev,', heritability:',h2))
		setwd("~/paper/heritability/ML_ver2/variousFam/")
		system(paste0("mkdir prev_",prev,"_h2_",h2))
		dataset = genNucFam(totalfam=50000,MAF=0.2,h2,ha2=0.05,prev,num_snp=1)
		write.table(dataset,paste0("prev_",prev,"_h2_",h2,"/dataset_add.txt"),row.names=F,quote=F)
	}
}


source('~/paper/heritability/ML_ver2/variousFam/gen_simuldata.r')
for(prev in c(0.1,0.2)){
	for(h2 in c(0.05)){
		print(paste0('prevalence:',prev,', heritability:',h2))
		setwd("~/paper/heritability/ML_ver2/variousFam/")
		system(paste0("mkdir prev_",prev,"_h2_",h2))
		dataset = genNucFam(totalfam=50000,MAF=0.2,h2,ha2=0.05,prev,num_snp=1)
		write.table(dataset,paste0("prev_",prev,"_h2_",h2,"/dataset_add.txt"),row.names=F,quote=F)
	}
}


# combine added dataset to existing dataset
for(prev in c(0.05)){
  for(h2 in c(0.05,0.2,0.4)){
    print(paste0('prevalence:',prev,', heritability:',h2))
    setwd("~/paper/heritability/ML_ver2/variousFam/")
    dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
    dat.add <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset_add.txt"),head=T,stringsAsFactor=F)
    new.dat <- rbind(dat,dat.add)
    write.table(new.dat,paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),row.names=F,quote=F)
  }
}


#### Convergence of heritability
## family size : 500
source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
num_snp=1;n.cores=32
for(prev in c(0.05,0.1,0.2)){
	for(h2 in c(0.05,0.2,0.4)){
		for(totalfam in c(500)){
			print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
			library(kinship2)
			setwd("~/paper/heritability/ML_ver2/variousFam/")
			fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset_add.txt"),head=T,stringsAsFactor=F)
			fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)

			model <- Y~std_snp-1
			out <- paste0("~/paper/heritability/ML_ver2/variousFam/prev_",prev,"_h2_",h2,"/Esth2_",totalfam,"_181030.txt")
			write.table(data.frame('obs','beta_std_snp','h2','n_iteration'),out,col.names=F,row.names=F,quote=F)
			Esth2 <- sapply(1:300,getEsth2.out,fin.dat=fin.dat,init_beta=1,init_h2=h2,totalfam=totalfam,assumed_prev=prev,model=model,n.cores=n.cores,out=out,seed=T)
		}
	}
}

##########################################################
######### Ascertained
##########################################################

#### Convergence of heritability
est.ft <- function(PREV,H2,N.fam,n.cores){
  source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
  for(prev in PREV){
  	for(h2 in H2){
  		for(totalfam in N.fam){
  			print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
  			library(kinship2)
  			setwd("~/paper/heritability/ML_ver2/variousFam/")
  			fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
  			fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
  			fin.dat <- fin.dat[fin.dat$FID%in%fin.dat$FID[fin.dat$ind==1 & fin.dat$Y==1],,drop=F]
  			PB <- 'ind'
  
  			model <- Y~std_snp-1
  			out <- paste0("~/paper/heritability/ML_ver2/variousFam/prev_",prev,"_h2_",h2,"/test.txt")
  			write.table(data.frame('obs','beta_std_snp','h2','n_iteration'),out,col.names=F,row.names=F,quote=F)
  			Esth2 <- sapply(1:300,getEsth2.asc,fin.dat=fin.dat,init_beta=0.1,init_h2=h2,totalfam=totalfam,assumed_prev=prev,model=model,n.cores=n.cores,out=out,seed=T,PB=PB)
  		}
  	}
  }
}

## n2, 20 cores, p 0.1 h2 0.05
est.ft(0.1,0.05,500,20)

## n14, 32 cores, p 0.1 h2 0.2
est.ft(0.1,0.2,500,32)

## n7, 24 cores, p 0.1 h2 0.4
est.ft(0.1,0.4,500,24)

## n10, 24 cores, p 0.2 h2 0.05
est.ft(0.2,0.05,500,24)

## n11, 24 cores, p 0.2 h2 0.2
est.ft(0.2,0.2,500,24)

## n4, 20 cores, p 0.2 h2 0.4
est.ft(0.2,0.4,500,20)

## n14, 32 cores, p 0.05 h2 0.05,0.2,0.4

for(h2 in c(0.05,0.2,0.4)) est.ft(0.05,h2,500,32)




#### new Convergence of heritability
est.ft <- function(PREV,H2,N.fam,n.cores){
  source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
  for(prev in PREV){
  	for(h2 in H2){
  		for(totalfam in N.fam){
  			print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
  			library(kinship2)
  			setwd("~/paper/heritability/ML_ver2/variousFam/")
  			fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
  			fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
  			fin.dat <- fin.dat[fin.dat$FID%in%fin.dat$FID[fin.dat$ind==1 & fin.dat$Y==1],,drop=F]
  			PB <- 'ind'
  
  			model <- Y~std_snp-1
  			out <- paste0("~/paper/heritability/ML_ver2/variousFam/prev_",prev,"_h2_",h2,"/Esth2_",totalfam,"_newasc.txt")
  			write.table(data.frame('obs','beta_std_snp','h2','n_iteration'),out,col.names=F,row.names=F,quote=F)
  			Esth2 <- sapply(301:600,getEsth2.asc,fin.dat=fin.dat,init_beta=0.1,init_h2=h2,totalfam=totalfam,assumed_prev=prev,model=model,n.cores=n.cores,out=out,seed=T,PB=PB)
  		}
  	}
  }
}

## n10, 24 cores, p 0.1 h2 0.05
est.ft(0.1,0.05,500,24)

## n11, 24 cores, p 0.2 h2 0.05
est.ft(0.2,0.05,500,24)





#########################################################
###### Increase the number of replicates to 2000 ########
#########################################################

###### Scenario 1
getLTMH.S1 <- function(PREV, H2, N.fam, n.cores){
  source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
  for(prev in PREV){
    for(h2 in H2){
      for(totalfam in N.fam){
        print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
        library(kinship2)
        setwd("~/paper/heritability/ML_ver2/variousFam/")
        fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
        fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
        
        model <- Y~std_snp-1
        out <- paste0("~/paper/heritability/ML_ver2/variousFam/prev_",prev,"_h2_",h2,"/Esth2_",totalfam,"_181114.txt")
        if(paste0("Esth2_",totalfam,"_181114.txt") %in% system(paste0("ls prev_",prev,"_h2_",h2),intern=T)){
          start <- as.numeric(strsplit(system(paste0("tail -1 ",out),intern=T)," ")[[1]][1])+1
        } else {
          start <- 301
          write.table(data.frame('obs','beta_std_snp','h2','n_iteration'),out,col.names=F,row.names=F,quote=F)
        }
        Esth2 <- sapply(start:2000,getEsth2.out,fin.dat=fin.dat,init_beta=1,init_h2=h2,totalfam=totalfam,assumed_prev=prev,model=model,n.cores=n.cores,out=out,seed=T)
      }
    }
  }
}

## n2, 20 cores, p 0.05
a <- getLTMH.S1(PREV=0.05,H2=c(0.2),N.fam=500,n.cores=20)

## n4, 24 cores, p 0.1
a <- getLTMH.S1(PREV=0.1,H2=c(0.2),N.fam=500,n.cores=24)

## n14, 32 cores
a <- getLTMH.S1(PREV=0.05,H2=c(0.4),N.fam=500,n.cores=32)

## n5, 24 cores
a <- getLTMH.S1(PREV=0.1,H2=c(0.4),N.fam=500,n.cores=24)


##### Scenario 2
getLTMH.S2 <- function(PREV,H2,N.fam,n.cores){
  source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
  for(prev in PREV){
    for(h2 in H2){
      for(totalfam in N.fam){
        print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
        library(kinship2)
        setwd("~/paper/heritability/ML_ver2/variousFam/")
        fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
        fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
        fin.dat <- fin.dat[fin.dat$FID%in%fin.dat$FID[fin.dat$ind==1 & fin.dat$Y==1],,drop=F]
        PB <- 'ind'
        
        model <- Y~std_snp-1
        out <- paste0("~/paper/heritability/ML_ver2/variousFam/prev_",prev,"_h2_",h2,"/Esth2_",totalfam,"_asc_181114.txt")
        if(paste0("Esth2_",totalfam,"_asc_181114.txt") %in% system(paste0("ls prev_",prev,"_h2_",h2),intern=T)){
          start <- as.numeric(strsplit(system(paste0("tail -1 ",out),intern=T)," ")[[1]][1])+1
        } else {
          start <- 301
          write.table(data.frame('obs','beta_std_snp','h2','n_iteration'),out,col.names=F,row.names=F,quote=F)
        }
        Esth2 <- sapply(start:2000,getEsth2.asc,fin.dat=fin.dat,init_beta=0.1,init_h2=h2,totalfam=totalfam,assumed_prev=prev,model=model,n.cores=n.cores,out=out,seed=T,PB=PB)
      }
    }
  }
}

## n9, 24 cores, p 0.05
a <- getLTMH.S2(PREV=0.05,H2=c(0.4),N.fam=500,n.cores=24)

## n10, 24 cores, p 0.1
a <- getLTMH.S2(PREV=0.1,H2=c(0.4),N.fam=500,n.cores=24)


## n6, 24 cores, h 0.05
getLTMH.S2 <- function(PREV,H2,N.fam,n.cores){
  source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
  for(prev in PREV){
    for(h2 in H2){
      for(totalfam in N.fam){
        print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
        library(kinship2)
        setwd("~/paper/heritability/ML_ver2/variousFam/")
        fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
        fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
        fin.dat <- fin.dat[fin.dat$FID%in%fin.dat$FID[fin.dat$ind==1 & fin.dat$Y==1],,drop=F]
        PB <- 'ind'
        
        model <- Y~std_snp-1
        out <- paste0("~/paper/heritability/ML_ver2/variousFam/prev_",prev,"_h2_",h2,"/Esth2_",totalfam,"_asc_181114.txt")
        if(paste0("Esth2_",totalfam,"_asc_181114.txt") %in% system(paste0("ls prev_",prev,"_h2_",h2),intern=T)){
          start <- as.numeric(strsplit(system(paste0("tail -1 ",out),intern=T)," ")[[1]][1])+1
        } else {
          start <- ifelse(prev==0.05,301,601)
          write.table(data.frame('obs','beta_std_snp','h2','n_iteration'),out,col.names=F,row.names=F,quote=F)
        }
        Esth2 <- sapply(start:2000,getEsth2.asc,fin.dat=fin.dat,init_beta=0.1,init_h2=h2,totalfam=totalfam,assumed_prev=prev,model=model,n.cores=n.cores,out=out,seed=T,PB=PB)
      }
    }
  }
}

a <- getLTMH.S2(PREV=c(0.05),H2=0.05,N.fam=500,n.cores=24)

## n11, 19 cores
a <- getLTMH.S2(PREV=c(0.2),H2=0.05,N.fam=500,n.cores=19)

## n13, 32 cores, h 0.05
getLTMH.S2 <- function(PREV,H2,N.fam,n.cores){
  source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
  for(prev in PREV){
    for(h2 in H2){
      for(totalfam in N.fam){
        print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
        library(kinship2)
        setwd("~/paper/heritability/ML_ver2/variousFam/")
        fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset_add.txt"),head=T,stringsAsFactor=F)
        fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
        fin.dat <- fin.dat[fin.dat$FID%in%fin.dat$FID[fin.dat$ind==1 & fin.dat$Y==1],,drop=F]
        PB <- 'ind'
        
        model <- Y~std_snp-1
        out <- paste0("~/paper/heritability/ML_ver2/variousFam/prev_",prev,"_h2_",h2,"/Esth2_",totalfam,"_asc_181120.txt")
        if(paste0("Esth2_",totalfam,"_asc_181120.txt") %in% system(paste0("ls prev_",prev,"_h2_",h2),intern=T)){
          start <- as.numeric(strsplit(system(paste0("tail -1 ",out),intern=T)," ")[[1]][1])+1
        } else {
          start <- 1
          write.table(data.frame('obs','beta_std_snp','h2','n_iteration'),out,col.names=F,row.names=F,quote=F)
        }
        Esth2 <- sapply(start:1000,getEsth2.asc,fin.dat=fin.dat,init_beta=0.1,init_h2=h2,totalfam=totalfam,assumed_prev=prev,model=model,n.cores=n.cores,out=out,seed=T,PB=PB)
      }
    }
  }
}
a <- getLTMH.S2(PREV=c(0.1),H2=0.05,N.fam=500,n.cores=32)

## n8, 24 cores
getLTMH.S2 <- function(PREV,H2,N.fam,n.cores){
  source('~/paper/heritability/ML_ver2/LTM_heritability_ML_ver2.R')
  for(prev in PREV){
    for(h2 in H2){
      for(totalfam in N.fam){
        print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
        library(kinship2)
        setwd("~/paper/heritability/ML_ver2/variousFam/")
        fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset_add.txt"),head=T,stringsAsFactor=F)
        fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
        fin.dat <- fin.dat[fin.dat$FID%in%fin.dat$FID[fin.dat$ind==1 & fin.dat$Y==1],,drop=F]
        PB <- 'ind'
        
        model <- Y~std_snp-1
        out <- paste0("~/paper/heritability/ML_ver2/variousFam/prev_",prev,"_h2_",h2,"/Esth2_",totalfam,"_asc_181120.txt")
        Esth2 <- sapply(1001:2000,getEsth2.asc,fin.dat=fin.dat,init_beta=0.1,init_h2=h2,totalfam=totalfam,assumed_prev=prev,model=model,n.cores=n.cores,out=out,seed=T,PB=PB)
      }
    }
  }
}
a <- getLTMH.S2(PREV=c(0.1),H2=0.05,N.fam=500,n.cores=24)


##################################
###### new standardization #######
##################################
##### Scenario 2

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
  output <- LTMH.asc(model=model,init_beta=init_beta,init_h2=init_h2,V=V,famid=famid,prev=assumed_prev,data=dataset,n.cores=n.cores,proband=proband)
  write.table(t(as.matrix(c(obs=i,unlist(output)))),out,col.names=F,row.names=F,quote=F,append=TRUE)
  return(c(i,output))
}

LTMH.asc <- function(model,init_beta=NULL,init_h2,V,famid,prev,data,max.iter=100,max.sub.iter=50,n.cores=1,proband,SAVE=F,out=NULL){
  ## This function depends on a package; tmvtnorm, parallel, mvtnorm
  ## model : Y~X 
  ## init_beta : initial value for beta (p+1 vector)
  ## init_h2 :  initial value for heritability
  ## V : genetic relationship matrix (theoritically, kinship coefficient matrix)
  ## famid : family id
  ## prev : disease prevalence
  ## data : data.frame containing variables specified at formula
  ## proband : proband indicator (1 = proband, 0 = non-proband)
  
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
  
  if(!require(mvtnorm))	{
    install.packages("mvtnorm")
    library(mvtnorm)
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
      return(list(beta=beta.old,h2=h2.old,n_iter=n.iter))
      break
    }
  }
}

getLTMH.S2 <- function(PREV,H2,N.fam,n.cores){
  for(prev in PREV){
    for(h2 in H2){
      for(totalfam in N.fam){
        print(paste0('prevalence:',prev,', heritability:',h2,', number of families:',totalfam))
        library(kinship2)
        setwd("~/paper/heritability/ML_ver2/variousFam/")
        fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
        fin.dat$std_snp <- (fin.dat$snp-mean(fin.dat$snp))/sd(fin.dat$snp)
        fin.dat <- fin.dat[fin.dat$FID%in%fin.dat$FID[fin.dat$ind==1 & fin.dat$Y==1],,drop=F]
        PB <- 'ind'
        
        model <- Y~snp-1
        out <- paste0("~/paper/heritability/ML_ver2/variousFam/prev_",prev,"_h2_",h2,"/Esth2_",totalfam,"_asc_181203.txt")
        if(paste0("Esth2_",totalfam,"_asc_181203.txt") %in% system(paste0("ls prev_",prev,"_h2_",h2),intern=T)){
          start <- as.numeric(strsplit(system(paste0("tail -1 ",out),intern=T)," ")[[1]][1])+1
        } else {
          start <- 1
          write.table(data.frame('obs','beta_std_snp','h2','n_iteration'),out,col.names=F,row.names=F,quote=F)
        }
        Esth2 <- sapply(start:2000,getEsth2.asc,fin.dat=fin.dat,init_beta=NULL,init_h2=h2,totalfam=totalfam,assumed_prev=prev,model=model,n.cores=n.cores,out=out,seed=T,PB=PB)
      }
    }
  }
}

a <- getLTMH.S2(PREV=c(0.1),H2=0.2,N.fam=500,n.cores=32)

