####### Table 2
#### Get a risk according to the age and sex
get.risk <- function(age,sex,mean.age,mean.sex,sd.age,sd.sex,thresh,out){
  aa <- read.table(out,head=T)
  std.age <- (age - mean.age)/sd.age
  std.sex <- (sex - mean.sex)/sd.sex
  new.thresh <- thresh - std.age*aa[nrow(aa),1] - std.sex*aa[nrow(aa),2]
  risk <- pnorm(new.thresh,lower.tail=F)
  return(data.frame(Age=age,Sex=sex,Risk=risk))
}
setwd("/home2/wjkim/project/CRC/20171020/1.model_based/2.mendel/0.data")
out <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/2.CRC/1.FDR/LTMH_FDR_agesex.txt"
fin.dat <- read.table("reduced_FDR_mendel.in",head=F,sep=",",stringsAsFactor=F)
fin.dat[fin.dat[,3]=='',3] <- 0
fin.dat[fin.dat[,4]=='',4] <- 0
colnames(fin.dat)[c(1:5,14,17)] <- c('FID','IID','PID','MID','SEX','DummyAge1','CRC')
fin.dat$SEX <- ifelse(fin.dat$SEX=='M',0,1)
prev=0.002326
thresh <- qnorm(prev,lower.tail=F)
M_under50 <- get.risk(age=0,sex=0,mean.age=mean(fin.dat$DummyAge1),mean.sex=mean(fin.dat$SEX),sd.age=sd(fin.dat$DummyAge1),sd.sex=sd(fin.dat$SEX),thresh=thresh,out=out)
FM_under50 <- get.risk(age=0,sex=1,mean.age=mean(fin.dat$DummyAge1),mean.sex=mean(fin.dat$SEX),sd.age=sd(fin.dat$DummyAge1),sd.sex=sd(fin.dat$SEX),thresh=thresh,out=out)
M_over50 <- get.risk(age=1,sex=0,mean.age=mean(fin.dat$DummyAge1),mean.sex=mean(fin.dat$SEX),sd.age=sd(fin.dat$DummyAge1),sd.sex=sd(fin.dat$SEX),thresh=thresh,out=out)
FM_over50 <- get.risk(age=1,sex=1,mean.age=mean(fin.dat$DummyAge1),mean.sex=mean(fin.dat$SEX),sd.age=sd(fin.dat$DummyAge1),sd.sex=sd(fin.dat$SEX),thresh=thresh,out=out)

a <- rbind(FM_under50,M_under50,FM_over50,M_over50)
a$RR <- a$Risk/a$Risk[1]
write.table(a,"~/paper/CRC/1.tables/Coef_effect_risk.txt",quote=F)


###### Figure 1
#### Get a risk according to the age and sex
get.famrisk <- function(Y,X,thresh,out){
  library(mvtnorm)
  aa <- read.table(out,head=T)
  betas <- t(as.matrix(aa[nrow(aa),1:2]))
  a <- ifelse(Y==0,-Inf,thresh)
  b <- ifelse(Y==0,thresh,Inf)
  a[3] <- thresh; b[3] <- Inf
  V <- matrix(0.5,4,4)
  V[1,2] <- V[2,1] <- 0
  S <- aa[nrow(aa),3]*V+(1-aa[nrow(aa),3])*diag(4)
  
  nom <- pmvnorm(lower=a,upper=b,mean=as.vector(X%*%betas),sigma=S)[1]
  denom <- pmvnorm(lower=a[-3],upper=b[-3],mean=as.vector(X%*%betas)[-3],sigma=S[-3,-3])[1]
  
  risk <- nom/denom
  return(risk)
}

kk <- function(j,thresh,out){
  ages <- list(c(0,0,0,0),c(1,1,0,0),c(1,1,1,1))  # 1) all under 50, 2) Parents over, 50 3) all over 50
  std.age <- (ages[[j]] - mean(fin.dat$DummyAge1))/sd(fin.dat$DummyAge1)
  std.sex <- (c(0,1,0,0) - mean(fin.dat$SEX))/sd(fin.dat$SEX)
  X <- cbind(std.age,std.sex)
  Y <- list(c(0,0,0,0),c(1,0,0,0),c(1,1,0,0),c(1,1,0,1))
  gg <- function(i) get.famrisk(Y[[i]],X,thresh,out)
  ggg <- sapply(1:4,gg)
  
  aa <- read.table(out,head=T)
  betas <- t(as.matrix(aa[nrow(aa),1:2]))
  baseline.risk <- pnorm(thresh,mean=X[3,,drop=F]%*%betas,sd=1,lower.tail=F)
  return(c(ggg,baseline.risk))
}

setwd("/home2/wjkim/project/CRC/20171020/1.model_based/2.mendel/0.data")
out <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/2.CRC/1.FDR/LTMH_FDR_agesex.txt"
fin.dat <- read.table("reduced_FDR_mendel.in",head=F,sep=",",stringsAsFactor=F)
fin.dat[fin.dat[,3]=='',3] <- 0
fin.dat[fin.dat[,4]=='',4] <- 0
colnames(fin.dat)[c(1:5,14,17)] <- c('FID','IID','PID','MID','SEX','DummyAge1','CRC')
fin.dat$SEX <- ifelse(fin.dat$SEX=='M',0,1)
prev=0.002326
thresh <- qnorm(prev,lower.tail=F)

age1 <- kk(1,thresh,out)
age2 <- kk(2,thresh,out)
age3 <- kk(3,thresh,out)

res <- RR <- rbind(age1,age2,age3)
colnames(res) <- colnames(RR) <- c('pheno1','pheno2','pheno3','pheno4','baseline.risk')
for(i in 1:3) RR[i,] <- RR[i,]/RR[i,5]

write.table(res,"~/paper/CRC/2.figures/Familial_risk.txt",quote=F)
write.table(RR,"~/paper/CRC/2.figures/Familial_RR.txt",quote=F)

options(stringsAsFactors=F)
res <- read.table("~/paper/CRC/2.figures/Familial_risk.txt",head=T)
png("~/paper/CRC/2.figures/CRC_familial_risk.png", width=1200, height=1200, units = "px", bg = "white", res = 200)
plot(1,1,type='n',xlim=c(0,3),ylim=c(0,max(res)),xlab="Number of affected FDRs",ylab="Risk probability",xaxt='n')
lines(x=0:3,y=res[1,1:4],lty=1)
lines(x=0:3,y=res[2,1:4],lty=2)
axis(side=1,at=c(0,1,2,3))
legend("topleft",lty=c(1,2),legend=c('Under 50','Over 50'),title='Age of parents',xpd=NA,bty="n")
dev.off()

options(stringsAsFactors=F)
res <- read.table("~/paper/CRC/2.figures/Familial_RR.txt",head=T)
png("~/paper/CRC/2.figures/CRC_familial_risk.png", width=1200, height=1200, units = "px", bg = "white", res = 200)
plot(1,1,type='n',xlim=c(0,2),ylim=c(0,max(res)),xlab="Number of affected FDRs",ylab="Risk probability",xaxt='n')
lines(x=0:3,y=res[1,1:4],lty=1)
lines(x=0:3,y=res[2,1:4],lty=2)
axis(side=1,at=c(0,1,2,3))
legend("topleft",lty=c(1,2),legend=c('Under 50','Over 50'),title='Age of parents',xpd=NA,bty="n")
dev.off()



#############Figure 2
######### Proportion of affected relatives
get.risk <- function(n.fdr,n.aff,fin.dat,out,thresh){
  aa <- read.table(out,head=T)
  betas <- t(as.matrix(aa[nrow(aa),1:2]))

  Y <- rep(0,n.fdr+1)
  Y[1:(1+n.aff)] <- 1
  std.age <- (rep(1,n.fdr+1) - mean(fin.dat$DummyAge1))/sd(fin.dat$DummyAge1)
  std.sex <- (rep(0,n.fdr+1) - mean(fin.dat$SEX))/sd(fin.dat$SEX)
  X <- cbind(std.age,std.sex)
  V <- matrix(0.5,n.fdr+1,n.fdr+1)
  S <- aa[nrow(aa),3]*V+(1-aa[nrow(aa),3])*diag(n.fdr+1)
  
  a <- ifelse(Y==0,-Inf,thresh)
  b <- ifelse(Y==0,thresh,Inf)
  nom <- pmvnorm(lower=a,upper=b,mean=as.vector(X%*%betas),sigma=S)[1]
  denom <- pmvnorm(lower=a[-1],upper=b[-1],mean=as.vector(X%*%betas)[-1],sigma=S[-1,-1])[1]
  risk <- nom/denom

  return(data.frame(n.FDR=n.fdr,n.aff=n.aff,Risk=risk))
}

setwd("/home2/wjkim/project/CRC/20171020/1.model_based/2.mendel/0.data")
out <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/2.CRC/1.FDR/LTMH_FDR_agesex.txt"
fin.dat <- read.table("reduced_FDR_mendel.in",head=F,sep=",",stringsAsFactor=F)
fin.dat[fin.dat[,3]=='',3] <- 0
fin.dat[fin.dat[,4]=='',4] <- 0
colnames(fin.dat)[c(1:5,14,17)] <- c('FID','IID','PID','MID','SEX','DummyAge1','CRC')
fin.dat$SEX <- ifelse(fin.dat$SEX=='M',0,1)
prev=0.002326
thresh <- qnorm(prev,lower.tail=F)

res.risk.2 <- do.call(rbind,lapply(2:6,get.risk,n.aff=2,fin.dat=fin.dat,out=out,thresh=thresh))
res.risk.2$RR <- res.risk.2$Risk/res.risk.2$Risk[1]
write.table(res.risk.2,"~/paper/CRC/1.tables/Risk_nFDR.txt",quote=F)

res <- read.table("~/paper/CRC/1.tables/Risk_nFDR.txt",head=T)
png("~/paper/CRC/2.figures/CRC_risk_nFDR.png", width=1200, height=1200, units = "px", bg = "white", res = 200)
plot(1,1,type='n',xlim=c(0,4),ylim=c(min(res[,4]),max(res[,4])),xlab="Number of unaffected FDRs",ylab="Relative Risk",xaxt='n')
lines(x=0:4,y=res[,4],lty=1)
points(x=0:4,y=res[,4],pch=20)
axis(side=1,at=c(0:4))
dev.off()



############ Figure 4
setwd('/home2/wjkim/project/CRC/20180426/3.coxph/2.output')
sep1 <- read.csv("sep1.csv",head=T)
png("~/paper/CRC/2.figures/Coxph_AIC.png", width=2400, height=1200, units = "px", bg = "white", res = 200)
par(mfrow=c(1,2))
plot(sep1$Cutoff,sep1$AIC,type='b',xlab="First cutoff",ylab="AIC",pch=1)
points(sep1$Cutoff[sep1$AIC==min(sep1$AIC)],sep1$AIC[sep1$AIC==min(sep1$AIC)],pch=16)
text(sep1$Cutoff[sep1$AIC==min(sep1$AIC)],sep1$AIC[sep1$AIC==min(sep1$AIC)],pos=4,labels="0.125")
abline(v=0.125,lty=2,col="gray")
mtext("a)",SNORTH<-3,line=2,adj=-0.15,cex=1.5)

sep1 <- read.csv("sep2.csv",head=T)
plot(sep1$prop2,sep1$AIC,type='b',xlab="Second cutoff",ylab="AIC",pch=1)
points(sep1$prop2[sep1$AIC==min(sep1$AIC)],sep1$AIC[sep1$AIC==min(sep1$AIC)],pch=16)
#text(sep1$prop2[sep1$AIC==min(sep1$AIC)],sep1$AIC[sep1$AIC==min(sep1$AIC)],pos=4,labels="0.125")
#abline(v=0.125,lty=2,col="gray")
mtext("b)",SNORTH<-3,line=2,adj=-0.15,cex=1.5)
dev.off()














get.famrisk <- function(Y,X,thresh,out){
  library(mvtnorm)
  aa <- read.table(out,head=T)
  betas <- t(as.matrix(aa[nrow(aa),1:2]))
  a <- ifelse(Y==0,-Inf,thresh)
  b <- ifelse(Y==0,thresh,Inf)
  a[3] <- thresh; b[3] <- Inf
  V <- matrix(0.5,4,4)
  V[1,2] <- V[2,1] <- 0
  S <- aa[nrow(aa),3]*V+(1-aa[nrow(aa),3])*diag(4)
  
  nom <- pmvnorm(lower=a,upper=b,mean=as.vector(X%*%betas),sigma=S)[1]
  denom <- pmvnorm(lower=a[-3],upper=b[-3],mean=as.vector(X%*%betas)[-3],sigma=S[-3,-3])[1]
  
  risk <- nom/denom
  return(risk)
}

kk <- function(j,prev,out){
  thresh <- qnorm(prev,lower.tail=F)
  ages <- list(c(0,0,0,0),c(1,1,0,0),c(1,1,1,1))  # 1) all under 50, 2) Parents over, 50 3) all over 50
  std.age <- (ages[[j]] - mean(fin.dat$DummyAge1))/sd(fin.dat$DummyAge1)
  std.sex <- (c(0,1,0,0) - mean(fin.dat$SEX))/sd(fin.dat$SEX)
  X <- cbind(std.age,std.sex)
  Y <- list(c(0,0,0,0),c(1,0,0,0),c(1,1,0,0),c(1,1,0,1))
  gg <- function(i) get.famrisk(Y[[i]],X,thresh,out)
  ggg <- sapply(1:4,gg)
  
  aa <- read.table(out,head=T)
  betas <- t(as.matrix(aa[nrow(aa),1:2]))
  baseline.risk <- pnorm(thresh,mean=X[3,,drop=F]%*%betas,sd=1,lower.tail=F)
  res<- c(ggg,baseline.risk)
  res.RR <- ggg/baseline.risk
  return(c(res,res.RR))
}

setwd("/home2/wjkim/project/CRC/20171020/1.model_based/2.mendel/0.data")
out <- "/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/2.CRC/1.FDR/LTMH_FDR_agesex.txt"
fin.dat <- read.table("reduced_FDR_mendel.in",head=F,sep=",",stringsAsFactor=F)
fin.dat[fin.dat[,3]=='',3] <- 0
fin.dat[fin.dat[,4]=='',4] <- 0
colnames(fin.dat)[c(1:5,14,17)] <- c('FID','IID','PID','MID','SEX','DummyAge1','CRC')
fin.dat$SEX <- ifelse(fin.dat$SEX=='M',0,1)

prevs <- c(0.001,0.002326,0.005,0.01,0.05,0.1)
g <- lapply(prevs,kk,j=3,out=out)
gg <- do.call(rbind,g)
colnames(gg) <- c(paste0('risk_n.aff_',0:3),'baseline.risk',paste0('RR_n.aff_',0:3))
fin.res <- data.frame(Prevalence=prevs,gg)
write.table(fin.res,"~/paper/CRC/2.figures/Familial_risk_RR_prevs.txt",row.names=F,quote=F)

options(stringsAsFactors=F)
res <- read.table("~/paper/CRC/2.figures/Familial_risk_RR_prevs.txt",head=T)
png("~/paper/CRC/2.figures/CRC_familial_RR_prevs.png", width=1200, height=1200, units = "px", bg = "white", res = 200)
plot(1,1,type='n',xlim=c(0,3),ylim=c(min(res[,7]),max(res[,10])),xlab="Number of affected FDRs",ylab="Relative risk",xaxt='n')
for(i in 1:6) lines(x=0:3,y=res[i,7:10],lty=i,col=i)
axis(side=1,at=c(0,1,2,3))
prevs <- c(0.001,0.002326,0.005,0.01,0.05,0.1)
legend("topleft",lty=c(1:6),col=1:6,legend=prevs,title='Prevalence',xpd=NA,bty="n")
dev.off()


