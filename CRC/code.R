##################################################
################ Paper review 1 ##################
##################################################

#### Cox ph model
library(survival)
library(kinship2)

setwd("/home2/wjkim/project/CRC/20180426/1.calCE/0.data")
ped <- read.table("/home2/wjkim/project/CRC/20180426/1.calCE/0.data/CRC_final.ped",head=T,stringsAsFactor=F)
ped$time <- ped$CurrentAge
ped$time[ped$CRC==1&!is.na(ped$OnsetAge)] <- ped$OnsetAge[ped$CRC==1&!is.na(ped$OnsetAge)] 
ped$event <- ifelse(ped$CRC==0,0,ifelse(ped$CRC==1&!is.na(ped$OnsetAge),1,2))
ped$type <- ifelse(ped$CRC==0,"right",ifelse(ped$CRC==1&!is.na(ped$OnsetAge),"counting","left"))

# Variables for Family History
total_ped <- pedigree(id=ped$IID,dadid=ped$FATHER,momid=ped$MOTHER,sex=ped$SEX,famid=ped$PID,missid=0)
V <- kinship(total_ped)*2
proband <- which(ped$PROBAND==1&!is.na(ped$OnsetAge))

get.FH <- function(i){
	kin <- V[i,]; names(kin) <- NULL
	FDR <- which(kin==0.5)
	FH1 <- ifelse(sum(ped$CRC[FDR])>0,1,0)
	FH2 <- ifelse(sum(ped$CRC[FDR])==1,1,ifelse(sum(ped$CRC[FDR])>1,2,0))
	n.FDR <- length(FDR)
	aff.FDR <- sum(ped$CRC[FDR]==1)
	prop.FDR <- aff.FDR/n.FDR

	onset <- ped$OnsetAge[which(kin==0.5&ped$CRC==1)]
	if(all(is.na(onset))){
		if(sum(ped$CRC[FDR])==0){
			return(c(n.FDR,aff.FDR,prop.FDR,FH1,FH2,0))
		} else {
			return(c(n.FDR,aff.FDR,prop.FDR,FH1,FH2,1))
		}
	} else {
		onset <- min(onset[!is.na(onset)])
		FH3 <- ifelse(sum(ped$CRC[FDR])==0,0,ifelse(is.na(onset),1,ifelse(onset>70,1,ifelse(onset>60&onset<=70,2,3))))
		return(c(n.FDR,aff.FDR,prop.FDR,FH1,FH2,FH3))
	}
	
}
FH <- t(sapply(proband,get.FH))
colnames(FH) <- c('n.FDR','aff.FDR','prop.FDR',"FH1","FH2","FH3")

ped.proband <- ped[proband,]
ped.proband <- cbind(ped.proband,FH)
ped.proband$CE_redFDR01 <- (ped.proband$CE_redFDR-min(ped.proband$CE_redFDR))/(max(ped.proband$CE_redFDR)-min(ped.proband$CE_redFDR))

## Cox ph model
# 0. reference model
res.FH0 <- coxph(Surv(time,event) ~ factor(SEX),data=ped.proband)

# 1. Family history of CRC
res.FH1 <- coxph(Surv(time,event) ~ factor(SEX)+factor(FH1),data=ped.proband)
anova(res.FH0,res.FH1)

# 2. No. of affected FDRs
res.FH2 <- coxph(Surv(time,event) ~ factor(SEX)+factor(FH2),data=ped.proband)
anova(res.FH0,res.FH2)

# 3. Age at diagnosis of affected FDR
res.FH3 <- coxph(Surv(time,event) ~ factor(SEX)+factor(FH3),data=ped.proband)
anova(res.FH0,res.FH3)

# 4. Proportion of affected FDR
res.FH4 <- coxph(Surv(time,event) ~ factor(SEX)+prop.FDR,data=ped.proband)
anova(res.FH0,res.FH4)

# 5. CE
res.FH5 <- coxph(Surv(time,event) ~ factor(SEX)+CE_redFDR01,data=ped.proband)
anova(res.FH0,res.FH5)

# 5. categorized CE
get.result3 <- function(ce){
	fh <- rep(0,nrow(FH))
	fh[ped.proband$CE_redFDR01>=ce] <- 1

	tmp.dat <- cbind(ped.proband,fh)
	res.FH <- coxph(Surv(time,event) ~ factor(SEX)+factor(fh),data=tmp.dat)
	return(summary(res.FH)$coef[2,])
}

res.cateCE2 <- do.call(rbind,lapply(seq(0.0001,0.002,by=0.00001),get.result3))
res.cateCE <- rbind(res.cateCE2,res.cateCE1)

setwd('/home2/wjkim/project/CRC/20180426/3.coxph/2.output')
png('categorizedCE.png')
x <- seq(0.0001,0.002,by=0.00001)
y <- -log10(res.cateCE2[,5])
lo <- loess(y~x)
plot(x,y)
lines(predict(lo), col='red', lwd=2)
dev.off()




prop.tab <- table(FH[,3])
props <- sort(as.numeric(names(prop.tab)))[-1]

get.result1 <- function(pp){
	FH4 <- rep(0,nrow(FH))
	FH4[FH[,3]>=pp] <- 1

	tmp.dat <- cbind(ped.proband,FH4)
	# 5. Categorized #4
	res.FH5 <- coxph(Surv(time,event) ~ factor(SEX)+factor(FH4),data=tmp.dat)
	return(cbind(summary(res.FH5)$coef[2,,drop=F],AIC=AIC(res.FH5)))
}

res1 <- t(sapply(props,get.result1))
res1 <- cbind(props,res1)
colnames(res1) <- c('Cutoff','Beta','exp(beta)','SE(beta)','Z-value','P-value','AIC')
setwd('/home2/wjkim/project/CRC/20180426/3.coxph/2.output')
write.csv(res1,'sep1.csv',row.names=F,quote=F)

get.result2 <- function(pp1,pp2){
	FH4 <- rep(0,nrow(FH))
	FH4[FH[,3]>=pp1 & FH[,3]<pp2] <- 1
	FH4[FH[,3]>=pp2] <- 2

	tmp.dat <- cbind(ped.proband,FH4)
	# 5. Categorized #4
	res.FH5 <- coxph(Surv(time,event) ~ factor(SEX)+factor(FH4),data=tmp.dat)
	return(cbind(prop1=pp1,prop2=pp2,summary(res.FH5)$coef[2:3,],Pvalue=anova(res.FH0,res.FH5)$P[2],AIC=AIC(res.FH5)))
}

props1 <- props[6:10]
props2 <- props[16:17]

res2.1 <- do.call(rbind,lapply(props1,get.result2,pp2=props2[1]))
res2.2 <- do.call(rbind,lapply(props1,get.result2,pp2=props2[2]))
res2 <- rbind(res2.1,res2.2)
setwd('/home2/wjkim/project/CRC/20180426/3.coxph/2.output')
write.csv(res2,'sep2.csv',row.names=F,quote=F)




get.result4 <- function(pp,nn){
	FH4 <- rep(0,nrow(FH))
	FH4[FH[,3]<pp & FH[,1]<nn] <- 1
	FH4[FH[,3]>=pp & FH[,1]<nn] <- 2
	FH4[FH[,3]>=pp & FH[,1]>=nn] <- 3

	tmp.dat <- cbind(ped.proband,FH4)
	# 5. Categorized #4
	res.FH5 <- coxph(Surv(time,event) ~ factor(SEX)+factor(FH4),data=tmp.dat)
	return(cbind(n_FDR=nn,prop=pp,summary(res.FH5)$coef[2:4,],Pvalue=anova(res.FH0,res.FH5)$P[2]))
}

res1 <- t(sapply(props,get.result1))
res1 <- cbind(props,res1)
setwd('/home2/wjkim/project/CRC/20180426/3.coxph/')
write.csv(res1,'2.output/sep1.csv',row.names=F,quote=F)
write.csv(ped.proband,'1.input/CRC_coxph.csv',row.names=F,quote=F)

######## coxph CRC ???? ????
setwd('/home2/wjkim/project/CRC/20180426/3.coxph/')
ped <- read.csv('1.input/CRC_coxph.csv',head=T,stringsAsFactor=F)
ped$FH4 <- 0
ped$FH4[ped$prop.FDR>=0.125 & ped$prop.FDR<0.3] <- 1
ped$FH4[ped$prop.FDR>=0.3] <- 2
write.csv(ped,'1.input/CRC_coxph.csv',row.names=F,quote=F)


####### coxph

source("/Users/wonjikim/WJ_codes/CRC/prior_ft.txt")
source("/Users/wonjikim/WJ_codes/CRC/REx_Coxph.R")

CRC <- read.csv('/Users/wonjikim/WJ_codes/CRC/CRC_coxph.csv',head=T)
setwd('/Users/wonjikim/WJ_codes/CRC/coxph')
out="./FH1.html";res.FH1 <- REx_Coxph(CRC, time1='time', event='event', cov_set2=c('SEX','FH1'), tdp_time=NULL, tdp_cov=NULL, tdp_var=NULL, vars=c('SEX','FH1'), confi=0.95, Select=FALSE, direct='forward', survival_plot=FALSE, a1_survival_plot=FALSE, CHazard_plot=FALSE, log_survival_plot=FALSE, expestim=TRUE, CI=TRUE, CI_level=0.95,out="./FH1.html");
out="./FH2.html";res.FH2 <- REx_Coxph(CRC, time1='time', event='event', cov_set2=c('SEX','FH2'), tdp_time=NULL, tdp_cov=NULL, tdp_var=NULL, vars=c('SEX','FH2'), confi=0.95, Select=FALSE, direct='forward', survival_plot=FALSE, a1_survival_plot=FALSE, CHazard_plot=FALSE, log_survival_plot=FALSE, expestim=TRUE, CI=TRUE, CI_level=0.95,out="./FH2.html");
out="./FH3.html";res.FH3 <- REx_Coxph(CRC, time1='time', event='event', cov_set2=c('SEX','FH3'), tdp_time=NULL, tdp_cov=NULL, tdp_var=NULL, vars=c('SEX','FH3'), confi=0.95, Select=FALSE, direct='forward', survival_plot=FALSE, a1_survival_plot=FALSE, CHazard_plot=FALSE, log_survival_plot=FALSE, expestim=TRUE, CI=TRUE, CI_level=0.95,out="./FH3.html");
out="./FH4.html";res.FH4 <- REx_Coxph(CRC, time1='time', event='event', cov_set2=c('SEX','FH4'), tdp_time=NULL, tdp_cov=NULL, tdp_var=NULL, vars=c('SEX','FH4'), confi=0.95, Select=FALSE, direct='forward', survival_plot=FALSE, a1_survival_plot=FALSE, CHazard_plot=FALSE, log_survival_plot=FALSE, expestim=TRUE, CI=TRUE, CI_level=0.95,out="./FH4.html");

new.CRC <- CRC[CRC$time<81,]
out="./FH1_under80.html";res.FH1 <- REx_Coxph(new.CRC, time1='time', event='event', cov_set2=c('SEX','FH1'), tdp_time=NULL, tdp_cov=NULL, tdp_var=NULL, vars=c('SEX','FH1'), confi=0.95, Select=FALSE, direct='forward', survival_plot=FALSE, a1_survival_plot=FALSE, CHazard_plot=FALSE, log_survival_plot=FALSE, expestim=TRUE, CI=TRUE, CI_level=0.95,out="./FH1_under80.html");
out="./FH2_under80.html";res.FH2 <- REx_Coxph(new.CRC, time1='time', event='event', cov_set2=c('SEX','FH2'), tdp_time=NULL, tdp_cov=NULL, tdp_var=NULL, vars=c('SEX','FH2'), confi=0.95, Select=FALSE, direct='forward', survival_plot=FALSE, a1_survival_plot=FALSE, CHazard_plot=FALSE, log_survival_plot=FALSE, expestim=TRUE, CI=TRUE, CI_level=0.95,out="./FH2_under80.html");
out="./FH3_under80.html";res.FH3 <- REx_Coxph(new.CRC, time1='time', event='event', cov_set2=c('SEX','FH3'), tdp_time=NULL, tdp_cov=NULL, tdp_var=NULL, vars=c('SEX','FH3'), confi=0.95, Select=FALSE, direct='forward', survival_plot=FALSE, a1_survival_plot=FALSE, CHazard_plot=FALSE, log_survival_plot=FALSE, expestim=TRUE, CI=TRUE, CI_level=0.95,out="./FH3_under80.html");
out="./FH4_under80.html";res.FH4 <- REx_Coxph(new.CRC, time1='time', event='event', cov_set2=c('SEX','FH4'), tdp_time=NULL, tdp_cov=NULL, tdp_var=NULL, vars=c('SEX','FH4'), confi=0.95, Select=FALSE, direct='forward', survival_plot=FALSE, a1_survival_plot=FALSE, CHazard_plot=FALSE, log_survival_plot=FALSE, expestim=TRUE, CI=TRUE, CI_level=0.95,out="./FH4_under80.html");

res.FH0 <- coxph(Surv(time,event) ~ factor(SEX),data=new.CRC)
res.FH2 <- coxph(Surv(time,event) ~ factor(SEX)+factor(FH2),data=new.CRC)





#### Tables

## table 2
crc$famsize <- 'medium'
crc$famsize[crc$n.FDR>=9] <- 'large'
crc$famsize[crc$n.FDR<4] <- 'small'

tab1 <- table(crc$famsize,crc$FH2)
tab2 <- do.call(rbind,lapply(c('large','medium','small'),function(i) data.frame(Mean=mean(crc$prop.FDR[crc$famsize==i]),SD=sd(crc$prop.FDR[crc$famsize==i]),row.names=i)))
tab <- cbind(as.data.frame.matrix(tab1),tab2)

## table 3
crc$DummyOnsetAge <- 0
crc$DummyOnsetAge[crc$OnsetAge <=50] <- 1

tab1 <- table(crc$DummyOnsetAge,crc$FH2)
tab2 <- do.call(rbind,lapply(0:1,function(i) data.frame(Mean=mean(crc$prop.FDR[crc$DummyOnsetAge==i]),SD=sd(crc$prop.FDR[crc$DummyOnsetAge==i]),row.names=i)))
tab <- cbind(as.data.frame.matrix(tab1),tab2)

## table 3-1 (???߿? ǥ??)
tab2 <- do.call(rbind,lapply(0:2,function(i) data.frame(Mean=mean(crc$OnsetAge[crc$FH4==i]),SD=sd(crc$OnsetAge[crc$FH4==i]),row.names=i)))
tab <- cbind(as.data.frame.matrix(tab1),tab2)

## table 5
get.result1 <- function(pp){
	temp.FH <- rep(0,nrow(new.CRC))
	if(pp==0){
		temp.FH[new.CRC[,'prop.FDR']>pp] <- 1
	} else {
		temp.FH[new.CRC[,'prop.FDR']>=pp] <- 1
	}

	tmp.dat <- cbind(new.CRC,temp.FH)
	# 5. Categorized #4
	res.FH5 <- coxph(Surv(time,event) ~ factor(SEX)+factor(temp.FH),data=tmp.dat)
	return(summary(res.FH5)$coef[2,])
}

res1 <- t(sapply(seq(0,0.3,0.025),get.result1))
res1 <- cbind(seq(0,0.3,0.025),res1)

## table 7
get.result2 <- function(pp1,pp2){
	temp.FH <- rep(0,nrow(new.CRC))
	temp.FH[new.CRC[,'prop.FDR']>=pp1 & new.CRC[,'prop.FDR']<pp2] <- 1
	temp.FH[new.CRC[,'prop.FDR']>=pp2] <- 2

	tmp.dat <- cbind(new.CRC,temp.FH)
	# 5. Categorized #4
	res.FH5 <- coxph(Surv(time,event) ~ factor(SEX)+factor(temp.FH),data=tmp.dat)
	return(cbind(prop1=pp1,prop2=pp2,summary(res.FH5)$coef[2:3,],Pvalue=anova(res.FH0,res.FH5)$P[2]))
}

res2 <- do.call(rbind,lapply(seq(0.2,0.325,0.025),get.result2,pp1=0.125))
setwd('/home2/wjkim/project/CRC/20180426/3.coxph/2.output')
write.csv(res2,'sep2.csv',row.names=F,quote=F)



#### Figures

hist(new.CRC$prop.FDR[new.CRC$prop.FDR!=0])


