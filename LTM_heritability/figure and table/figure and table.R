################ Figures and tables

######## Tables ########

##### Table 1 (Descriptive statistics for scenario 1) - Need to be updated
get.res.S1 <- function(prev,h2){
	setwd("~/paper/heritability/ML_ver2/variousFam/")
	print(c(prev,h2))
	LTMH <- read.table(paste0("prev_",prev,"_h2_",h2,"/Esth2_500_181030.txt"),head=T)
	if("Esth2_500_181114.txt" %in% system(paste0("ls prev_",prev,"_h2_",h2),intern=T)){
		tmp <- read.table(paste0("prev_",prev,"_h2_",h2,"/Esth2_500_181114.txt"),head=T)
		LTMH <- rbind(LTMH,tmp)
	}

	GCTA <- read.table(paste0("otherSW/1.GCTA/prev_",prev,"_h2_",h2,"/Esth2_500.txt"),head=T)
	tmp <- read.table(paste0("otherSW/1.GCTA/prev_",prev,"_h2_",h2,"/Esth2_500_181114.txt"),head=T)
	GCTA <- rbind(GCTA,tmp)

	fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
	LTMH$beta <- LTMH$beta_std_snp*sd(fin.dat$snp)		
	# LTMH$beta <- LTMH$beta_std_snp/sd(fin.dat$snp)		
		
	res.beta <- data.frame(beta=paste0(round(mean(LTMH$beta),4)," (",round(sd(LTMH$beta),4),")"))
	res.h2.LTMH <- data.frame(h2_LTMH=paste0(round(mean(LTMH$h2),4)," (",round(sd(LTMH$h2),4),")"))
	res.h2.GCTA <- data.frame(h2_GCTA=paste0(round(mean(GCTA$h2),4)," (",round(sd(GCTA$h2),4),")"))
	res <- cbind(res.beta,res.h2.LTMH,res.h2.GCTA)
	rownames(res) <- paste0("prev_",prev,"_h2_",h2)
	return(res)
}

fin.res <- lapply(c(0.05,0.2,0.4),function(hh) do.call(rbind,lapply(c(0.05,0.1,0.2),get.res.S1,h2=hh)))
fin.res.I <- do.call(rbind,fin.res)
write.csv(fin.res.I,"figure_and_table/table1.csv",quote=F)


##### Table 2 (Descriptive statistics for scenario 2)
get.res.S2 <- function(prev,h2){
	setwd("~/paper/heritability/ML_ver2/variousFam/")
	
	if(prev==0.1 & h2==0.05){
		LTMH <- read.table(paste0("prev_",prev,"_h2_",h2,"/Esth2_500_asc_181120.txt"),head=T)
	} else {
		LTMH <- read.table(paste0("prev_",prev,"_h2_",h2,"/Esth2_500_asc.txt"),head=T)
		if("Esth2_500_newasc.txt" %in% system(paste0("ls prev_",prev,"_h2_",h2),intern=T)){
			tmp <- read.table(paste0("prev_",prev,"_h2_",h2,"/Esth2_500_newasc.txt"),head=T)
			LTMH <- rbind(LTMH,tmp)
		}
		if("Esth2_500_asc_181114.txt" %in% system(paste0("ls prev_",prev,"_h2_",h2),intern=T)){
			tmp <- read.table(paste0("prev_",prev,"_h2_",h2,"/Esth2_500_asc_181114.txt"),head=T)
			LTMH <- rbind(LTMH,tmp)
		}
	}

	GCTA <- read.table(paste0("otherSW/1.GCTA/prev_",prev,"_h2_",h2,"/Esth2_500_asc.txt"),head=T)
	tmp <- read.table(paste0("otherSW/1.GCTA/prev_",prev,"_h2_",h2,"/Esth2_500_asc_181114.txt"),head=T)
	GCTA <- rbind(GCTA,tmp)

	fin.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/dataset.txt"),head=T,stringsAsFactor=F)
	LTMH$beta <- LTMH$beta_std_snp*sd(fin.dat$snp)		
		
	res.beta <- data.frame(beta=paste0(round(mean(LTMH$beta),4)," (",round(sd(LTMH$beta),4),")"))
	res.h2.LTMH <- data.frame(h2_LTMH=paste0(round(mean(LTMH$h2),4)," (",round(sd(LTMH$h2),4),")"))
	res.h2.GCTA <- data.frame(h2_GCTA=paste0(mean(GCTA$h2)," (",sd(GCTA$h2),")"))
	res <- cbind(res.beta,res.h2.LTMH,res.h2.GCTA)
	rownames(res) <- paste0("prev_",prev,"_h2_",h2)
	return(res)
}

fin.res <- lapply(c(0.05,0.2,0.4),function(hh) do.call(rbind,lapply(c(0.05,0.1,0.2),get.res.S2,h2=hh)))
fin.res.I <- do.call(rbind,fin.res)
write.csv(fin.res.I,"figure_and_table/table2.csv",quote=F)


#### Table 3 (Empirical size & power table for h2)
setwd("~/paper/heritability/ML_ver2/variousFam/CEST/1.h2")

CEST.h2 <- function(prev,h2){
	setwd("~/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
	LTMH.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/CEST_h2_500_all.txt"),head=T,stringsAsFactor=F)
	tmp.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/CEST_h2_500_all_181116.txt"),head=T,stringsAsFactor=F)
	LTMH.dat <- cbind(LTMH.dat,tmp.dat)
	LTMH.res <- sapply(c(0.01,0.05,0.1),function(i) sum(LTMH.dat$Pvalue<i)/nrow(LTMH.dat))
	
	return(LTMH.res)
}

fin.res <- lapply(c(0,0.2,0.4),function(hh) do.call(rbind,lapply(c(0.05,0.1,0.2),CEST.h2,h2=hh)))
fin.res.I <- do.call(rbind,fin.res)
colnames(fin.res.I) <- paste('siglevel_',c(0.01,0.05,0.1))
fin.res.I <- cbind(h2=rep(c(0,0.2,0.4),each=3),prev=rep(c(0.05,0.1,0.2),3),fin.res.I)

setwd("~/paper/heritability/ML_ver2/variousFam/figure_and_table")
write.csv(fin.res.I,"table3_h2_CEST_S1.csv",quote=F,row.names=F)



#### Table 4-1 (Empirical size table for beta)
setwd("~/paper/heritability/ML_ver2/variousFam/CEST/2.beta")

CEST.beta.size <- function(prev,h2){
	setwd("~/paper/heritability/ML_ver2/variousFam/CEST/2.beta")
	dat <- read.table(paste0("prev_",prev,"_h2_",h2,"_size/CEST_beta_500_tau.txt"),head=T,stringsAsFactor=F)
	res <- matrix(sapply(c(0.01,0.05,0.1),function(i) sum(dat$Pvalue<i)/nrow(dat)),ncol=1)
	return(res)
}

fin.res <- lapply(c(0.2,0.4),function(hh) do.call(cbind,lapply(c(0.1,0.2),CEST.beta.size,h2=hh)))
fin.res.I <- do.call(rbind,fin.res)

colnames(fin.res.I) <- paste0("Prev_",c(0.1,0.2))
fin.res.I <- cbind(h2=rep(c(0.2,0.4),each=3),sig_level=rep(c(0.01,0.05,0.1),2),fin.res.I)
setwd("~/paper/heritability/ML_ver2/variousFam/figure_and_table")
write.csv(fin.res.I,"table5_CEST_beta_size.csv",quote=F)


#### Table 4-2 (Empirical power table for beta)
setwd("~/paper/heritability/ML_ver2/variousFam/CEST/2.beta")

CEST.beta.power <- function(prev,h2){
	setwd("~/paper/heritability/ML_ver2/variousFam/CEST/2.beta")
	dat <- read.table(paste0("prev_",prev,"_h2_",h2,"_power/CEST_beta_500_tau.txt"),head=T,stringsAsFactor=F)
	res <- matrix(sapply(c(0.01,0.05,0.1),function(i) sum(dat$Pvalue<i)/nrow(dat)),ncol=1)
	return(res)
}

fin.res <- lapply(c(0.2,0.4),function(hh) do.call(cbind,lapply(c(0.1,0.2),CEST.beta.power,h2=hh)))
fin.res.I <- do.call(rbind,fin.res)

colnames(fin.res.I) <- paste0("Prev_",c(0.1,0.2))
fin.res.I <- cbind(h2=rep(c(0.2,0.4),each=3),sig_level=rep(c(0.01,0.05,0.1),2),fin.res.I)
setwd("~/paper/heritability/ML_ver2/variousFam/figure_and_table")
write.csv(fin.res.I,"table6_CEST_beta_power.csv",quote=F)


#### Table 5 Empirical size and power for scenario 2
setwd("~/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
CEST.h2 <- function(prev,h2){
	setwd("~/paper/heritability/ML_ver2/variousFam/CEST/1.h2")
	LTMH.dat <- read.table(paste0("prev_",prev,"_h2_",h2,"/CEST_h2_500_asc_181120.txt"),head=T,stringsAsFactor=F)
	LTMH.res <- sapply(c(0.01,0.05,0.1),function(i) sum(LTMH.dat$Pvalue<i)/nrow(LTMH.dat))
	
	return(LTMH.res)
}

fin.res <- lapply(c(0,0.2,0.4),function(hh) do.call(rbind,lapply(c(0.05,0.1,0.2),CEST.h2,h2=hh)))
fin.res.I <- do.call(rbind,fin.res)
colnames(fin.res.I) <- paste('siglevel_',c(0.01,0.05,0.1))
fin.res.I <- cbind(h2=rep(c(0,0.2,0.4),each=3),prev=rep(c(0.05,0.1,0.2),3),fin.res.I)

setwd("~/paper/heritability/ML_ver2/variousFam/figure_and_table")
write.csv(fin.res.I,"table5_h2_CEST_S2.csv",quote=F,row.names=F)




#### Table 6 (Demographic information for T2D)
fam <- read.table("/data/kare/allmerge/allmerge_withrela.fam",head=F,stringsAsFactor=F)
new.fam <- fam[-grep("^KNIH",fam$V1),]
pheno <- read.csv("/home2/wjkim/project/kare+snu/snu/snu_pheno.csv",head=T,stringsAsFactor=F)
rela <- read.table("/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/1.SNUH/SNUH_rela_age.txt",head=T,stringsAsFactor=F)

case <- cbind(pheno[pheno$GWAS.ID!="",c('GWAS.ID','AGE')],5)
rela <- rela[,c('IID','AGE','RELATIONSHIP')]
colnames(case) <- c('IID','AGE','RELATIONSHIP')
age.info <- rbind(case,rela)

mer <- merge(new.fam,age.info,by.x='V2',by.y='IID',sort=F,all=F)
ind <- which(mer$V6!=-9 & !is.na(mer$AGE))
dataset <- mer[ind,,drop=F]
PB.candi <- grep('PKS',dataset$V2)
PB.fam <- dataset$V1[PB.candi]
PB <- PB.candi[!duplicated(PB.fam)]
proband <- rep(0,nrow(dataset))
proband[PB] <- 1

T2D <- table(dataset$V6,proband)
T2D.perc <- apply(T2D,2,function(jj) jj/sum(jj)*100) 
res.T2D <- paste0(paste0(T2D," (",round(T2D.perc,2),"%)"),collapse=" ")

sex <- table(dataset$V5,proband)
sex.perc <- apply(sex,2,function(jj) jj/sum(jj)*100) 
res.sex <- paste0(paste0(sex," (",round(sex.perc,2),"%)"),collapse=" ")

age <- sapply(c(0,1),function(jj) mean(dataset$AGE[proband==jj]))
sd.age <- sapply(c(0,1),function(jj) sd(dataset$AGE[proband==jj]))
res.age <- paste(paste0(round(age,2)," (",round(sd.age,2),")"),collapse=" ")

get.demo <- cbind(NPB=c('T2D','Sex','Age'),PB=c(res.T2D,res.sex,res.age))

rela.dat <- dataset[proband==0,]
rela.dat[rela.dat$RELATIONSHIP==5,'RELATIONSHIP'] <- c(7,6,4,7,4,4,4,7,7,4,7,7,6,4,4,4,7,4,7,4,6,7,6,7,6,7,4,4,4)
RT <- as.matrix(table(rela.dat$RELATIONSHIP))
RT <- cbind(c('Sibling','Parents','Offspring'),RT)

setwd("~/paper/heritability/ML_ver2/variousFam/figure_and_table")
write.table(get.demo,"table7_1.csv",quote=F,row.names=F)
write.table(RT,"table7_2.csv",quote=F,row.names=F,col.names=F)


#### Table 8 (Familial relationship)
fam <- read.table("/data/kare/allmerge/allmerge_withrela.fam",head=F,stringsAsFactor=F)
new.fam <- fam[-grep("^KNIH",fam$V1),]
pheno <- read.csv("/home2/wjkim/project/kare+snu/snu/snu_pheno.csv",head=T,stringsAsFactor=F)
rela <- read.table("/home2/wjkim/paper/heritability/ML_ver2/variousFam/realdata/1.SNUH/SNUH_rela_age.txt",head=T,stringsAsFactor=F)

case <- cbind(pheno[pheno$GWAS.ID!="",c('GWAS.ID','AGE')],5)
rela <- rela[,c('IID','AGE','RELATIONSHIP')]
colnames(case) <- c('IID','AGE','RELATIONSHIP')
age.info <- rbind(case,rela)

mer <- merge(new.fam,age.info,by.x='V2',by.y='IID',sort=F,all=F)
ind <- which(mer$V6!=-9 & !is.na(mer$AGE))
dataset <- mer[ind,,drop=F]
dataset[,7] <- sample(c('A','C'),nrow(dataset),replace=T)
dataset[,8] <- sample(c('A','C'),nrow(dataset),replace=T)
dataset <- dataset[,c(2,1,3:8)]
setwd("~/paper/heritability/ML_ver2/variousFam/figure_and_table")
write.table(dataset,"temp_SNUH.ped",row.names=F,col.names=F,quote=F)
write.table(data.frame(1,'SNP1',0,1),"temp_SNUH.map",row.names=F,col.names=F,quote=F)

system("plink --file temp_SNUH --recode vcf --out temp_SNUH")
system("plink --file temp_SNUH --make-bed --out temp_SNUH")
a <- read.table("temp_SNUH.fam",head=F,stringsAsFactor=F)
a[!a[,3]%in%a[,2],3] <- 0
a[!a[,4]%in%a[,2],4] <- 0
write.table(a,"temp_SNUH.fam",row.names=F,col.names=F,quote=F)
system("onetool --fam temp_SNUH.fam --vcf temp_SNUH.vcf --out SNUH_pedinfo")



######## Figures ########
## Fig 1 : KKT conditions
setwd("~/paper/heritability/ML_ver2/variousFam/")
png("figure_and_table/fig1_KKT.png",width=1900,height=650,res=200)
par(mfrow=c(1,3))

f1 <- function(x) -5*(x+0.5)^2+5
f1_p <- function(x) -10*(x+0.5)
jup <- function(x,y) f1_p(x)*(y-x)+f1(x)

xdat <- seq(-2,2,0.01)
ydat <- sapply(xdat, f1)

plot(xdat,ydat,type='n',xlim=c(-1,2),ylim=c(-10,10),xlab=expression(italic(h)^2),ylab=expression(paste("Q(",italic(h)^2,")")),main=expression(paste(argmax[italic(h)^2],' Q(',italic(h)^2,') < 0')))
polygon(c(0,1,1,0),c(-10.8,-10.8,10.78,10.78),col='light gray',border=NA)
lines(xdat,ydat)
abline(v=c(0,1),lty=2)
lines(c(-0.5,0.5),c(jup(0,-0.5),jup(0,0.5)),col='red',lwd=2)
lines(c(0.7,1.3),c(jup(1,0.7),jup(1,1.3)),col='red',lwd=2)
points(c(-0.5,0,1),c(5,f1(0),f1(1)),pch=c(2,19,19))
mtext("A",SNORTH<-3,line=2,adj=0.01,cex=1)


f1 <- function(x) -5*(x-1.5)^2+5
f1_p <- function(x) -10*(x-1.5)
jup <- function(x,y) f1_p(x)*(y-x)+f1(x)

xdat <- seq(-2,2,0.01)
ydat <- sapply(xdat, f1)

plot(xdat,ydat,type='n',xlim=c(-1,2),ylim=c(-10,10),xlab=expression(italic(h)^2),ylab=expression(paste("Q(",italic(h)^2,")")),main=expression(paste(argmax[italic(h)^2],' Q(',italic(h)^2,') > 1')))
polygon(c(0,1,1,0),c(-10.8,-10.8,10.78,10.78),col='light gray',border=NA)
lines(xdat,ydat)
abline(v=c(0,1),lty=2)
lines(c(-0.7,0.3),c(jup(0,-0.7),jup(0,0.3)),col='red',lwd=2)
lines(c(0.7,1.3),c(jup(1,0.7),jup(1,1.3)),col='red',lwd=2)
points(c(1.5,0,1),c(5,f1(0),f1(1)),pch=c(2,19,19))
mtext("B",SNORTH<-3,line=2,adj=0.01,cex=1)

f1 <- function(x) -5*(x-0.5)^2+5
f1_p <- function(x) -10*(x-0.5)
jup <- function(x,y) f1_p(x)*(y-x)+f1(x)

xdat <- seq(-2,3,0.01)
ydat <- sapply(xdat, f1)

plot(xdat,ydat,type='n',xlim=c(-1,2),ylim=c(-10,10),xlab=expression(italic(h)^2),ylab=expression(paste("Q(",italic(h)^2,")")),main=expression(paste(argmax[italic(h)^2],' Q(',italic(h)^2,') in [0,1]')))
polygon(c(0,1,1,0),c(-10.8,-10.8,10.78,10.78),col='light gray',border=NA)
lines(xdat,ydat)
abline(v=c(0,1),lty=2)
lines(c(-0.3,0.3),c(jup(0,-0.3),jup(0,0.3)),col='red',lwd=2)
lines(c(0.7,1.3),c(jup(1,0.7),jup(1,1.3)),col='red',lwd=2)
points(c(0.5,0,1),c(5,f1(0),f1(1)),pch=c(2,19,19))
mtext("C",SNORTH<-3,line=2,adj=0.01,cex=1)
dev.off()

#### Figure 2 - color
library(ggplot2)
library(gridExtra)
setwd("~/paper/heritability/ML_ver2/variousFam/")
png("figure_and_table/fig2_boxplot.png",width=1200,height=1200,res=200)
get.plot <- function(h2){
	setwd("~/paper/heritability/ML_ver2/variousFam/")
	get.dat <- function(h2,prev){
		LTMH <- read.table(paste0("prev_",prev,"_h2_",h2,"/Esth2_500_181030.txt"),head=T)
		if("Esth2_500_181114.txt" %in% system(paste0("ls prev_",prev,"_h2_",h2),intern=T)){
			tmp <- read.table(paste0("prev_",prev,"_h2_",h2,"/Esth2_500_181114.txt"),head=T)
			LTMH <- rbind(LTMH,tmp)
		}

		GCTA <- read.table(paste0("otherSW/1.GCTA/prev_",prev,"_h2_",h2,"/Esth2_500.txt"),head=T)
		tmp <- read.table(paste0("otherSW/1.GCTA/prev_",prev,"_h2_",h2,"/Esth2_500_181114.txt"),head=T)
		GCTA <- rbind(GCTA,tmp)
		res <- data.frame(h2=c(LTMH$h2,GCTA$h2),Prevalence=rep(prev,nrow(LTMH)+nrow(GCTA)),Method=c(rep("LTMH",nrow(LTMH)),rep("GCTA",nrow(GCTA))))
		res$Method <- factor(res$Method,level=c("LTMH","GCTA"))
		return(res)
	}
	
	dat <- do.call(rbind,lapply(c(0.05,0.1,0.2),get.dat,h2=h2))
	dat$Prevalence <- factor(dat$Prevalence)
	p <- ggplot(dat,aes(x=Prevalence,y=h2,fill=Method))+geom_boxplot()+theme_classic()+geom_hline(yintercept=h2,linetype='dashed',color='gray')+annotate("text",x=factor(0.2),y=(max(dat$h2)-0.05),label=paste0("italic(h)^2==",h2),parse=T)
	if(h2==0.4){
		p <- p+labs(y=expression(hat(italic(h))^2))
	} else {
		p <- p+labs(y=expression(hat(italic(h))^2),x="")
	}	
	return(p)
}
myPlot <- lapply(c(0.05,0.2,0.4),get.plot)
grid.arrange(grobs=myPlot,mrow=3)
dev.off()


#### Figure 2 - black & white
library(ggplot2)
library(gridExtra)
setwd("~/paper/heritability/ML_ver2/variousFam/")
png("figure_and_table/fig2_boxplot.png",width=1200,height=1200,res=200)
get.plot <- function(h2){
	setwd("~/paper/heritability/ML_ver2/variousFam/")
	get.dat <- function(h2,prev){
		LTMH <- read.table(paste0("prev_",prev,"_h2_",h2,"/Esth2_500_181030.txt"),head=T)
		if("Esth2_500_181114.txt" %in% system(paste0("ls prev_",prev,"_h2_",h2),intern=T)){
			tmp <- read.table(paste0("prev_",prev,"_h2_",h2,"/Esth2_500_181114.txt"),head=T)
			LTMH <- rbind(LTMH,tmp)
		}

		GCTA <- read.table(paste0("otherSW/1.GCTA/prev_",prev,"_h2_",h2,"/Esth2_500.txt"),head=T)
		tmp <- read.table(paste0("otherSW/1.GCTA/prev_",prev,"_h2_",h2,"/Esth2_500_181114.txt"),head=T)
		GCTA <- rbind(GCTA,tmp)
		res <- data.frame(h2=c(LTMH$h2,GCTA$h2),Prevalence=rep(prev,nrow(LTMH)+nrow(GCTA)),Method=c(rep("LTMH",nrow(LTMH)),rep("GCTA",nrow(GCTA))))
		res$Method <- factor(res$Method,level=c("LTMH","GCTA"))
		return(res)
	}
	
	dat <- do.call(rbind,lapply(c(0.05,0.1,0.2),get.dat,h2=h2))
	dat$Prevalence <- factor(dat$Prevalence)
	p <- ggplot(dat,aes(x=Prevalence,y=h2,fill=Method))+geom_boxplot()+theme_classic()+geom_hline(yintercept=h2,linetype='dashed',color='gray')+annotate("text",x=factor(0.2),y=(max(dat$h2)-0.05),label=paste0("italic(h)^2==",h2),parse=T)
	if(h2==0.4){
		p <- p+labs(y=expression(hat(italic(h))^2))
	} else {
		p <- p+labs(y=expression(hat(italic(h))^2),x="")
	}	
	return(p)
}
myPlot <- lapply(c(0.05,0.2,0.4),get.plot)
grid.arrange(grobs=myPlot,mrow=3)
dev.off()


#### Figure 3
V <- matrix(1/2,4,4)
diag(V) <- 1
V[1,2] <- V[2,1] <- 0

h2=0.2944
beta.hat=matrix(0.8,1,1)
thres <- qnorm(0.109,lower.tail=F)

get.risk <- function(Y,age,V,h2,beta.hat,thres,ind){
	library(mvtnorm)
	X <- matrix((c(age+29,age+26,age,age-3)-48.62664)/15.7037,ncol=1)

	# Nominator
	a <- ifelse(Y==0,-Inf,thres)
	b <- ifelse(Y==0,thres,Inf)
	a[ind] <- thres; b[ind] <- Inf
	S <- h2*V+(1-h2)*diag(length(Y))
	nom <- pmvnorm(lower=a,upper=b,mean=as.vector(X%*%beta.hat),sigma=S)[1]
	
	# Denominator
	a <- a[-ind]; b <- b[-ind]
	S <- S[-ind,-ind]
	X <- X[-ind,]
	denom <- pmvnorm(lower=a,upper=b,mean=as.vector(X%*%beta.hat),sigma=S)[1]
	
	res <- data.frame(age=age,prob=nom/denom)
	return(res)
}

get.risk.uni <- function(age,beta.hat,thres){
	library(mvtnorm)
	X <- (age-48.62664)/15.7037
	res.1 <- pnorm(thres,mean=X*beta.hat,sd=1,lower.tail=F)
	res <- data.frame(age=age,prob=res.1)
	return(res)
}
	
library(parallel)
# 1. Y=c(0,0,0,0)
Y=c(0,0,0,0)
Res <- mclapply(seq(20,80,by=1),get.risk,Y=Y,V=V,h2=h2,beta.hat=beta.hat,thres=thres,ind=3,mc.cores=24)
Res.I <- do.call(rbind,Res)
Res.I <- Res.I[order(Res.I$age),]
	
# 2. Y=c(1,0,0,0)
Y=c(1,0,0,0)
Res <- mclapply(seq(20,80,by=1),get.risk,Y=Y,V=V,h2=h2,beta.hat=beta.hat,thres=thres,ind=3,mc.cores=24)
Res.II <- do.call(rbind,Res)
Res.II <- Res.II[order(Res.II$age),]

# 3. Y=c(1,1,0,0)
Y=c(1,1,0,0)
Res <- mclapply(seq(20,80,by=1),get.risk,Y=Y,V=V,h2=h2,beta.hat=beta.hat,thres=thres,ind=3,mc.cores=24)
Res.III <- do.call(rbind,Res)
Res.III <- Res.III[order(Res.III$age),]

# 4. Y=c(1,1,0,1)
Y=c(1,1,0,1)
Res <- mclapply(seq(20,80,by=1),get.risk,Y=Y,V=V,h2=h2,beta.hat=beta.hat,thres=thres,ind=3,mc.cores=24)
Res.IV <- do.call(rbind,Res)
Res.IV <- Res.IV[order(Res.IV$age),]

# 5. random sample
Res <- lapply(seq(20,80,by=1),get.risk.uni,beta.hat=beta.hat,thres=thres)
Res.V <- do.call(rbind,Res)

# Fig 3A
setwd("~/paper/heritability/ML_ver2/variousFam/figure_and_table/")
png("fig3_Riskscore.png",width=1200,height=580,res=100)
par(mfrow=c(1,2))
par(oma=c(0, 0, 0, 5))
plot(1,1,type='n',xlim=c(20,80),ylim=c(Res.I[1,'prob'],Res.IV[61,'prob']),xlab="Age of proband",ylab="Probability of being affected")
lines(seq(20,80,by=1),Res.I$prob,lty=1)
lines(seq(20,80,by=1),Res.II$prob,lty=2)
lines(seq(20,80,by=1),Res.III$prob,lty=3)
lines(seq(20,80,by=1),Res.IV$prob,lty=4)
mtext(text="(A)", side=3, at=10, cex=1.2, line=1)

# Fig 3B
Res.I$RR <- Res.I$prob/Res.V$prob
Res.II$RR <- Res.II$prob/Res.V$prob
Res.III$RR <- Res.III$prob/Res.V$prob
Res.IV$RR <- Res.IV$prob/Res.V$prob

plot(1,1,type='n',xlim=c(20,80),ylim=c(0,max(Res.IV[,'RR'])),xlab="Age of proband",ylab="Relative Risk")
lines(seq(20,80,by=1),Res.I$RR,lty=1)
lines(seq(20,80,by=1),Res.II$RR,lty=2)
lines(seq(20,80,by=1),Res.III$RR,lty=3)
lines(seq(20,80,by=1),Res.IV$RR,lty=4)
legend("topright", inset=c(-0.3,0.07), legend=paste0(0:3), lty=1:4, title="Number of\naffected\nrelatives",xpd=NA,bty="n")
mtext(text="(B)", side=3, at=10, cex=1.2, line=1)
dev.off()


#### Figure 4
Res.I$OR <- (Res.I$prob/(1-Res.I$prob))/(Res.V$prob/(1-Res.V$prob))
Res.II$OR <- (Res.II$prob/(1-Res.II$prob))/(Res.V$prob/(1-Res.V$prob))
Res.III$OR <- (Res.III$prob/(1-Res.III$prob))/(Res.V$prob/(1-Res.V$prob))
Res.IV$OR <- (Res.IV$prob/(1-Res.IV$prob))/(Res.V$prob/(1-Res.V$prob))

setwd("~/paper/heritability/ML_ver2/variousFam/figure_and_table/")
png("fig4_OR.png",width=600,height=580,res=100)
par(oma=c(0, 0, 0, 5))
plot(1,1,type='n',xlim=c(20,80),ylim=c(0,max(Res.IV[,'OR'])),xlab="Age of proband",ylab="Odds Ratio")
lines(seq(20,80,by=1),Res.I$OR,lty=1)
lines(seq(20,80,by=1),Res.II$OR,lty=2)
lines(seq(20,80,by=1),Res.III$OR,lty=3)
lines(seq(20,80,by=1),Res.IV$OR,lty=4)
legend("topright", inset=c(-0.3,0.07), legend=paste0(0:3), lty=1:4, title="Number of\naffected\nrelatives",xpd=NA,bty="n")
dev.off()


#### Figure 5
setwd('/home2/wjkim/paper/heritability/ML_ver2/variousFam/CEST/3.realdata/1.LAM/2.Whole_CHR/')
a <- sapply(1:22,function(chr) return(ifelse('GRM_181024.txt'%in%system(paste0('ls chr',chr),intern=T),chr,'')))
aa <- a[a!='']
dat <- c()
for(i in aa){
	bb <- read.table(paste0('chr',i,'/GRM_181024.txt'),head=T)
	dat <- rbind(dat,bb)
	print(paste0('chr',i,', ',nrow(bb),' SNPs'))
}
plotdat <- dat[,c('CHR','SNP','Pvalue')]
colnames(plotdat) <- c('CHR','GENE','p.val')
write.table(plotdat,'LAM_plotdata.txt',row.names=F,quote=F)
system('Rscript ~/Output_plot.r LAM_plotdata.txt LAM_GWAS')

library("GenABEL")
lambda <- estlambda(plotdat$p.val,method="median")	#1.075946

library(gap)
p.val.gc <- pchisq(gcontrol2(plotdat$p.val)$y/lambda$estimate,df=1,lower.tail=F)
lambda.gc <- estlambda(p.val.gc,method="median")
plotdat.gc <- plotdat
plotdat.gc$p.val <- p.val.gc
write.table(plotdat.gc,'LAM_plotdata_gc.txt',row.names=F,quote=F)
setwd('/home2/wjkim/paper/heritability/ML_ver2/variousFam/')
system('Rscript ~/Output_plot.r CEST/3.realdata/1.LAM/2.Whole_CHR/LAM_plotdata_gc.txt figure_and_table/LAM_GWAS_gc')


