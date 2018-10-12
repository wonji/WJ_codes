################ Revision #####################

############# New CLR excluding 2 cases based on the MDS plots ##################

#### Draw MDS plots excluding 2 pairs
#setwd("C:\\Users\\rewki\\Documents\\projects\\LAM\\Final_analysis\\figure and table")
setwd("/home2/wjkim/project/LAM/final_analysis/19.MDSplot/")
pc <- read.table("1.PCs/PCs.pca.res",head=T)
fam <- read.table("0.data/Disc_1000GP.fam",head=F)
obs <- 1:nrow(pc)
pc <- cbind(pc,obs)
pc <- pc[!((pc$PC1< -4| pc$PC2 >2) & fam$V6==2),]

png("MDS_Disc_1000GP_wo_Outlier.png", width=1800, height=650, units = "px", bg = "white", res = 200)
phe <- fam[,6]
col <- ifelse(phe==1,"dark blue",ifelse(phe==2,"dark red","light gray"))
par(mfrow=c(1,3))
# PC1 vs PC2
with(pc,plot(PC1,PC2,col=col,main="MDS plot : PC1 and PC2",type='n'))
points(pc$PC1[phe==-9],pc$PC2[phe==-9],col=col[phe==-9])
points(pc$PC1[phe==1],pc$PC2[phe==1],col=col[phe==1])
points(pc$PC1[phe==2],pc$PC2[phe==2],col=col[phe==2])

# PC1 vs PC3
with(pc,plot(PC1,PC3,col=col,main="MDS plot : PC1 and PC3",type='n'))
points(pc$PC1[phe==-9],pc$PC3[phe==-9],col=col[phe==-9])
points(pc$PC1[phe==1],pc$PC3[phe==1],col=col[phe==1])
points(pc$PC1[phe==2],pc$PC3[phe==2],col=col[phe==2])

# PC2 vs PC3
with(pc,plot(PC2,PC3,col=col,main="MDS plot : PC2 and PC3",type='n'))
points(pc$PC2[phe==-9],pc$PC3[phe==-9],col=col[phe==-9])
points(pc$PC2[phe==1],pc$PC3[phe==1],col=col[phe==1])
points(pc$PC2[phe==2],pc$PC3[phe==2],col=col[phe==2])
dev.off()



########### Method 1 :  Do matching newly #############
#### Get list of two subjects
setwd("/home2/wjkim/project/LAM/final_analysis/19.MDSplot/")
pc <- read.table("1.PCs/PCs.pca.res",head=T)
fam <- read.table("0.data/Disc_1000GP.fam",head=F)

rm.list <- pc[(pc$PC1< -4| pc$PC2 >2) & fam$V6==2,1:2]
setwd("/home2/wjkim/project/LAM/final_analysis/3.clogit/1.revision/0.data")
write.table(rm.list,"outliers.txt",col.names=F,row.names=F,quote=F) # two cases


#### Do CLR
##Matcing##
library(Matching)
setwd("/home2/wjkim/project/LAM/final_analysis/")
Phen <- read.table("0.data/phenofile_withage.txt",head=T,stringsAsFactor=F)
rm.list <- read.table("3.clogit/1.revision/0.data/outliers.txt",head=F,stringsAsFactor=F)
phen <- Phen[!Phen$IID%in%rm.list$V2,]

res <- Match(Tr=phen$pheno01, X=phen[,c("PC1","PC2","Age")],replace=F,M=2)
cont <- res$index.control; case <- res$index.treat
matching.index <- cbind(case,cont)

dat1 <- matching.index[seq(1,nrow(matching.index),2),]
dat2 <- matching.index[seq(2,nrow(matching.index),2),]

matching.index <- temp <- merge(dat1,dat2,by="case",sort=F)
matching.index[,1] <- phen[temp[,1],'IID']
matching.index[,2] <- phen[temp[,2],'IID']
matching.index[,3] <- phen[temp[,3],'IID']

colnames(matching.index) <- c("case","control1","control2")
write.table(matching.index,"3.clogit/1.revision/0.data/matching_index_1:2.txt",row.names=F,quote=F)


#make pair variable#
phe <- read.table("0.data/phenofile_withage.txt",head=T,stringsAsFactor=F)
match <- read.table("3.clogit/1.revision/0.data/matching_index_1:2.txt",head=T,stringsAsFactor=F)
pair <- rep(NA,nrow(phe))
for(i in 1:nrow(match)) pair[which(phe$IID%in%match[i,])] <- i
phe$pair_revision <- pair
write.table(phe,"/home2/wjkim/project/LAM/final_analysis/3.clogit/1.revision/0.data/phenofile_withage.txt",row.names=F,quote=F)



##conditional logit##
library(survival)
library(parallel)

setwd("/home2/wjkim/project/LAM/final_analysis/0.data")
bim <- read.table("final_clean_withage.bim",head=F,stringsAsFactor=F)

for(i in 1:22){
  print(paste("chromosome ",i,sep=""))
  setwd("/home2/wjkim/project/LAM/final_analysis/3.clogit")
  raw <- read.table(paste("chr",i,"/lam_chr",i,".raw",sep=""),head=T)
  phe <- read.table("1.revision/0.data/phenofile_withage.txt",head=T)
  n.snp <- ncol(raw)-6
  name.snp <- colnames(raw)[-c(1:6)]
  name.snp <- gsub("_[ATGC]","",name.snp)
  
  pair.data <- phe[which(!is.na(phe$pair_revision)),]
  pair.raw <- raw[which(!is.na(phe$pair_revision)),]
  
  con.logit <- function(ii){
    if(ii%%100==0) print(ii)
    snp <- pair.raw[,6+ii]
    dat <- cbind(pair.data,snp)
    fit.lam <- clogit(pheno01 ~ snp + strata(pair), method="exact", data=dat)
    res.coef <- summary(fit.lam)$coef
    res <- c(i,name.snp[ii],res.coef)

    return(res)
  }
  outI <- mclapply(1:n.snp,con.logit,mc.cores=24)
  outII <- do.call(rbind,outI)
  colnames(outII) <- c("chr","snp","beta","exp(beta)","se(beta)","z","p_value")
  setwd("/home2/wjkim/project/LAM/final_analysis/3.clogit/1.revision")
  write.table(outII,paste("chr",i,"/clogit_chr",i,"_181012.txt",sep=""),row.names=F,quote=F)
}

clogit <- c()
for(i in 1:22){
  print(i)
  setwd("/home2/wjkim/project/LAM/final_analysis/3.clogit/1.revision")
  res <- read.table(paste("chr",i,"/clogit_chr",i,"_181012.txt",sep=""),head=T)
  dat <- res[,c(1,2,7)]
  clogit <- rbind(clogit,dat)
}
colnames(clogit) <- c("CHR","GENE","p.val")
setwd("/home2/wjkim/project/LAM/final_analysis/3.clogit/1.revision")
write.table(clogit,"plotdata_clogit_181012.txt",row.names=F,quote=F,sep="\t")

system("Rscript ~/Output_plot.r plotdata_clogit_181012.txt clogit_181012")



################# Method 2 - Use previous matching result ######################
##conditional logit##
library(survival)
library(parallel)

setwd("/home2/wjkim/project/LAM/final_analysis/")
bim <- read.table("0.data/final_clean_withage.bim",head=F,stringsAsFactor=F)
phe <- read.table("0.data/phenofile_withage.txt",head=T)
rm.list <- read.table("3.clogit/1.revision/0.data/outliers.txt",head=F,stringsAsFactor=F)
phe$pair_age[which(phe$pair_age%in%phe[phe$IID%in%rm.list$V2,'pair_age'])] <- NA


for(i in 1:22){
  print(paste("chromosome ",i,sep=""))
  setwd("/home2/wjkim/project/LAM/final_analysis/3.clogit")
  raw <- read.table(paste("chr",i,"/lam_chr",i,".raw",sep=""),head=T)
  n.snp <- ncol(raw)-6
  name.snp <- colnames(raw)[-c(1:6)]
  name.snp <- gsub("_[ATGC]","",name.snp)
  
  pair.data <- phe[which(!is.na(phe$pair_age)),]
  pair.raw <- raw[which(!is.na(phe$pair_age)),]
  
  con.logit <- function(ii){
    if(ii%%100==0) print(ii)
    snp <- pair.raw[,6+ii]
    dat <- cbind(pair.data,snp)
    fit.lam <- clogit(pheno01 ~ snp + strata(pair), method="exact", data=dat)
    res.coef <- summary(fit.lam)$coef
    res <- c(i,name.snp[ii],res.coef)
    
    return(res)
  }
  outI <- mclapply(1:n.snp,con.logit,mc.cores=32)
  outII <- do.call(rbind,outI)
  colnames(outII) <- c("chr","snp","beta","exp(beta)","se(beta)","z","p_value")
  setwd("/home2/wjkim/project/LAM/final_analysis/3.clogit/1.revision")
  write.table(outII,paste("chr",i,"/clogit_chr",i,"_181012_method2.txt",sep=""),row.names=F,quote=F)
}

clogit <- c()
for(i in 1:22){
  print(i)
  setwd("/home2/wjkim/project/LAM/final_analysis/3.clogit/1.revision")
  res <- read.table(paste("chr",i,"/clogit_chr",i,"_181012_method2.txt",sep=""),head=T)
  dat <- res[,c(1,2,7)]
  clogit <- rbind(clogit,dat)
}
colnames(clogit) <- c("CHR","GENE","p.val")
setwd("/home2/wjkim/project/LAM/final_analysis/3.clogit/1.revision")
write.table(clogit,"plotdata_clogit_181012_method2.txt",row.names=F,quote=F,sep="\t")

system("Rscript ~/Output_plot.r plotdata_clogit_181012_method2.txt clogit_181012_method2")
