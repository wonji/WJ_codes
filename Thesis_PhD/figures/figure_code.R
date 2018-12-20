
plot_QQ <- function(h2,ha2,prev,method,ylim){

for(type in c("size","power")){
	png(paste("new_QQ_",method,"_",type,"_h2_",h2,"_ha2_",ha2,"_prev_",prev,".png",sep=""), width=1500, height=1000, units = "px", bg = "white", res = 200)
	par(mfrow=c(2,3))

	resultII <- read.table(paste("h2_",h2,"_ha2_",ha2,"_prev_",prev,"/GWAS_",type,"_",method,"_result.txt",sep=""),head=T)

	p.val <- resultII[,c(4,8,12,16,20)]


	i=1
	y <- -log(p.val[,i],10)
	v <- -log10(0.05/sum(is.finite(y)))
	
	# QQ plot
	o.y <- sort(y[is.finite(y)],decreasing=T)
	xx<- (1:length(o.y))/length(o.y)
	x <- -log(xx,10)
	ifelse(max(o.y[is.finite(o.y)])>=x[1],YY<-o.y,YY<-x)

	plot(YY, YY, main='S1',type = "n", xlab = expression(paste("Expected ",-log[10],"(p-value)")), ylab = expression(paste("Observed ",-log[10],"(p-value)")),ylim=c(0,ylim),xlim=c(0,ylim))
	N=length(o.y)
	c95 <- rep(0,N)
	c05 <- rep(0,N)
	for(i in 1:N){
		c95[i] <- qbeta(0.95,i,N-i+1)
		c05[i] <- qbeta(0.05,i,N-i+1)
	}
	#abline(h = c(0:max(max(x),max(o.y[is.finite(o.y)]))), v =c(0:max(max(x),max(o.y[is.finite(o.y)]))), col = "darkgray", lty=3)
	polygon(c(x,sort(x)),c(-log(c95,10),sort(-log(c05,10))),col=c("gray"),border=NA)
	abline(a=0,b=1,col="black", lty=1)
	points(x, o.y, cex = .5, col = "dark red")

	title <- c("S2","S3","S4","S5")
	
	for(i in 2:5){
		y <- -log(p.val[,i],10)
		v <- -log10(0.05/sum(is.finite(y)))
	
		# QQ plot
		o.y <- sort(y[is.finite(y)],decreasing=T)
		xx<- (1:length(o.y))/length(o.y)
		x <- -log(xx,10)
		ifelse(max(o.y[is.finite(o.y)])>=x[1],YY<-o.y,YY<-x)

		plot(YY, YY, main=title[i-1],type = "n", xlab = expression(paste("Expected ",-log[10],"(p-value)")), ylab = expression(paste("Observed ",-log[10],"(p-value)")),ylim=c(0,ylim),xlim=c(0,ylim))
		N=length(o.y)
		c95 <- rep(0,N)
		c05 <- rep(0,N)
		for(i in 1:N){
			c95[i] <- qbeta(0.95,i,N-i+1)
			c05[i] <- qbeta(0.05,i,N-i+1)
		}
		#abline(h = c(0:max(max(x),max(o.y[is.finite(o.y)]))), v =c(0:max(max(x),max(o.y[is.finite(o.y)]))), col = "darkgray", lty=3)
		polygon(c(x,sort(x)),c(-log(c95,10),sort(-log(c05,10))),col=c("gray"),border=NA)
		abline(a=0,b=1,col="black", lty=1)
		points(x, o.y, cex = .5, col = "dark red")
	}
	dev.off()
}
}

# Figure 3.2
setwd('/home2/wjkim/paper/famhis/prob/simulation/nuc_fam')
plot_QQ(0.2,0.005,0.1,"moment",4.5)

# Figure 3.3
setwd('/home2/wjkim/paper/famhis/prob/simulation/ext_fam')
plot_QQ(0.2,0.005,0.1,"moment",4.5)

# Figure 3.4
setwd('/home2/wjkim/paper/famhis/prob/simulation/')
plot_QQ(0.2,0.005,0.1,"moment",4)



# Figure 3.6

#Horizontally - paper#
setwd("~/paper/famhis/prob/realdata")
png(paste("/home2/wjkim/paper/famhis/prob/realdata/new_kare_all_QQ_horizontal.png",sep=""), width=1400, height=500, units = "px", bg = "white", res = 200)
par(mfrow=c(1,3))

path <- c("plotdata_all_indiv_0925.txt","plotdata_case_1000_cont_4000_rs.txt","plotdata_case_1000_cont_4000_rp.txt")
title <- c("All subjects","GWAS using S1","GWAS using S3")

for(i in 1:3){
  print(i)
  
  resultII <- read.table(path[i],head=T)
  
  p.val <- resultII[,3]
  
  y <- -log(p.val,10)
  v <- -log10(0.05/sum(is.finite(y)))
  
  # QQ plot
  o.y <- sort(y[is.finite(y)],decreasing=T)
  xx<- (1:length(o.y))/length(o.y)
  x <- -log(xx,10)
  ifelse(max(o.y[is.finite(o.y)])>=x[1],YY<-o.y,YY<-x)
  
  plot(YY, YY, main=title[i],type = "n", xlab = expression(paste("Expected ",-log[10],"(p-value)")), ylab = expression(paste("Observed ",-log[10],"(p-value)")),ylim=c(0,18),xlim=c(0,18))
  N=length(o.y)
  c95 <- rep(0,N)
  c05 <- rep(0,N)
  for(ii in 1:N){
    c95[ii] <- qbeta(0.95,ii,N-ii+1)
    c05[ii] <- qbeta(0.05,ii,N-ii+1)
  }
  #abline(h = c(0:max(max(x),max(o.y[is.finite(o.y)]))), v =c(0:max(max(x),max(o.y[is.finite(o.y)]))), col = "darkgray", lty=3)
  polygon(c(x,sort(x)),c(-log(c95,10),sort(-log(c05,10))),col=c("gray"),border=NA)
  abline(a=0,b=1,col="black", lty=1)
  points(x, o.y, cex = .5, col = "dark red")
  
}

dev.off()


# Figure 3.7
#Horizontally - paper#
setwd("~/paper/famhis/prob/realdata/permutation")
png(paste("/home2/wjkim/paper/famhis/prob/realdata/permutation/new_kare_all_QQ_horizontal.png",sep=""), width=1400, height=500, units = "px", bg = "white", res = 200)
par(mfrow=c(1,3))

path <- c("plotdata_case_1000_cont_4000_All.txt","plotdata_case_1000_cont_4000_rs.txt","plotdata_case_1000_cont_4000_rp.txt")
title <- c("GWAS using all subjects","GWAS using S1","GWAS using S3")

for(i in 1:3){
  print(i)
  
  resultII <- read.table(path[i],head=T)
  
  p.val <- resultII[,3]
  
  y <- -log(p.val,10)
  v <- -log10(0.05/sum(is.finite(y)))
  
  # QQ plot
  o.y <- sort(y[is.finite(y)],decreasing=T)
  xx<- (1:length(o.y))/length(o.y)
  x <- -log(xx,10)
  ifelse(max(o.y[is.finite(o.y)])>=x[1],YY<-o.y,YY<-x)
  
  plot(YY, YY, main=title[i],type = "n", xlab = expression(paste("Expected ",-log[10],"(p-value)")), ylab = expression(paste("Observed ",-log[10],"(p-value)")))
  N=length(o.y)
  c95 <- rep(0,N)
  c05 <- rep(0,N)
  for(ii in 1:N){
    c95[ii] <- qbeta(0.95,ii,N-ii+1)
    c05[ii] <- qbeta(0.05,ii,N-ii+1)
  }
  #abline(h = c(0:max(max(x),max(o.y[is.finite(o.y)]))), v =c(0:max(max(x),max(o.y[is.finite(o.y)]))), col = "darkgray", lty=3)
  polygon(c(x,sort(x)),c(-log(c95,10),sort(-log(c05,10))),col=c("gray"),border=NA)
  abline(a=0,b=1,col="black", lty=1)
  points(x, o.y, cex = .5, col = "dark red")
  
}

dev.off()




#Horizontally - paper#
setwd("~/paper/famhis/prob/realdata")

path <- c("plotdata_all_indiv_0925.txt","plotdata_case_1000_cont_4000_rs.txt","plotdata_case_1000_cont_4000_rp.txt")
title <- c("All subjects","GWAS using S1","GWAS using S3")

All <- read.table(path[1],head=T,stringsAsFactor=F)
sigs <- All[,2][which(All[,3]<1.872e-7)]

dats <- lapply(2:3,function(i) read.table(path[i],head=T,stringsAsFactor=F))
pvals <- cbind(dats[[1]],dats[[2]][,3])
colnames(pvals)[3:4] <- c('pval_s1','pval_s3')
SIGS <- which(pvals[,2]%in%sigs)

png(paste("/home2/wjkim/paper/famhis/prob/realdata/Pairwise_Pvalue.png",sep=""), width=900, height=1000, units = "px", bg = "white", res = 200)
plot(-log10(pvals[,3]),-log10(pvals[,4]),xlab=expression(paste("S1 ",-log[10],"(p-value)")),ylab=expression(paste("S3 ",-log[10],"(p-value)")),xlim=c(0,ceiling(-log10(min(pvals[,3:4])))),ylim=c(0,ceiling(-log10(min(pvals[,3:4])))))
abline(a=0,b=1,col="black", lty=1)
points(-log10(pvals[SIGS,3]),-log10(pvals[SIGS,4]),col='red')
dev.off()




for(i in 1:3){
  print(i)
  
  resultII <- read.table(path[i],head=T)
  
  p.val <- resultII[,3]
  
  y <- -log(p.val,10)
  v <- -log10(0.05/sum(is.finite(y)))
  
  # QQ plot
  o.y <- sort(y[is.finite(y)],decreasing=T)
  xx<- (1:length(o.y))/length(o.y)
  x <- -log(xx,10)
  ifelse(max(o.y[is.finite(o.y)])>=x[1],YY<-o.y,YY<-x)
  
  plot(YY, YY, main=title[i],type = "n", xlab = expression(paste("Expected ",-log[10],"(p-value)")), ylab = expression(paste("Observed ",-log[10],"(p-value)")),ylim=c(0,18),xlim=c(0,18))
  N=length(o.y)
  c95 <- rep(0,N)
  c05 <- rep(0,N)
  for(ii in 1:N){
    c95[ii] <- qbeta(0.95,ii,N-ii+1)
    c05[ii] <- qbeta(0.05,ii,N-ii+1)
  }
  #abline(h = c(0:max(max(x),max(o.y[is.finite(o.y)]))), v =c(0:max(max(x),max(o.y[is.finite(o.y)]))), col = "darkgray", lty=3)
  polygon(c(x,sort(x)),c(-log(c95,10),sort(-log(c05,10))),col=c("gray"),border=NA)
  abline(a=0,b=1,col="black", lty=1)
  points(x, o.y, cex = .5, col = "dark red")
  
}

dev.off()


