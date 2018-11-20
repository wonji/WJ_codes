#독립표본T검정

REx_twoTtest <- function(dataset, res_var, group_var, alternative = "two.sided", CI=F, conf.level = 0.95, levenes.test=T, var.equal=F, digit=4, mu=0)  {

	## Variables 
	# res_var		: 변수선택 탭의 반응변수.
	# group_var		: 변수선택 탭의 집단변수.
	# mu			: 분석옵션 탭의 모평균의 차이.
	# alternative		: 분석옵션 탭의 검정방법. 'two-sided' or 'less' or 'greater'
	# levenes.test		: 분석옵션 탭의 등분산 검정. a levene's test to test equalness of variance b/w groups
	# val.equal		: 분석옵션 탭의 등분산 가정. t-test with equal or unequal variance assumption
	# CI			: 분석옵션 탭의 신뢰구간 출력
	# conf.level		: 분석옵션 탭의 유의수준
	# digit 		: # of digits to be printed

	# Required packages : R2HTML, car
	load.pkg(c("R2HTML", "car"))
	options(digits=digit)

	html.output <- capture.output({
		# HTML file title
		R2HTML::HTML(R2HTML::as.title("Two Sample T-test"), HR = 1, file = stdout(), append = F)
		dataset <- dataset[complete.cases(dataset[,c(res_var,group_var),drop=F]),]
	
		# warnings
		if(!is.numeric(dataset[,res_var])) warn.msg1 <-  "\a Error : Response variable should be numeric. Analysis has been stopped."
		if(length(unique(dataset[,group_var]))!=2) warn.msg2 <- paste0("\a Error : Dataset should consist of two groups. (Observed groups: ",paste(as.character(unique(dataset[,group_var])),collapse=", "),") Analysis has been stopped.")
		if(is.numeric(dataset[,group_var])) dataset[,group_var] <- factor(dataset[,group_var])
		if(exists('warn.msg1')|exists('warn.msg2')){
			R2HTML::HTML(R2HTML::as.title("Warnings"), HR = 2, file = stdout())
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1, file=stdout(), innerBorder = 1, align = "left") 
			if(exists('warn.msg2')) R2HTML::HTML(warn.msg2, file=stdout(), innerBorder = 1, align = "left") 
		} else {

			# Data structure
			group1 <- as.character(unique(dataset[,group_var])[1]); n.group1 <- sum(dataset[,group_var]==group1)
			group2 <- as.character(unique(dataset[,group_var])[2]); n.group2 <- sum(dataset[,group_var]==group2)

			data_str <- matrix(c("Number of observations", paste("Number of observations in group ",group1,sep=""),paste("Number of observations in group ",group2,sep=""), nrow(dataset), n.group1,n.group2), ncol=2)
			R2HTML::HTML(R2HTML::as.title("Data Structure"), HR = 2, file = stdout())
			R2HTML::HTML(data_str, file=stdout(), row.names = FALSE, col.names = TRUE, innerBorder = 1, align = "left") ;
			
			# Variable list
			R2HTML::HTML(R2HTML::as.title("Variable List"), HR = 2, file = stdout())
			VL <- matrix(c('Response variable','Group variable',res_var,group_var), ncol=2)
			R2HTML::HTML(VL, file = stdout(), row.names = FALSE, col.names = TRUE, innerBorder = 1, align = "left")
			
			# Analysis Description
			R2HTML::HTML(R2HTML::as.title("Analysis Description"), HR = 2, file = stdout())
			t.test.res <- t.test(formula(paste(res_var,'~',group_var)),data=dataset, alternative = alternative, var.equal=var.equal, conf.level = conf.level, mu=mu)
			H1 <- capture.output(t.test.res)[grep("alternative",capture.output(t.test.res))]
			H1 <- gsub(".+: ","",H1)
			AD <- matrix(c('Alternative hypothesis',H1,'Test for homogeneity variance',levenes.test,'Assumption of homogeneity variance',var.equal),ncol=2,byrow=T)
			if(CI) AD <- rbind(AD,c('Significance level for CI',conf.level))
			R2HTML::HTML(AD, file = stdout(), row.names = FALSE, col.names = TRUE, innerBorder = 1, align = "left")

			# Descriptive Statistics
			R2HTML::HTML(R2HTML::as.title("Descriptive Statistics"), HR = 2, file = stdout())
			DS <- t(sapply(c(group1,group2),function(i) return(c(length(dataset[dataset[,group_var]==i,res_var]),mean(dataset[dataset[,group_var]==i,res_var]),sd(dataset[dataset[,group_var]==i,res_var])))))
			DS <- as.data.frame(DS)
			colnames(DS) <- c('# of obs.','Mean','SD')
			rownames(DS) <- paste('Group',rownames(DS))
			R2HTML::HTML(DS, file = stdout(), row.names = T, col.names = TRUE, innerBorder = 1, align = "left",digits=6)

			# Levene's Test for Homogeneity of Variance between Groups
			if (levenes.test) {	
				R2HTML::HTML(R2HTML::as.title("Levene's Test for Homogeneity of Variance between Groups"), HR = 2, file = stdout())
				res.lev <- leveneTest(y=dataset[,res_var],group=dataset[,group_var])
				res.lev.F <- res.lev$F[1]
				res.lev.DF <- c(res.lev[1])[[1]]
				res.lev.P <- res.lev$"Pr(>F)"[1]
				LT <- data.frame(res.lev.F,res.lev.DF[1],res.lev.DF[2],res.lev.P)
				colnames(LT) <- c('F-value','DF1','DF2','P-value')
				rownames(LT) <- res_var
				R2HTML::HTML(LT, file = stdout(),col.names = TRUE, innerBorder = 1, align = "left")
			}

			# Test results
			R2HTML::HTML(R2HTML::as.title("Results of Two Sample T-Test"), HR = 2, file =stdout())
			t.test.res.I <- as.data.frame(matrix(with(t.test.res,c(estimate[1],estimate[2],estimate[1]-estimate[2],statistic,parameter,p.value)),nrow=1))
			colnames(t.test.res.I) <- c(names(t.test.res$estimate),'Difference of the means','T-value','DF','P-value')
			rownames(t.test.res.I) <- res_var
			if(CI) {
				t.test.res.I <- cbind(t.test.res.I,t.test.res$conf.int[1],t.test.res$conf.int[2])
				data.frame(t.test.res.I,t.test.res$conf.int[1],t.test.res$conf.int[2])
				colnames(t.test.res.I)[7:8] <- 	c(paste0('Lower bound of ',conf.level*100,'% CI'),paste0('Upper bound of ',conf.level*100,'% CI'))
			}
			R2HTML::HTML(t.test.res.I, file = stdout(), col.names = TRUE, innerBorder = 1, align = "left",digits=5)
		}

		#### End of output
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Two Sample T-test",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})
	
	return(html.output)
}
