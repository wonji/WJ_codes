# 짝진표본T검정

REx_pairedTtest <- function(dataset, res_var1, res_var2, alternative = "two.sided", CI=F, conf.level = 0.95, levenes.test=T, var.equal=F, digit=4, mu=0)  {
	## Variables 
	# res_var1		: 변수선택 탭의 첫번째 반응변수.
	# res_var2		: 변수선택 탭의 두번째 반응변수.
	# mu			: 분석옵션 탭의 모평균의 차이.
	# alternative		: 분석옵션 탭의 검정방법. 'two-sided' or 'less' or 'greater'
	# CI			: 분석옵션 탭의 신뢰구간 출력
	# conf.level		: 분석옵션 탭의 유의수준
	# digit 		: # of digits to be printed

	# Required packages : R2HTML
	load.pkg("R2HTML")
	options(digits=digit)

	html.output <- capture.output({
		# HTML file title
		R2HTML::HTML(R2HTML::as.title("Paired Sample T-test"), HR = 1, file = stdout(), append = F)

		n.sample <- nrow(dataset)
		n.miss.var1 <- sum(!complete.cases(dataset[,res_var1,drop=F]))
		n.miss.var2 <- sum(!complete.cases(dataset[,res_var2,drop=F]))
		dataset <- dataset[complete.cases(dataset[,c(res_var1,res_var2),drop=F]),]
	
		# warnings
		if(!is.numeric(dataset[,res_var1])|!is.numeric(dataset[,res_var2])) {
			R2HTML::HTML(R2HTML::as.title("Warnings"), HR = 2, file = stdout())
			warn.msg1 <-  "\a Error : Response variable should be numeric. Analysis has been stopped."
			R2HTML::HTML(warn.msg1, file=stdout(), innerBorder = 1, align = "left") 
		} else {
			# Data structure
			#data_str <- matrix(c("Number of observations", nrow(dataset)), ncol=2)
			data_str <- matrix(c("Total number of paired observations", n.sample, paste0("Number of missing observations of ",res_var1),n.miss.var1,paste0("Number of missing observations of ",res_var2),n.miss.var2), nrow=3,byrow=T)

			R2HTML::HTML(R2HTML::as.title("Data Structure"), HR = 2, file = stdout())
			R2HTML::HTML(data_str, file=stdout(), row.names = FALSE, col.names = TRUE, innerBorder = 1, align = "left") ;
			
			# Analysis Description
			R2HTML::HTML(R2HTML::as.title("Analysis Description"), HR = 2, file = stdout())
			t.test.res <- t.test(x=dataset[,res_var1],y=dataset[,res_var2],paired=T,alternative = alternative,mu=mu)
			H1 <- capture.output(t.test.res)[grep("alternative",capture.output(t.test.res))]
			H1 <- gsub(".+: ","",H1)
			AD <- matrix(c('Response variables',paste(c(res_var1,res_var2),collapse=', '),'Alternative hypothesis',H1),ncol=2,byrow=T)
			if(CI) AD <- rbind(AD,c('Significance level for CI',conf.level))
			R2HTML::HTML(AD, file = stdout(), row.names = FALSE, col.names = TRUE, innerBorder = 1, align = "left")

			# Descriptive Statistics
			R2HTML::HTML(R2HTML::as.title("Descriptive Statistics"), HR = 2, file = stdout())
			R2HTML::HTML(R2HTML::as.title("Mean and Standard Deviation"), HR = 3, file = stdout())
			DS <- t(sapply(c(res_var1,res_var2),function(i) return(c(mean(dataset[,i]),sd(dataset[,i])))))
			DS <- as.data.frame(DS)
			colnames(DS) <- c('Mean','SD')
			rownames(DS) <- c(res_var1,res_var2)
			R2HTML::HTML(DS, file = stdout(), row.names = T, col.names = TRUE, innerBorder = 1, align = "left",digits=6)

			R2HTML::HTML(R2HTML::as.title("Sample Correlation Matrix"), HR = 3, file = stdout())
			Cor <- cor(dataset[,c(res_var1,res_var2)])
			R2HTML::HTML(Cor, file = stdout(), row.names = T, col.names = TRUE, innerBorder = 1, align = "left",digits=6)

			# Test results
			R2HTML::HTML(R2HTML::as.title("Results of Paired Sample T-Test"), HR = 2, file =stdout())
			t.test.res.I <- as.data.frame(matrix(with(t.test.res,c(estimate,statistic,parameter,p.value)),nrow=1))
			colnames(t.test.res.I) <- c('Mean of the differences','T-value','DF','P-value')
			rownames(t.test.res.I) <- paste(res_var1,res_var2,sep=" & ")
			if(CI) {
				t.test.res.I <- cbind(t.test.res.I,t.test.res$conf.int[1],t.test.res$conf.int[2])
				data.frame(t.test.res.I,t.test.res$conf.int[1],t.test.res$conf.int[2])
				colnames(t.test.res.I)[5:6] <- 	c(paste0('Lower bound of ',conf.level*100,'% CI'),paste0('Upper bound of ',conf.level*100,'% CI'))
			}
			R2HTML::HTML(t.test.res.I, file = stdout(), col.names = TRUE, innerBorder = 1, align = "left",digits=5)
		}

		#### End of output
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Paired Sample T-test",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})
	
	return(html.output)
}
