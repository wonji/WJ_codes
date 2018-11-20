REx_MANOVA <- function(dataset,res_var,qual_var=NULL,vars=NULL,noint=FALSE,ss="I",test_mode="Pillai"){

	#### 변수 설명 ####
	## dataset : 데이터셋이름
	## res_var : 종속변수, 필수로 입력되어야 하며 2개 이상 입력가능
	## qual_var : 질적설명변수(independent variables) (반드시 한개 이상 입력. 2개 이상의 변수가 있는 경우 c(,,)로 구분)
	## vars : 최종모형에 선택된 변수들 (2개 이상의 변수가 있는 경우 c(,,)로 구분,교호작용의 경우 :를 이용하여 입력, 필수사항)	#수정
	## noint : 모형에 절편한 포함 여부. TRUE or FALSE.
	## ss: 제곱합 유형 지정. I,II,III 유형 지정가능
	## test_mode : c("Pillai", "Hotelling-Lawley", "Wilks", "Roy") 의 값을 가질 수 있음.
	## ref : http://users.stat.umn.edu/~helwig/notes/mvlr-Notes.pdf, http://ibgwww.colorado.edu/~carey/p7291dir/handouts/manova1.pdf
	###################
	#### Required packages : R2HTML, car, MASS ####
	###################
	load.pkg(c("R2HTML", "car", "MASS"))

	global_op<-options()	#수정
	options(contrasts=c("contr.sum", "contr.poly"))	#수정

	html.output <- capture.output({
		# Title
		R2HTML::HTML(R2HTML::as.title("Multivariate Analysis of Variance"),HR=1,file=stdout(),append=FALSE)

		## Warnings
		if(any(sapply(res_var,function(i) !is.numeric(dataset[,i])))) warn.msg1 <- '\a Error : Response variable should be numeric. Analysis has been stopped.'
		if(is.null(qual_var)) warn.msg4 <- '\a Error : At least 1 group variable should be selected. Analysis has been stopped.'	#수정
		if(!is.null(qual_var)) {
			is.num <- sapply(qual_var,function(i) is.numeric(dataset[,i]))
			if(any(is.num)) warn.msg5 <- paste0("\a Warning : The type of variable '",paste(qual_var[is.num],collapse=', '),"' is numeric but selected as the qualitative variable. It was coreced into character.")
		}

		if(exists('warn.msg1')|exists('warn.msg4')){
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists('warn.msg4')) R2HTML::HTML(warn.msg4,file=stdout())
		
		} else {

			## Model formula
			y<-noquote(paste("cbind(",paste(res_var,collapse=','),")"))
			form.0 <- paste0("cbind(",paste(res_var,collapse=','),") ~ ",paste(vars,collapse=' + '))	#추가수정
			form <- as.formula(paste(form.0,ifelse(noint,'-1',''),sep=''))	#수정
			form.1 <- gsub('cbind','',form.0)
			newdat <- dataset[,c(res_var,qual_var),drop=F]
			newdat <- newdat[complete.cases(newdat),,drop=F]
			if(!is.null(qual_var)) for(i in qual_var) newdat[,i] <- factor(newdat[,i])

			# Data Structure
			R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout())
			DS <- matrix(c('Number of observations','Number of response variables','Number of group variables',nrow(newdat),length(res_var),length(qual_var)),ncol=2)
			R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")
			if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())

			# Varibale list
			R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout())
			varlist <- matrix(c('Response variables',paste(res_var,collapse=', '),'Group variable',paste(qual_var,collapse=', ')),ncol=2,byrow=T)
			R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")
			if(exists('warn.msg5')) R2HTML::HTML(warn.msg5,file=stdout())

			# Analysis Description
			R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout())
			AD <- matrix(c('Fitted model',ifelse(noint,paste0(gsub(" - 1","",form.1),' (Intercept is not included)'),form.1),'Test statistics',paste(test_mode,collapse=', '),'Sums of Squares',ss),ncol=2,byrow=T)
			R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")

			# Scatter plot matrix
			R2HTML::HTML(R2HTML::as.title("Scatter Plots"),HR=2,file=stdout()) 
			n <- ifelse(length(res_var)<4,4,length(res_var))
			for(i in qual_var){
				REx_ANA_PLOT(w=130*n,h=130*n)
				form.2 <- as.formula(paste0('~',paste(res_var,collapse='+'),'|',i))
				check.ellipse <- function(j) apply(newdat[newdat[,i]==j,res_var],2,function(ii) length(unique(ii))==1)
				if(any(sapply(levels(newdat[,i]),check.ellipse))) {					
					scatterplotMatrix(form.2,data=newdat, smooth=FALSE, reg.line=FALSE, by.groups=TRUE, diagonal="none",main=i)
					cap <- paste0("<div style='text-align:left'> \a Due to constant values in some groups in variable '",i,"' there is no way to draw an ellipse.")
				} else {
					scatterplotMatrix(form.2,data=newdat, smooth=FALSE, reg.line=FALSE, ellipse=TRUE, by.groups=TRUE, diagonal="none",main=i)
					cap <- ""
				}
				REx_ANA_PLOT_OFF(cap)
			}
	
			## Model Fitting
			R2HTML::HTML(R2HTML::as.title("Results of Multivariate Analysis of Variance"),HR=2,file=stdout())
			
			# Coefficient estimates
			R2HTML::HTML(R2HTML::as.title("Coefficient Estimates"),HR=3,file=stdout())
			mlm.fit <- lm(form, data=newdat)
			for(i in 1:length(res_var)){
				R2HTML::HTML(R2HTML::as.title(paste0(res_var[i])),HR=4,file=stdout())
				CE <- as.data.frame(summary(mlm.fit)[[i]]$coef)
				CE[,ncol(CE)] <- format(CE[,ncol(CE)],scientific=T,digits=4)
				R2HTML::HTML(CE,file=stdout(), innerBorder = 1,align="left",digits=4)
			}
				
			# ANOVA Table
			R2HTML::HTML(R2HTML::as.title(paste0("Multivariate ANOVA Table with Type ",ss," SS")),HR=3,file=stdout())
			if(ss=='I'){
				res<-try(manova(form, data = newdat),s=T)
				if(class(res)[1]!='try-error'){
					ANOVA <- summary(res,intercept=!noint)$stat[,-2]	#수정
					for(i in test_mode) ANOVA <- cbind(summary(res,intercept=!noint, test=i)$stat[,2,drop=F],ANOVA)	#수정
					if(nrow(ANOVA)==2){ANOVA <- as.data.frame(t(ANOVA[-nrow(ANOVA),]));rownames(ANOVA)<-qual_var}	#수정
					if(nrow(ANOVA)>2){ANOVA <- as.data.frame(ANOVA[-nrow(ANOVA),])}		#수정
					ANOVA[,ncol(ANOVA)] <- format(ANOVA[,ncol(ANOVA)],scientific=T,digits=4)
					R2HTML::HTML(ANOVA,file=stdout(), innerBorder = 1,align="left",digits=4)
				} else {
					warn.msg6 <- "\a Error : Fail to fit the multivariate Analysis of Variance."
					R2HTML::HTML(warn.msg6,file=stdout())
				}
			}

			if(ss %in% c('II','III')){
				fit<-try(lm(form,data=newdat),silent=T)
				if(ss=='II') res<-try(car::Manova(fit,type='II'))
				if(ss=='III'){res<-try(car::Manova(fit,type='III'))}
				if(class(res)[1]!='try-error'){
					outtests <- car:::print.Anova.mlm
					body(outtests)[[16]] <- quote(invisible(tests))
					body(outtests)[[15]] <- NULL
					ANOVA <- outtests(car::Manova(fit, type=ss, test.statistic="Pillai"))[,-2]	
					for(i in test_mode) {
						ANOVA <- cbind(outtests(car::Manova(fit, type=ss, test.statistic=i))[,2,drop=F],ANOVA)
						colnames(ANOVA)[1] <- i
					}
					ANOVA[,ncol(ANOVA)] <- format(ANOVA[,ncol(ANOVA)],scientific=T,digits=4)
					R2HTML::HTML(ANOVA,file=stdout(), innerBorder = 1,align="left",digits=4)
				} else {
					warn.msg7 <- "\a Error : Fail to fit the multivariate Analysis of Variance."
					R2HTML::HTML(warn.msg7,file=stdout())
				}
			}
		}

		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Multivariate Analysis of Variance",sep=""),stdout())
		R2HTML::HTMLhr(file=stdout())
	})

	options(global_op)	#수정
	return(html.output)
}
