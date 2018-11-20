REx_TSLS <- function(dataset,res_var,quan_var=NULL,qual_var=NULL,exp_var=NULL,ins_var=NULL,noint=FALSE,CI=TRUE,confint.level=0.95,ANOVA=FALSE,Predict=FALSE,Resid=FALSE,stdResid=FALSE){
	#### 변수 설명 ####
	## res_var : 종속변수 (필수 1개 선택되어야 함. 숫자형 변수만 가능)
	## qual_var : 질적변수 (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님)
	## quan_var : 양적변수 (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님)
	## exp_var : 설명변수 (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님)
	## ins_var : 도구변수 (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님)
	## noint : 모델설정 탭의 '상수항 포함하지 않음' 옵션, 이 옵션이 활성화되면 양적변수, 질적변수, 도구변수는 모두 선택되어야 함.
	## CI : 신뢰구간
	## confint.level : 신뢰수준
	## ANOVA : 출력옵션 탭의 'ANOVA Table' 옵션
	## Predict : 출력옵션 탭의 '예측값' 옵션
	## Resid : 출력옵션 탭의 '잔차' 옵션
	## stdResid: 출력옵션 탭의 '표준화잔차' 옵션
	###################
	#### Required packages : AER, R2HTML ####
	###################
	load.pkg("AER")
	options("scipen"=999,"digits"=4)

	html.output <- capture.output({
		# Title
		R2HTML::HTML(R2HTML::as.title("Two-Stage Least Square Regression"),HR=1,file=stdout(),append=FALSE)
	
		## Warnings
		# Response variable type
		if(!is.numeric(dataset[,res_var])) warn.msg1 <- '\a Error : Response variable should be numeric. Analysis has been stopped.'

		vars <- c(exp_var,ins_var)
		if(is.null(vars)) {
			Vars <- c()
		} else {
			Vars <-  unique(unlist(strsplit(vars,":")))
		}
		qual_var <- qual_var[qual_var%in%Vars]
		quan_var <- c(quan_var[quan_var%in%Vars],res_var)

		# explanatory & instrumental variable type
		if(!is.null(quan_var)) {
			is.nom <- sapply(quan_var,function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg2 <- paste0("\a Warning : The type of variable '",paste(quan_var[is.nom],collapse=', '),"' is not continuous but selected as the quantitative variable. It was excluded from the analysis.")
				ec <- quan_var[is.nom]
				# explanatory variables
				for(jj in ec) exp_var <- exp_var[-grep(jj,exp_var)]
				if(length(exp_var)==0) exp_var <- NULL
				# instrumental variables
				for(jj in ec) ins_var <- ins_var[-grep(jj,ins_var)]
				if(length(ins_var)==0) ins_var <- NULL
			}
		}
		if(!is.null(qual_var)) {
			is.num <- sapply(qual_var,function(i) is.numeric(dataset[,i]))
			if(any(is.num)) warn.msg3 <- paste0("\a Warning : The type of variable '",paste(qual_var[is.num],collapse=', '),"' is continuous but selected as the qualitative variable. It was coreced into character.")
			for(i in qual_var) dataset[,i] <- factor(dataset[,i])
		}

		# About explanatory & instrumental variables
		dim.exp <- ifelse(is.null(exp_var),0,ncol(model.matrix(formula(paste0(res_var,"~",paste(exp_var,collapse="+"))),data=dataset)))
		dim.ins <- ifelse(is.null(ins_var),0,ncol(model.matrix(formula(paste0(res_var,"~",paste(ins_var,collapse="+"))),data=dataset)))
		if(dim.exp>dim.ins) warn.msg4 <- paste("\aError : The number of explanatory variables should be less than that of instruments. Your analysis has ",dim.exp," explanatory variables and ",dim.ins," instruments. (All qualitative variables are transformed to dummy variables.)",sep="")

		# no intercept & no explanatory variable : error
		if(noint & is.null(exp_var)) warn.msg5 <- "\aError : With no intercept, at least 1 explanatory variable should be selected. Analysis has been stopped."

		# warnings
		if(exists('warn.msg1')|exists('warn.msg4')|exists('warn.msg5')){
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists('warn.msg5')){
				if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
				R2HTML::HTML(warn.msg5,file=stdout())
			}
			if(exists('warn.msg4')) R2HTML::HTML(warn.msg4,file=stdout())
		} else {
			### Model formula
			if(is.null(exp_var)){
				if(is.null(ins_var)){
					form.0 <- paste0(res_var,' ~ 1 | 1')
				} else {
					form.0 <- paste0(res_var,' ~ 1 | ',paste(ins_var,collapse=' + '))
				}
			} else {
				form.0 <- paste0(res_var,' ~ ',paste(exp_var,collapse=' + '),' | ',paste(ins_var,collapse=' + '))
			}
			if(noint) form.0 <- paste(gsub('\\|','-1 \\|',form.0),-1)
			form.1 <- as.formula(form.0)

			### Default output
			res_TSLS <- AER::ivreg(form.1,x=TRUE,data=dataset)

			## Data Structure
			R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout())
			total.var <- ncol(dataset)
			used.var <- ifelse(is.null(c(exp_var,ins_var)),0,length(unique(unlist(strsplit(c(exp_var,ins_var),":")))))
			DS <- matrix(c('Number of observations','Number of total variables','Number of used variables',nrow(dataset),total.var,used.var+1),ncol=2)
			R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")

			## Varibale list
			R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout())
			if(is.null(c(exp_var,ins_var))) {
				Vars <- c()
			} else {
				Vars <-  unique(unlist(strsplit(c(exp_var,ins_var),":")))
			}
			qual <- qual_var[qual_var%in%Vars]
			quan <- c(quan_var[quan_var%in%Vars],res_var)
			varlist <- matrix(c('Quantitative variable',paste(quan,collapse=', ')),ncol=2,byrow=T)
			if(length(qual)!=0) varlist <- rbind(varlist,c('Qualitative variable',paste(qual,collapse=', ')))
			R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")
			if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
			if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())

			## Analysis Description
			R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout())
			AD <- matrix(c('Response variable',res_var),ncol=2,byrow=T)
			if(!is.null(exp_var)) AD <- rbind(AD,c('Explanatory variable',paste(exp_var,collapse=', ')))
			if(!is.null(ins_var)) AD <- rbind(AD,c('Instrumental variable',paste(ins_var,collapse=', ')))
			AD <- rbind(AD,matrix(c('Intercept included',!noint),ncol=2,byrow=T))
			R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")

			## Analysis Results
			R2HTML::HTML(R2HTML::as.title("Results of Two-Stage Least Square Regression"),HR=2,file=stdout())

			# Coefficient Estimate
			R2HTML::HTML(R2HTML::as.title("Coefficients"),HR=3,file=stdout())
			if(!is.null(c(exp_var,ins_var))){
				CE <- summary(res_TSLS)$coef
				if(CI){
					CItab <- confint(res_TSLS, level=confint.level)
					CItab <- CItab[match(rownames(CE),rownames(CItab)),]
					CE <- cbind(CE,CItab)
					colnames(CE)[(ncol(CE)-1):ncol(CE)] <- paste0(c('Lower','Upper'),' bound of ',confint.level*100,'% CI')
				}
				R2HTML::HTML(round(CE,4),file=stdout(), innerBorder = 1,align="left")
			} else {
				R2HTML::HTML(res_TSLS$coef,file=stdout(), innerBorder = 1,align="left")
				R2HTML::HTML("\aWarning : Only estimate for intercept is provided for for intercept only model.",file=stdout())
			}
		
			# Anova table
			if(ANOVA){
				R2HTML::HTML(R2HTML::as.title("ANOVA Table"),HR=3,file=stdout())
				if(!noint){
					if(!is.null(c(exp_var,ins_var))){
						cc <- rownames(model.matrix(res_TSLS))		# Complete cases
						temp.dataset <- dataset[rownames(dataset)%in%cc,]
						Anova <- as.matrix(anova(AER::ivreg(formula(paste(res_var,'~1|1',sep='')),x=TRUE,data=temp.dataset),res_TSLS))
						A1 <- round(rbind(Anova[2,4:3], Anova[2:1,c(2,1)]),5)
						MS <- c(round(A1[1:2,1]/A1[1:2,2],5)," ")
						A2 <- data.frame(round(A1,5),MS=MS,F=c(round(Anova[2,5],5)," "," "),Pvalue=c(format(Anova[2,6],scientific=T)," "," "))
						colnames(A2)[1] <- 'SS'
						rownames(A2) <- c('SSR','SSE','SST')
						R2HTML::HTML(as.matrix(A2),file=stdout(), innerBorder = 1,align="left")
					} else {
						warn.msg6 <- "\a Warning : ANOVA is not supported for intercept only model."
						R2HTML::HTML(warn.msg6,file=stdout())
					}
				} else {
					warn.msg7 <- "\a Warning : ANOVA is not supported for no intercept model."
					R2HTML::HTML(warn.msg7,file=stdout())
				}
			}

			# 다음의 값들은 엑셀 시트에 저장
			if(Predict) {
				temp.0 <- data.frame(Fitted_TSLS=fitted(res_TSLS))
				temp <- data.frame(Fitted_TSLS=rep(NA,nrow(dataset)))
				row.dataset <- rownames(dataset)
				row.output <- rownames(temp.0)
				for(i in row.dataset) if(i%in%row.output) temp[row.dataset==i,1] <- temp.0[row.output==i,1]
				O <- temp
			}

			if(Resid) {
				temp.0 <- data.frame(Resid_TSLS=resid(res_TSLS))
				temp <- data.frame(Resid_TSLS=rep(NA,nrow(dataset)))
				row.dataset <- rownames(dataset)
				row.output <- rownames(temp.0)
				for(i in row.dataset) if(i%in%row.output) temp[row.dataset==i,1] <- temp.0[row.output==i,1]
				if(exists('O')){
					O <- cbind(O,temp)
				} else {
					O <- temp
				}
			}

			if(stdResid) {
				X <- as.matrix(res_TSLS$x$regressor)
				is.singular <- class(try(solve(t(X)%*%X),silent=T))=="try-error"
				if(is.singular){
					warn.msg6 <- "\a Error : Cannnot calculate standardized residuals because of the (near) linear dependences in explatory variables. Please check multicollinearity problem."
					R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
					R2HTML::HTML(warn.msg6,file=stdout())
				} else {
					H <- X%*%solve(t(X)%*%X)%*%t(X)
					oo <- res_TSLS$resid/(res_TSLS$sigma*sqrt(1-diag(H)))
					temp.0 <- data.frame(stdResid_TSLS=oo)
					temp <- data.frame(stdResid_TSLS=rep(NA,nrow(dataset)))
					row.dataset <- rownames(dataset)
					row.output <- rownames(temp.0)
					for(i in row.dataset) if(i%in%row.output) temp[row.dataset==i,1] <- temp.0[row.output==i,1]
					if(exists('O')){
						O <- cbind(O,temp)
					} else {
						O <- temp
					}
				}
			}
		}

		# Analysis end time
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Two-Stage Least Square Regression",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})

	if(exists('O')){
		return(list(html=html.output,Output=O))
	} else {
		return(html.output)
	}
}
