# 편최소제곱

REx_PLS <- function(dataset,dep_numeric_var=NULL,dep_cat_var=NULL,dep_cat_baseline=NULL,
indep_numeric_var=NULL,indep_cat_var=NULL,indep_cat_baseline=NULL,dimnum=2,crosval=TRUE,
Xscore=FALSE,Yscore=FALSE,Pred=FALSE,Resid=FALSE){		
			
	#### 변수 설명 ####
	## 종속변수는 missing값이 있으면 안됨(dep_numeric_var,dep_car_var)
	## dep_numeric_var : 질적종속변수(dependent variables)(2개 이상의 변수가 있을수 있음)
	## dep_cat_var : 양적종속변수(dependent variables)(2개 이상의 변수가 있을수 있음)
	## dep_cat_baseline : 기저범주
	## indep_numeric_var : 질적설명변수(independent variables) (2개 이상의 변수가 있는 경우 +로 구분, 필수사항아님,단,상수항옵션이(noint=TRUE)인경우 돌아갈수없음)
	## indep_cat_var : 양적설명변수(independent variables) (자동factor처리,2개 이상의 변수가 있는 경우 +로 구분, 필수사항아님,단,상수항옵션이(noint=TRUE)인경우 돌아갈수없음)
	## indep_cat_baseline : 기저범주
	## dimnum : 출력옵션에서 dimension 수 선택 최소 2
	## Xscore : 출력옵션(X-scroed (T-componets), default=FALSE)
	## Yscore : 출력옵션(Y-scroed (U-componets), default=FALSE)
	## Pred : 출력옵션(Y-predicted, default=FALSE)
	## Resid : 출력옵션(Residulas, default=FALSE)
	## 도표관련 옵션은 함수안에 코드로 구성되어있지만 UI구현 안함(Plot=FALSE이면 활성화 안됨)
	## Plot : 출력옵션 (Plot, default=FALSE)
	## plot_option : 출력옵션 (Plot=TRUE인경우, "observations" 또는 "variables"(default)중 선택할수 있음)
	## comp          : 출력옵션 (Plot=TRUE인경우,  X vs Y 축에 표현할 dimension의 숫자 선택가능 (숫자는 dimnum이상이 될 수 없음) )
	## comp (comp_x) : x축의 dimension의 숫자 (default=1<=dimnum)
	## comp (comp_y) : y축의 dimension의 숫자 (default=2<=dimnum)
	## where      	    : 출력옵션 (Plot=TRUE인경우,  X vs Y 축에 표현할 components의 조합)
	## where (wh_options) :  Possible options are: c("t","u") for using xy components
     	##                                     c("t","t"), for using x components (default)
      	##                                     c("u","u") for using y components. 
	## show.names (sn_options) : 출력옵션 (Plot=TRUE인경우, observation name으로 점 출력할 것인가? TRUE(default)/FALSE)
	###################
	#### Required packages : R2HTML , rms, MASS, plsdepot  ####
	###################
	load.pkg(c("R2HTML", "rms", "MASS", "plsdepot"))

	html.output <- capture.output({
		# Title
		R2HTML::HTML(R2HTML::as.title("Partial Least Squares (PLS) Data Analysis Methods"),HR=1,file=stdout(),append=FALSE)

		## Warnings


		## dependent variables
		## numeric
		numeric_var <- c(dep_numeric_var,indep_numeric_var)
		if(!is.null(numeric_var)) {
			is.nom <- sapply(numeric_var,function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg2 <- paste0("\a Warning : The type of variable '",paste(numeric_var[is.nom],collapse=', '),"' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")
				ec <- numeric_var[is.nom]
				dep_numeric_var <- dep_numeric_var[!dep_numeric_var%in%ec]
				indep_numeric_var <- indep_numeric_var[!indep_numeric_var%in%ec]
			}
		}
		if(length(dep_numeric_var)==0) dep_numeric_var <- NULL
		if(length(indep_numeric_var)==0) indep_numeric_var <- NULL

		##category
		cat_var <- c(dep_cat_var,indep_cat_var)
		if(!is.null(cat_var)) {
			is.num <- sapply(cat_var,function(i) is.numeric(dataset[,i]))
			if(any(is.num)) warn.msg3 <- paste0("\a Warning : The type of variable '",paste(cat_var[is.num],collapse=', '),"' is numeric but selected as the qualitative variable. It was coreced into character.")
			for(i in cat_var) dataset[,i] <- factor(dataset[,i])
		}

		if(is.null(dep_numeric_var)&is.null(dep_cat_var)) warn.msg1 <- '\a Error : At least 1 dependent variable should be selected. Analysis has been stopped.'
		if(is.null(indep_numeric_var)&is.null(indep_cat_var)) warn.msg4 <- '\a Error : At least 1 independent variable should be selected. Analysis has been stopped.'

		if(!crosval&is.null(dimnum)) warn.msg7 <- '\a Error : TRUE must be selected for cross validation when number of dimension is NULL. Analysis has been stopped.'
		
		# warnings
		if(exists('warn.msg1')|exists('warn.msg4')|exists('warn.msg7')){
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists('warn.msg4')) R2HTML::HTML(warn.msg4,file=stdout())
			if(exists('warn.msg7')) R2HTML::HTML(warn.msg7,file=stdout())
		} else {
			## Data Processing ## 
			var_info <- c(dep_numeric_var,dep_cat_var,indep_numeric_var, indep_cat_var)
			is.na <- sapply(var_info,function(i) is.na(dataset[,i]))
			raw.dat <- dataset
			oo <- which(rowSums(is.na)>0)
			if(length(oo)>0) dataset <- dataset[-oo,]


			if(!is.null(dep_cat_var)) {
				dep_o1 <- gsub("::.+","",dep_cat_baseline)
				dep_o2 <- gsub(".+::","",dep_cat_baseline)
				names(dep_o2) <- dep_o1
				sapply(dep_o1,function(i) relevel(factor(dataset[,i]),ref=dep_o2[i]))
			}
			if(!is.null(indep_cat_var)){
				indep_o1 <- gsub("::.+","",indep_cat_baseline)
				indep_o2 <- gsub(".+::","",indep_cat_baseline)
				names(indep_o2) <- indep_o1
				sapply(indep_o1,function(i) relevel(factor(dataset[,i]),ref=indep_o2[i]))
			}

			### 분석을 위해서 종속변수/독립변수 행렬 및 벡터는 무조건 numeric 형태여야함
			### 범주형변수는 numeric 0/1벡터로 바꿔서 범주가 n개이면 새로생성되는 벡터(변수)는 n-1개
			f <- function(i,dataset,cat,cat_baseline){
				level <- levels(dataset[,cat[i]])
				nn <- length(level)
				resout <- as.data.frame(matrix(0,nrow(dataset),nn))
				colnames(resout) <- paste(cat[i],"_",level,sep="")
				oo <- which(level==cat_baseline[i])		
				for(j in 1:nn){
					oo1 <- which(dataset[,cat[i]]==level[j])
					resout[oo1,j] <- 1
				}
				if(nn>2){
					resout <- resout[,-oo]
				}else{	
					resout <- as.data.frame(resout[,-oo])
					colnames(resout) <- paste(cat[i],"_",level[-oo],sep="")
				}
				resout
			}

			## 종속변수에 범주형 없을때
			if(is.null(dep_cat_var)){
				## 독립변수에 범주형 없을때
				if(is.null(indep_cat_var)){
					responses  <- as.data.frame(dataset[,dep_numeric_var])
					colnames(responses) <- dep_numeric_var
					predictors <- as.data.frame(dataset[,indep_numeric_var])				
					colnames(predictors) <- indep_numeric_var

				## 독립변수에 범주형 있을때
				}else{
					responses <- as.data.frame(dataset[,dep_numeric_var])
					colnames(responses) <- dep_numeric_var

					nn1 <- length(indep_cat_var)					
					if(nn1==1){
						indep_cat_new <- f(1,dataset,indep_o1,indep_o2)
					}else{
						indep_cat_new <- f(1,dataset,indep_o1,indep_o2)
						for(k in 2:nn1){
							indep_cat <- f(k,dataset,indep_o1,indep_o2)
							indep_cat_new <- cbind(indep_cat_new,indep_cat)
						}
					}
				
					## 독립변수에 연속형변수 여부에따라서 predictors 구성
					if(!is.null(indep_numeric_var)){
						predictors <- cbind(dataset[,indep_numeric_var],indep_cat_new)
						colnames(predictors) <- c(indep_numeric_var,colnames(indep_cat_new))
					}else{
						predictors <- indep_cat_new
					}
				}	
			## 종속변수에 범주형 있을때
			}else{
				nn2 <- length(dep_cat_var)
				if(nn2==1){
						dep_cat_new <- f(1,dataset,dep_o1,dep_o2)
				}else{
						dep_cat_new <- f(1,dataset,dep_o1,dep_o2)
						for(k in 2:nn2){
							dep_cat <- f(k,dataset,dep_o1,dep_o2)
							dep_cat_new <- cbind(dep_cat_new,dep_cat)
						}
				}
				## 종속변수에 연속형변수 여부에따라서 responses 구성
				if(!is.null(indep_numeric_var)){
					responses <- cbind(dataset[,dep_numeric_var],dep_cat_new)
					colnames(responses) <- c(dep_numeric_var,colnames(dep_cat_new))
				}else{
					responses <- dep_cat_new
				}

				## 독립변수에 범주형 없을때
				if(is.null(indep_cat_var)){
					predictors <- as.data.frame(dataset[,indep_numeric_var])				
					colnames(predictors) <- indep_numeric_var					
				## 독립변수에 범주형 있을때
				}else{
					nn3 <- length(indep_cat_var)					
					if(nn3==1){
						indep_cat_new <- f(1,dataset,indep_o1,indep_o2)
					}else{
						indep_cat_new <- f(1,dataset,indep_o1,indep_o2)
						for(k in 2:nn3){
							indep_cat <- f(k,dataset,indep_o1,indep_o2)
							indep_cat_new <- cbind(indep_cat_new,indep_cat)
						}
					}

				## 독립변수에 연속형변수 여부에따라서 predictors 구성
					if(!is.null(indep_numeric_var)){
						predictors <- cbind(dataset[,indep_numeric_var],indep_cat_new)
						colnames(predictors) <- c(indep_numeric_var,colnames(indep_cat_new))
					}else{
						predictors <- indep_cat_new
					}
				}
			}


			## responses : 종속변수 dataframe
			## predictors : 독립변수 dataframe

			## model
			form.0 <- paste0(paste(colnames(responses),collapse=' , '))
			#form <- ifelse(noint, paste0(form.0,' - 1'),form.0)
			form.1 <- paste0(paste(colnames(predictors),collapse=' , '))

			## Fitting PLS
			if(ncol(responses)==1){
				command_str	<- paste0("res_PLS_1<- plsreg1(predictors,responses,comps=dimnum,crosval=crosval)")
			}else{
				command_str	<- paste0("res_PLS_1<- plsreg2(predictors,responses,comps=dimnum,crosval=crosval)")
			}
			eval(parse(text=command_str)) ;
			#res_PLS_1 <- plsreg1(predictors,responses,comps=dimnum,crosval=crosval)
			#res_PLS_1 <- plsreg2(predictors,responses,comps=dimnum,crosval=crosval)

			### Default output
			# Data Structure
			R2HTML::HTML(R2HTML::as.title("Data structure"),HR=2,file=stdout())
			DS <- matrix(c('Number of observations','Number of total variables','Number of used variables',nrow(raw.dat),ncol(dataset),ncol(responses)+ncol(predictors)),ncol=2)
			R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")

			# Varaible List
			R2HTML::HTML(R2HTML::as.title("Variable list"),HR=2,file=stdout())
			VL <- matrix(0,1,2)
			if(!is.null(c(dep_numeric_var,indep_numeric_var))) VL <- rbind(VL,c('Quantitative variable',paste(c(dep_numeric_var,indep_numeric_var),collapse=', ')))
			if(!is.null(c(dep_cat_var,indep_cat_var))) VL <- rbind(VL,c('Quantitative variable',paste(c(dep_cat_var,indep_cat_var),collapse=', ')))
			VL <- VL[-1,]
			R2HTML::HTML(VL,file=stdout(), innerBorder = 1,align="left")
			if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
			if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())
			
			# Analysis description
			R2HTML::HTML(R2HTML::as.title("Analysis description"),HR=2,file=stdout())
			AD <- matrix(c('Reponse variable','Explanatory variable',form.0,form.1),ncol=2)
			R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")

			#### Results
			R2HTML::HTML(R2HTML::as.title("Results of Partial Least Squares Data Analysis Methods"),HR=2,file=stdout())

			#################### response의 열이 하나인경우 출력결과(plsreg1) ##############
			if(ncol(responses)==1){
				####### Results of Models
				Xscores <- as.data.frame(res_PLS_1$x.scores)
				Xload   <- as.data.frame(res_PLS_1$x.loads)
				Yscores <- as.data.frame(res_PLS_1$y.scores)
				rownames(Yscores) <- rownames(Xscores)
				Yload   <- as.data.frame(res_PLS_1$y.loads)

				corxyt <- as.data.frame(res_PLS_1$cor.xyt)
			
				raw_weight <- as.data.frame(res_PLS_1$raw.wgs)
				mod_weight <- as.data.frame(res_PLS_1$mod.wgs)

				std <- as.data.frame(res_PLS_1$std.coefs)
				regcoeff <- as.data.frame(res_PLS_1$reg.coefs)	

				pred <- as.data.frame(res_PLS_1$y.pred)
				rownames(pred) <- rownames(Xscores)
				resid <- as.data.frame(res_PLS_1$resid)
				rownames(resid) <- rownames(Xscores)
				R2 <- as.data.frame(res_PLS_1$R2)
				R2xy <- as.data.frame(res_PLS_1$R2Xy)
				T2 <- as.data.frame(res_PLS_1$T2)
				Q2 <- as.data.frame(res_PLS_1$Q2)

				R2HTML::HTML(R2HTML::as.title("X-Loadings"),HR=3,file=stdout())
				R2HTML::HTML(Xload,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Y-Loadings"),HR=3,file=stdout())
				R2HTML::HTML(Yload,file=stdout(), innerBorder = 1,align="left")	

				R2HTML::HTML(R2HTML::as.title("Score correlation"),HR=3,file=stdout())
				R2HTML::HTML(corxyt,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Raw Weights"),HR=3,file=stdout())
				R2HTML::HTML(raw_weight,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Modified Weights"),HR=3,file=stdout())
				R2HTML::HTML(mod_weight,file=stdout(), innerBorder = 1,align="left")	
	
				R2HTML::HTML(R2HTML::as.title("Regular Coefficients"),HR=3,file=stdout())
				R2HTML::HTML(regcoeff,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Standard Coefficients"),HR=3,file=stdout())
				R2HTML::HTML(std,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("R-squared"),HR=3,file=stdout())
				R2HTML::HTML(R2,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Explained Variance of X-Y by T"),HR=3,file=stdout())
				R2HTML::HTML(R2xy,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("T2 hotelling"),HR=3,file=stdout())
				R2HTML::HTML(T2,file=stdout(), innerBorder = 1,align="left")

				if(crosval==TRUE){
					R2HTML::HTML(R2HTML::as.title("Q2 cross validation"),HR=3,file=stdout())
					R2HTML::HTML(Q2,file=stdout(), innerBorder = 1,align="left")
				}


			#################### response의 열이 하나인경우 출력결과(plsreg2) ##############
			}else{			
				####### Results of Models
				Xscores <- as.data.frame(res_PLS_1$x.scores)
				Xload   <- as.data.frame(res_PLS_1$x.loads)
				Yscores <- as.data.frame(res_PLS_1$y.scores)
				rownames(Yscores) <- rownames(Xscores)
				Yload   <- as.data.frame(res_PLS_1$y.loads)

				corxt <- as.data.frame(res_PLS_1$cor.xt)
				coryt <- as.data.frame(res_PLS_1$cor.yt)
				corxu <- as.data.frame(res_PLS_1$cor.xu)
				coryu <- as.data.frame(res_PLS_1$cor.yu)
				cortu <- as.data.frame(res_PLS_1$cor.tu)			

				raw_weight <- as.data.frame(res_PLS_1$raw.wgs)
				mod_weight <- as.data.frame(res_PLS_1$mod.wgs)

				std <- as.data.frame(res_PLS_1$std.coefs)
				regcoeff <- as.data.frame(res_PLS_1$reg.coefs)	

				pred <- as.data.frame(res_PLS_1$y.pred)
				rownames(pred) <- rownames(Xscores)
				resid <- as.data.frame(res_PLS_1$resid)
				rownames(resid) <- rownames(Xscores)			
				exp_var <- as.data.frame(res_PLS_1$expvar)
				vip <- as.data.frame(res_PLS_1$VIP)
				q2 <- as.data.frame(res_PLS_1$Q2)
				q2sum <- as.data.frame(res_PLS_1$Q2cum)

				R2HTML::HTML(R2HTML::as.title("X-Loadings"),HR=3,file=stdout())
				R2HTML::HTML(Xload,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Y-Loadings"),HR=3,file=stdout())
				R2HTML::HTML(Yload,file=stdout(), innerBorder = 1,align="left")	

				R2HTML::HTML(R2HTML::as.title("(X,T) correlation"),HR=3,file=stdout())
				R2HTML::HTML(corxt,file=stdout(), innerBorder = 1,align="left")
	
				R2HTML::HTML(R2HTML::as.title("(Y,T) correlation"),HR=3,file=stdout())
				R2HTML::HTML(coryt,file=stdout(), innerBorder = 1,align="left")	
	
				R2HTML::HTML(R2HTML::as.title("(X,U) correlation"),HR=3,file=stdout())
				R2HTML::HTML(corxu,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("(Y,U) correlation"),HR=3,file=stdout())
				R2HTML::HTML(coryu,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("(T,U) correlation"),HR=3,file=stdout())
				R2HTML::HTML(cortu,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Raw Weights"),HR=3,file=stdout())
				R2HTML::HTML(raw_weight,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Modified Weights"),HR=3,file=stdout())
				R2HTML::HTML(mod_weight,file=stdout(), innerBorder = 1,align="left")	
	
				R2HTML::HTML(R2HTML::as.title("Regular Coefficients"),HR=3,file=stdout())
				R2HTML::HTML(regcoeff,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Standard Coefficients"),HR=3,file=stdout())
				R2HTML::HTML(std,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Explained Variance"),HR=3,file=stdout())
				colnames(exp_var) <- c("X variance","Cumulative X variance","Y variance","Cummulative Y variance (R-square)")
				R2HTML::HTML(exp_var,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Variable Importance for Projection(VIP)"),HR=3,file=stdout())
				R2HTML::HTML(vip,file=stdout(), innerBorder = 1,align="left")
				if(crosval==TRUE){
					R2HTML::HTML(R2HTML::as.title("Q2 Index"),HR=3,file=stdout())
					colnames(q2) <- c(colnames(responses),"Total")
					R2HTML::HTML(q2,file=stdout(), innerBorder = 1,align="left")
					R2HTML::HTML(R2HTML::as.title("Cummulated Q2"),HR=3,file=stdout())
					colnames(q2sum) <- c(colnames(responses),"Total")
					R2HTML::HTML(q2sum,file=stdout(), innerBorder = 1,align="left")
				}

			}

			# 다음의 값들은 엑셀 시트에 저장
			if(Xscore) {
				temp.0 <- Xscores
				temp <- as.data.frame(matrix(NA,nrow(raw.dat),ncol(Xscores)))
				colnames(temp) <- paste("Xscore_",colnames(temp.0),sep="")
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(dataset)
				for(i in row.dataset) if(i%in%row.output) temp[i,] <- temp.0[row.output==i,]
				if(exists('O')) {
					O <- cbind(O,temp)
				} else {
					O <- temp
				}
			}

			if(Yscore) {
				temp.0 <- Yscores
				temp <- as.data.frame(matrix(NA,nrow(raw.dat),ncol(Yscores)))
				colnames(temp) <- paste("Yscore_",colnames(temp.0),sep="")
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(dataset)
				for(i in row.dataset) if(i%in%row.output) temp[i,] <- temp.0[row.output==i,]
				if(exists('O')) {
					O <- cbind(O,temp)
				} else {
					O <- temp
				}
			}

			if(Pred) {
				temp.0 <- pred
				temp <- as.data.frame(matrix(NA,nrow(raw.dat),ncol(pred)))
				colnames(temp) <- paste("Pred_",colnames(responses),sep="")
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(dataset)
				for(i in row.dataset) if(i%in%row.output) temp[i,] <- temp.0[row.output==i,]
				if(exists('O')) {
					O <- cbind(O,temp)
				} else {
					O <- temp
	
				}
			}
			if(Resid) {
				temp.0 <- resid
				temp <- as.data.frame(matrix(NA,nrow(raw.dat),ncol(resid)))
				colnames(temp) <- paste("Resid_",colnames(responses),sep="")
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(dataset)
				for(i in row.dataset) if(i%in%row.output) temp[i,] <- temp.0[row.output==i,]
				if(exists('O')) {
					O <- cbind(O,temp)
				} else {
					O <- temp
	
				}
			}
		}
		# Analysis end time
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Partial Least Squares Analysis.",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})

	if(exists('O')){
		return(list(html=html.output,Output=O))
	} else {
		return(html.output)
	}
}
