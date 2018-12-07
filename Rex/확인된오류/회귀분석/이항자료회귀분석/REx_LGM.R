# 이항자료회귀분석
REx_LGM <- function(dataset,dep_var,dep_ref=NULL,indep_cat_var=NULL,indep_numeric_var=NULL,vars=NULL,indep_ref=NULL,noint=FALSE,link='logit',
			Valid.method=c('Partition','Cross'),Part.method=c('all','percent','variable'),train.perc=70,Part.var=NULL, Cross.method=c('LOOCV','KFOLD'),k=10,Pred.var=NULL, Cutoff=0.5,
			CI=TRUE,confint.level=0.95,odds=FALSE,VIF=FALSE,ANOVA=TRUE,ss="III",classtab=FALSE,GOF=FALSE,VIP=FALSE,ROCcurve=FALSE,
			Predict_prob_train=FALSE, Predict_CI_train=FALSE,confint.level_train=0.95,Predict_g_train=FALSE,
			Resid=FALSE,stdResid=FALSE,studResid=FALSE,cook_distance=FALSE,linear_pred=FALSE,hat_value=FALSE,
			Predict_prob_test=FALSE,Predict_CI_test=FALSE,confint.level_test=0.95,Predict_g_test=FALSE,
			Predict_CI_pred=FALSE,confint.level_pred=0.95,Part_index=FALSE,select=FALSE,direct="forward",keep_var=NULL){
		
	load.pkg(c("R2HTML", "rms", "MASS", "caret", "ResourceSelection", "MKmisc", "lmtest", "VGAM", "ROCR", "ggplot2", "ggfortify","pscl"))

	orig.dataset <- dataset 
	if(length(dep_var) > 1) 
		Predict_prob_train <- Predict_CI_train <- Predict_g_train <- Resid <- stdResid <- studResid <- cook_distance <- linear_pred <- hat_value <- Predict_prob_test <- Predict_CI_test <- Predict_g_test <- Predict_CI_pred <- Part_index <- FALSE

	Dep_var <- dep_var
	FIN.RESULT <- rep(NA, length(Dep_var))
	
	for(dep_var in Dep_var){
		dataset <- orig.dataset
		if(length(Dep_var) > 1) 
			dep_ref <- NULL

		html.output <- capture.output({
			# Title
			R2HTML::HTML(R2HTML::as.title("Generalized Linear Regression for Binomial Data"),HR=1,file="./test.html",append=FALSE)

			## Warnings
			# Response variable type
			dataset[,dep_var] <- as.factor(dataset[,dep_var])
			if(is.null(dep_ref)) dep_ref <- levels(dataset[,dep_var])[1]
			if(length(levels(factor(dataset[,dep_var])))!=2) warn.msg1 <- '<li> Error : Dependent variable should be binary. Analysis has been stopped.'		
			
			# explanatory variable type
			if(is.null(vars)) {
				Vars <- c()
			} else {
				Vars <-  unique(unlist(strsplit(vars,":")))
			}
			indep_cat_var <- indep_cat_var[indep_cat_var%in%Vars]
			indep_numeric_var <- indep_numeric_var[indep_numeric_var%in%Vars]
			if(length(indep_numeric_var)==0) indep_numeric_var <- NULL
			if(length(indep_cat_var)==0) indep_cat_var <- NULL

			if(!is.null(indep_numeric_var)) {
				is.nom <- sapply(indep_numeric_var,function(i) !is.numeric(dataset[,i]))
				if(any(is.nom)){
					warn.msg2 <- paste0("<li> Warning : The type of variable '",paste(indep_numeric_var[is.nom],collapse=', '),"' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")
					ec <- indep_numeric_var[is.nom]
					for(jj in ec) vars <- vars[-grep(jj,vars)]
					if(length(vars)==0) vars <- NULL

				}
			}
			if(!is.null(indep_cat_var)) {
				is.num <- sapply(indep_cat_var,function(i) is.numeric(dataset[,i]))
				if(any(is.num)) warn.msg3 <- paste0("<li> Warning : The type of variable '",paste(indep_cat_var[is.num],collapse=', '),"' is numeric but selected as the qualitative variable. It was coerced into character.")
				for(i in indep_cat_var[is.num]) dataset[,i] <- as.factor(dataset[,i])
			}

			if(is.null(vars)) {
				Vars <- c()
			} else {
				Vars <-  unique(unlist(strsplit(vars,":")))
			}
			indep_cat_var <- c(indep_cat_var[indep_cat_var%in%Vars],dep_var)
			indep_numeric_var <- indep_numeric_var[indep_numeric_var%in%Vars]
			if(length(indep_numeric_var)==0) indep_numeric_var <- NULL
			if(length(indep_cat_var)==0) indep_cat_var <- NULL

			# no intercept & no explanatory variable : error
			if(noint&is.null(vars)) warn.msg4 <- '<li> Error : With no intercept, at least 1 independent variable should be selected. Analysis has been stopped.'

			# Warning in VS
			if(select) {
				if(is.null(vars)) {
					select <- FALSE
					warn.VS <- '<li> Warning : Variable selection is not supported for intercept only model.'
				}
				keep_var <- keep_var[keep_var%in%vars]
				if(length(keep_var)==0) keep_var <- NULL
			}
			
			## Data processing
			original.dataset <- dataset
			
			# predict variable
			if(!is.null(Pred.var)){
				pred.level <- na.omit(unique(dataset[,Pred.var]))
				# No training & testing dataset
				if(!2%in%pred.level) {
					warn.DP1 <- paste0("<li> Error : '2' is not observed in the split variable for prediction, '",Pred.var,"'. That is, no observations are assigned to the training & testing dataset. Analysis has been stopped.")
				} else {
					if(!1%in%pred.level) {	
						# Only 2 is obaserved
						warn.DP3 <- paste0("<li> Warning : '1' is not observed in the split variable for prediction, '",Pred.var,"'. That is, no observations are assigned to the prediction dataset. Prediction is not supported in this analysis.")
						Predict_CI_pred <- FALSE
					} else if(!all(pred.level%in%c(1,2))) {
						# Value other than 1 and 2
						warn.DP2 <- paste0("<li> Warning : Values other than '1' and '2' are observed in the split variable for prediction, '",Pred.var,"'. Observations with the value other than '1' and '2' are not included the analysis.")
						pred.dataset <- dataset[dataset[,Pred.var]==1,]	# prediction dataset
						dataset <- dataset[dataset[,Pred.var]==2,]	# training & testing dataset
					}
				}
			}
			
			# training & testing
			if(Valid.method=='Partition'){
				if(Part.method=='percent'){
					data.case <- dataset[dataset[,dep_var]!=dep_ref & !is.na(dataset[,dep_var]),,drop=F]
					data.cont <- dataset[dataset[,dep_var]==dep_ref & !is.na(dataset[,dep_var]),,drop=F]
					n.train.case <- round(nrow(data.case)*train.perc/100)
					n.train.cont <- round(nrow(data.cont)*train.perc/100)
					
					if(train.perc==100 | (n.train.case==nrow(data.case) & n.train.cont==nrow(data.cont))){
						warn.DP4 <- "<li> Warning : No observations are assigned to the testing dataset due to too high percent for the training dataset. Validation using testing dataset is not supported in this analysis."
						Predict_prob_test <- Predict_CI_test <- Predict_g_test <- FALSE
					} else if(n.train.case==0 | n.train.cont==0){
						warn.DP5 <- "<li> Error : No observations are assigned to the training dataset due to too low percent for the training dataset. To secure sufficient number of observations for the training dataset, increase the percent. Analysis has been stopped."
					} else {
						train.case <- sample(seq(nrow(data.case)),round(nrow(data.case)*train.perc/100))
						train.cont <- sample(seq(nrow(data.cont)),round(nrow(data.cont)*train.perc/100))

						# training set
						dataset <- rbind(data.case[train.case,],data.cont[train.cont,])
						
						# testing set
						test.dataset <- rbind(data.case[-train.case,],data.cont[-train.cont,])
						if(nrow(original.dataset)!=(nrow(dataset)+nrow(test.dataset)))	warn.DP9 <- "<li> Warning : Observations with missing dependent variable were not assigned to training dataset or testing dataset."
						
						# ordering
						dataset <- dataset[order(as.numeric(rownames(dataset))),,drop=F]
						test.dataset <- test.dataset[order(as.numeric(rownames(test.dataset))),,drop=F]

					}
				}
				if(Part.method=='variable'){
					part.level <- na.omit(unique(dataset[,Part.var]))
					# No training & test dataset
					if(!1%in%part.level) {
						warn.DP6 <- paste0("<li> Error : '1' is not observed in the split variable for validation, '",Part.var,"'. That is, no observations are assigned to the training dataset. Analysis has been stopped.")
					} else if(!2%in%part.level) {
						warn.DP7 <- paste0("<li> Warning : '2' is not observed in the split variable for validation, '",Part.var,"'. That is, no observations are assigned to the testing dataset. Validation using testing dataset is not supported in this analysis.")
						Predict_prob_test <- Predict_CI_test <- Predict_g_test <- FALSE
					} else {
						# Value other than 1 and 2
						if(!all(part.level%in%c(1,2))) {
							warn.DP8 <- paste0("<li> Warning : Values other than '1' and '2' are observed in the split variable for validation, '",Part.var,"'. Observations with the value other than '1' and '2' are not included the analysis.")
						}
						test.dataset <- dataset[dataset[,Part.var]==2,]	# testing dataset
						dataset <- dataset[dataset[,Part.var]==1,]	# training dataset
					}
				}
			}
			#### Then, original.dataset : original dataset, dataset : training dataset, test.dataset : testing dataset, pred.dataset : prediction dataset.

			# warnings
			if(exists('warn.msg1')|exists('warn.msg4')|exists('warn.DP1')|exists('warn.DP5')|exists('warn.DP6')){
				R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file="./test.html")
				if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file="./test.html")
				if(exists('warn.msg4')){
					if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file="./test.html")
					R2HTML::HTML(warn.msg4,file="./test.html")
				}
				if(exists('warn.DP1')) R2HTML::HTML(warn.DP1,file="./test.html")
				if(exists('warn.DP5')) R2HTML::HTML(warn.DP5,file="./test.html")
				if(exists('warn.DP6')) R2HTML::HTML(warn.DP6,file="./test.html")
			} else {
				## Data Processing
				var_info <- c(dep_var,indep_numeric_var, indep_cat_var)
				temp.dat <- dataset[complete.cases(dataset[,var_info,drop=F]),,drop=F]
				## To define baseline for dep/indep cat
				if(!is.null(dep_ref)) {
					temp.dat[,dep_var] <- relevel(factor(temp.dat[,dep_var]),ref=as.character(dep_ref))
					if(exists('test.dataset')) test.dataset[,dep_var] <- relevel(factor(test.dataset[,dep_var]),ref=as.character(dep_ref))
					if(exists('pred.dataset')) pred.dataset[,dep_var] <- relevel(factor(pred.dataset[,dep_var]),ref=as.character(dep_ref))
				}
				## model formula
				form.0 <- ifelse(is.null(vars),paste0(dep_var,' ~ 1'),paste0(dep_var,' ~ ',paste(vars,collapse=' + ')))
				if(noint) form.0 <- paste(form.0,-1)
				form.1 <- as.formula(form.0)
				
				## Fitting LGM
				command_str	<- paste0("res_LGM_1<- try(glm(formula(form.1), data=temp.dat, family=binomial(",link,")))")
				eval(parse(text=command_str)) ;
				res_LGM <- res_LGM_1
				res_LGM <<- res_LGM
				if('try-error'%in%class(res_LGM)){
					R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file="./test.html")
					R2HTML::HTML("<li> Warning : Generalized linear regression for binomial data was failed to fit. Analysis has been stopped.",file="./test.html")
				} else if(!res_LGM$conv) {
					R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file="./test.html")
					warn.msg.conv <- "<li> Warning : Generalized linear regression for binomial data did not converge. Analysis has been stopped."
					R2HTML::HTML(warn.msg.conv,file="./test.html")
				} else {
					## confidence & prediction interval
					getIntervals <- function(LGM_obj,newdat,con.level){
						# CI for linear predictor
						linpred <- predict(LGM_obj,newdata=newdat,se.fit=T)
						lower <- linpred$fit-qnorm((1-con.level)/2,lower.tail=F)*linpred$se.fit
						upper <- linpred$fit+qnorm((1-con.level)/2,lower.tail=F)*linpred$se.fit

						# Get CI for fitted value by transforming the CI for linear predictor
						LinkFt <- function(x,link){
							tf <- switch(link, logit=exp(x)/(1+exp(x)), probit=pnorm(x,lower.tail=T), cauchit=(atan(x)+pi/2)/pi, cloglog=1-exp(-exp(x)))
							return(tf)
						}
						t.lower <- LinkFt(lower,link); t.lower[t.lower<0] <- 0
						t.upper <- LinkFt(upper,link); t.upper[t.upper>1] <- 1
						t.CI <- cbind(t.lower,t.upper)
						return(t.CI)
					}
			
					## 다음의 값들은 엑셀 시트에 저장
					# Training dataset
					if(Predict_prob_train) {
						temp.0 <- data.frame(Fitted_LGM=predict(res_LGM,newdata=temp.dat,type='response'))
						temp <- data.frame(Fitted_LGM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
						mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
						mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
						O <- data.frame(Fitted_train_LGM=mer.1[,3])
					}
					if(Predict_CI_train){
						temp.0 <- as.data.frame(getIntervals(LGM_obj=res_LGM,newdat=temp.dat,con.level=confint.level_train))
						temp <- data.frame(temp=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
						mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
						mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),c((ncol(mer)-1),ncol(mer))]
						colnames(mer.1) <- paste0('Fitted_',round(confint.level_train*100,0),'CI_',c('Lower','Upper'),'_train_LGM')
					
						if(exists('O')) {
							O <- cbind(O,mer.1)
						} else {
							O <- mer.1
						}
					}
					if(Predict_g_train) {
						temp.0 <- data.frame(Fitted_group_LGM=ifelse(predict(res_LGM,newdata=temp.dat,type='response')<Cutoff,levels(temp.dat[,dep_var])[1],levels(temp.dat[,dep_var])[2]))
						temp.0[,1] <- as.character(temp.0[,1])
						temp <- data.frame(Fitted_group_LGM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
						mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
						mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
						if(exists('O')) {
							O <- cbind(O,Fitted_group_train_LGM=mer.1[,3],stringsAsFactors=F)
						} else {
							O <- data.frame(Fitted_group_train_LGM=mer.1[,3],stringsAsFactors=F)
						}
					}
					if(Resid) {
						temp.0 <- data.frame(Resid_LGM=resid(res_LGM))
						temp <- data.frame(Resid_LGM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
						mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
						mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
						if(exists('O')) {
							O <- cbind(O,unstdResid_train_LGM=mer.1[,3])
						} else {
							O <- data.frame(unstdResid_train_LGM=mer.1[,3])
						}
					}
					if(stdResid) {
						temp.0 <- data.frame(stdResid_LGM=stdres(res_LGM))
						is.singular <- "try-error"%in%class(try(temp.0[,1],silent=T))
						if(is.singular){
							warn.res1 <- "<li> Error : Cannnot calculate standardized residuals because of the (near) linear dependences in explatory variables. Please check multicollinearity problem."
						} else {
							temp <- data.frame(stdResid_LGM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
							mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
							mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
							if(exists('O')) {
								O <- cbind(O,stdResid_train_LGM=mer.1[,3])
							} else {
								O <- data.frame(stdResid_train_LGM=mer.1[,3])
							}
						}
					}
					if(studResid) {
						temp.0 <- data.frame(studResid_LGM=studres(res_LGM))
						is.singular <- "try-error"%in%class(try(temp.0[,1],silent=T))
						if(is.singular){
							warn.res2 <- "<li> Error : Cannnot calculate studentized residuals because of the (near) linear dependences in explatory variables. Please check multicollinearity problem."
						} else {
							temp <- data.frame(studResid_LGM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
							mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
							mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
							if(exists('O')) {
								O <- cbind(O,studResid_train_LGM=mer.1[,3])
							} else {
								O <- data.frame(studResid_train_LGM=mer.1[,3])
							}
						}
					}
					if (cook_distance) {
						temp.0 <- data.frame(CookDist_LGM=cooks.distance(res_LGM))
						temp <- data.frame(CookDist_LGM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
						mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
						mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
						if(exists('O')) {
							O <- cbind(O,CookDist_train_LGM=mer.1[,3])
						} else {
							O <- data.frame(CookDist_train_LGM=mer.1[,3])
						}
					}
					if (linear_pred) {
						temp.0 <- data.frame(LinearPred_LGM=res_LGM$linear.predictors)
						temp <- data.frame(LinearPred_LGM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
						mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
						mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
						if(exists('O')) {
							O <- cbind(O,LinearPred_train_LGM=mer.1[,3])
						} else {
							O <- data.frame(LinearPred_train_LGM=mer.1[,3])
						}
					}
					if (hat_value) {
						temp.0 <- data.frame(HatValue_LGM=hatvalues(res_LGM))
						temp <- data.frame(HatValue_LGM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
						mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
						mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
						if(exists('O')) {
							O <- cbind(O,HatValue_train_LGM=mer.1[,3])
						} else {
							O <- data.frame(HatValue_train_LGM=mer.1[,3])
						}
					}
					
					# test dataset
					if(exists('test.dataset') & Predict_prob_test) {
						temp.0 <- data.frame(Fitted_LGM=predict(res_LGM,newdata=test.dataset,type='response'))
						temp <- data.frame(Fitted_LGM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
						mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
						mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
						if(exists('O')) {
							O <- cbind(O,Predicted_testing_LGM=mer.1[,3])
						} else {
							O <- data.frame(Predicted_testing_LGM=mer.1[,3])
						}
					}
					if(exists('test.dataset') & Predict_CI_test){
						temp.0 <- as.data.frame(getIntervals(LGM_obj=res_LGM,newdat=test.dataset,con.level=confint.level_test))
						temp <- data.frame(temp=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
						mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
						mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),c((ncol(mer)-1),ncol(mer))]
						colnames(mer.1) <- paste0('Predicted_',round(confint.level_test*100,0),'CI_',c('Lower','Upper'),'_testing_LGM')
					
						if(exists('O')) {
							O <- cbind(O,mer.1)
						} else {
							O <- mer.1
						}
					}
					if(exists('test.dataset') & Predict_g_test) {
						temp.0 <- data.frame(Fitted_group_LGM=ifelse(predict(res_LGM,newdata=test.dataset,type='response')<Cutoff,levels(temp.dat[,dep_var])[1],levels(temp.dat[,dep_var])[2]))
						temp.0[,1] <- as.character(temp.0[,1])
						temp <- data.frame(Fitted_group_LGM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
						mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
						mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
						if(exists('O')) {
							O <- cbind(O,Predicted_group_testing_LGM=mer.1[,3],stringsAsFactors=F)
						} else {
							O <- data.frame(Predicted_group_testing_LGM=mer.1[,3],stringsAsFactors=F)
						}
					}
					
					# prediction dataset
					if(exists('pred.dataset')){
						# Probability
						temp.0 <- data.frame(Predicted_LGM=predict(res_LGM,newdata=pred.dataset,type='response'))
						temp <- data.frame(Predicted_LGM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
						mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
						mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
						if(exists('O')) {
							O <- cbind(O,Predicted_pred_LGM=mer.1[,3])
						} else {
							O <- data.frame(Predicted_pred_LGM=mer.1[,3])
						}
						if(Predict_CI_pred){
							temp.0 <- as.data.frame(getIntervals(LGM_obj=res_LGM,newdat=pred.dataset,con.level=confint.level_pred))
							temp <- data.frame(temp=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
							mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
							mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),c((ncol(mer)-1),ncol(mer))]
							colnames(mer.1) <- paste0('Predicted_',round(confint.level_pred*100,0),'CI_',c('Lower','Upper'),'_pred_LGM')
						
							if(exists('O')) {
								O <- cbind(O,mer.1)
							} else {
								O <- mer.1
							}
						}
						# group
						temp.0 <- data.frame(Predicted_group_LGM=ifelse(predict(res_LGM,newdata=pred.dataset,type='response')<Cutoff,levels(temp.dat[,dep_var])[1],levels(temp.dat[,dep_var])[2]))
						temp.0[,1] <- as.character(temp.0[,1])
						temp <- data.frame(Predicted_group_LGM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
						mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
						mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
						if(exists('O')) {
							O <- cbind(O,Predicted_group_pred_LGM=mer.1[,3],stringsAsFactors=F)
						} else {
							O <- data.frame(Predicted_group_pred_LGM=mer.1[,3],stringsAsFactors=F)
						}
					}
					
					if(Part_index){
						temp <- rep('None',nrow(original.dataset))
						temp[rownames(original.dataset)%in%rownames(dataset)] <- 'Training'
						if(exists('test.dataset')) temp[rownames(original.dataset)%in%rownames(test.dataset)] <- 'Testing'
						if(exists('pred.dataset')) temp[rownames(original.dataset)%in%rownames(pred.dataset)] <- 'Prediction'
						if(exists('O')) {
							O <- cbind(O,Partition_idx_LGM=temp,stringsAsFactors=F)
						} else {
							O <- data.frame(Partition_idx_LGM=temp,stringsAsFactors=F)
						}
					}

					if(exists('O'))	O[is.na(O)] <- ''

					### Default output
					# Data Structure
					R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file="./test.html")
					total.var <- ncol(original.dataset)
					used.var <- ifelse(is.null(vars),0,length(unique(unlist(strsplit(vars,":")))))+length(c(Part.var,Pred.var))
					none.n <- nrow(original.dataset)-nrow(dataset)-ifelse(exists("test.dataset"),nrow(test.dataset),0)-ifelse(exists("pred.dataset"),nrow(pred.dataset),0)
					# total.n <- paste0(nrow(original.dataset)," (Training: ",nrow(dataset),ifelse(exists("test.dataset"),paste0(", Testing: ",nrow(test.dataset)),""),ifelse(exists("pred.dataset"),paste0(", Prediction: ",nrow(pred.dataset)),""),ifelse(none.n!=0,paste0(", None: ",none.n),""),")")
					total.n <- paste0(nrow(original.dataset)," (Missing: ",sum(!complete.cases(original.dataset[,c(var_info,Part.var,Pred.var),drop=F])),")")
					
					DS <- matrix(c('Number of observations','Number of total variables','Number of used variables',total.n,total.var,used.var+1),ncol=2)
					R2HTML::HTML(DS,file="./test.html",align="left")
					if(exists('warn.DP4')) R2HTML::HTML(warn.DP4,file="./test.html")
					if(exists('warn.DP7')) R2HTML::HTML(warn.DP7,file="./test.html")
					if(exists('warn.DP8')) R2HTML::HTML(warn.DP8,file="./test.html")
					if(exists('warn.DP9')) R2HTML::HTML(warn.DP9,file="./test.html")
					if(exists('warn.DP2')) R2HTML::HTML(warn.DP2,file="./test.html")
					if(exists('warn.DP3')) R2HTML::HTML(warn.DP3,file="./test.html")

					# Varibale list
					R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file="./test.html")
					if(is.null(vars)) {
						Vars <- c()
					} else {
						Vars <-  unique(unlist(strsplit(vars,":")))
					}

					qual <- c(indep_cat_var[indep_cat_var%in%Vars],dep_var,Pred.var,Part.var)
					quan <- indep_numeric_var[indep_numeric_var%in%Vars]
					varlist.qual <- matrix(c('Qualitative variable',paste(qual,collapse=', ')),ncol=2,byrow=T)
					if(length(quan)!=0) {
						varlist.quan <- matrix(c('Quantitative variable',paste(quan,collapse=', ')),ncol=2,byrow=T)
					} else {
						varlist.quan <- NULL
					}
					varlist <- rbind(varlist.quan,varlist.qual)
					R2HTML::HTML(varlist,file="./test.html",align="left")
					if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file="./test.html")
					if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file="./test.html")

					# Analysis Description
					R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file="./test.html")
					Alllevel <- levels(temp.dat[,dep_var])
					Alllevel[Alllevel==dep_ref] <- paste0(dep_ref,"(baseline)")
					AD <- matrix(c('Dependent variable',dep_var,'Levels of dependent variable',paste0(Alllevel,collapse=', ')),ncol=2,byrow=T)
					if(!is.null(vars)) AD <- rbind(AD,c('Explanatory variable',paste(vars,collapse=', ')))
					AD <- rbind(AD,matrix(c('Intercept included',!noint,'Link function',link,'Variable selection',select),ncol=2,byrow=T))
					if(select) {
						direction <- switch(direct,forward="Forward selection",backward="Backward elimination",both="Stepwise regression")
						AD <- rbind(AD,matrix(c('Variable selection method',direction),ncol=2,byrow=T))
						if(!is.null(keep_var)) AD <- rbind(AD,c('Fixed variable for variable selection',paste(keep_var,collapse=', ')))
					}
					if(Valid.method=='Partition'){
						VM <- ifelse(!exists('test.dataset'),"Validation using training dataset","Validation using training/testing dataset")
						if(exists('test.dataset')){
							if(Part.method=='percent'){
								PM <- paste0("Randomly split by percent (training: ",train.perc,"%, testing: ",100-train.perc,"%)")
							} else if(Part.method=='variable'){
								PM <- paste0("Split by variable '",Part.var,"' (1: training, 2: testing)")
							}
						}
					} else {
						VM <- ifelse(Cross.method=='LOOCV',"Leave-One-Out cross validation",paste0(k,"-fold cross validation"))
					}
					AD <- rbind(AD,matrix(c('Validation method',VM),ncol=2,byrow=T))
					if(exists('PM')) AD <- rbind(AD,matrix(c('Data splitting method for validataion',PM),ncol=2,byrow=T))
					if(exists('pred.dataset')){
						AD <- rbind(AD,matrix(c('Prediction using new dataset','TRUE','Data splitting method for prediction',paste0("Split by variable '",Pred.var,"' (1: prediction, 2: training/testing)")),ncol=2,byrow=T))
					} else {
						AD <- rbind(AD,matrix(c('Prediction using new dataset','FALSE'),ncol=2,byrow=T))
					}
					R2HTML::HTML(AD,file="./test.html",align="left")
					if(exists('warn.VS')) R2HTML::HTML(warn.VS,file="./test.html")
					if(exists('warn.DP4')) R2HTML::HTML(warn.DP4,file="./test.html")
					if(exists('warn.DP7')) R2HTML::HTML(warn.DP7,file="./test.html")
					if(exists('warn.DP3')) R2HTML::HTML(warn.DP3,file="./test.html")

					#### Results
					R2HTML::HTML(R2HTML::as.title("Results of Generalized Linear Regression for Binomial Data"),HR=2,file="./test.html")
					# Coefficient Estimate
					R2HTML::HTML(R2HTML::as.title("Coefficient Estimates"),HR=3,file="./test.html")
					VS1 <- capture.output(summary(res_LGM))
					blank1 <- which(VS1=="")
					IF1 <- paste(VS1[blank1[length(blank1)]-1],collapse="")
					IF1 <- gsub("[ ]{2,}"," ",IF1)
					IF2 <- paste(VS1[blank1[length(blank1)]-3],collapse="")
					IF2 <- gsub("[ ]{2,}"," ",IF2)
					IF3 <- paste(VS1[blank1[length(blank1)]-4],collapse="")
					IF3 <- gsub("[ ]{2,}"," ",IF3)
					IF4 <- paste(VS1[blank1[length(blank1)]-5],collapse="")
					IF4 <- gsub("[ ]{2,}"," ",IF4)
					IF5 <- paste(VS1[blank1[length(blank1)]-7],collapse="")
					IF5 <- gsub("[ ]{2,}"," ",IF5)

					IF <- rbind(IF5,IF4,IF3,IF2,IF1)
					rownames(IF) <- c("","","","","")
					
					CE <- as.data.frame(summary(res_LGM)$coef)
					colnames(CE) <- c('Estimate','SE','Z-value','P-value')
					if(CI){
						tmp <- merge(CE,confint.default(res_LGM, level=confint.level),by="row.names",all=TRUE,sort=F)
						rownames(tmp) <- tmp[,1]
						CE <- tmp[,-1]
						colnames(CE)[(ncol(CE)-1):ncol(CE)] <- paste0(c('Lower','Upper'),' bound of<br>',confint.level*100,'% CI for<br>Estimate')
					}

					if(odds) {
						ORs <- exp(CE[,1,drop=F])
						colnames(ORs) <- 'exp(Estimate)'
						CE <- cbind(CE,ORs)
						CE <- CE[,c(1,ncol(CE),2:(ncol(CE)-1))]
						if(CI){
							tmp <- merge(CE,exp(confint.default(res_LGM, level=confint.level)),by="row.names",all=TRUE,sort=F)
							rownames(tmp) <- tmp[,1]
							CE <- tmp[,-1]
							colnames(CE)[(ncol(CE)-1):ncol(CE)] <- paste0(c('Lower','Upper'),' bound of<br>',confint.level*100,'% CI for<br>exp(Estimate)')
						}
					}

					if(VIF)	{
						if(is.null(vars)){
							warn.vif <- '<li> Warning : VIF is not supported for intercept-only model.'
						} else {
							VIF.1 <- rms::vif(res_LGM_1)
							if(!noint) {
							  VIF.1 <- c(NA,VIF.1)
							  names(VIF.1)[1] <- "(Intercept)"
							}
							
							tmp <- merge(CE,VIF.1,by='row.names',all=TRUE)
							rownames(tmp) <- tmp[,1]
							colnames(tmp)[ncol(tmp)] <- "VIF"
							CE <- tmp[,-1]
						}
					}
					if(VIP){
						if(!is.null(vars)){
							vip <- varImp(res_LGM)
							if(!noint) {
							  vip <- rbind(NA,vip)
							  rownames(vip)[1] <- "(Intercept)"
							}

							tmp <- merge(CE,vip,by='row.names',all=TRUE)
							rownames(tmp) <- tmp[,1]
							colnames(tmp)[ncol(tmp)] <- "VIP"
							CE <- tmp[,-1]
						} else {
							warn.msg9 <- "<li> Warning : Variable importance table is not supported for intercept only model."
						}
					}
					R2HTML::HTML(Digits(CE),file="./test.html",align="left",digits=15)
					if(exists('warn.msg9')) R2HTML::HTML(warn.msg9,file="./test.html")
					if(exists('warn.vif')) R2HTML::HTML(warn.vif,file="./test.html")
					R2HTML::HTML(IF,file="./test.html", align="left",digits=15)
				
					#### Anova table, Goodness of fit, R-squared and classfication table
					# ANOVA=TRUE;GOF=TRUE;ng=10;R2=TRUE;classtab=TRUE;cutpoint=.5
					if(ANOVA){
						R2HTML::HTML(R2HTML::as.title("Analysis-of-Deviance Table"),HR=3,file="./test.html")
						if(length(vars)==0){
							warn.msg8 <- "<li> Warning : ANOVA table is not supported for intercept only model."
							R2HTML::HTML(warn.msg8,file="./test.html")
						} else if(noint){
							warn.AT1 <- "<li> Warning : ANOVA table is not supported for the model without intercept."
							R2HTML::HTML(warn.AT1,file="./test.html")
						} else {
							R2HTML::HTML(R2HTML::as.title("Model Effect (Goodness of Fit Test)"),HR=4,file="./test.html")
							command_str<- paste0("null.fit<- try(glm(formula(gsub('~.+','~ 1',form.0)), data=temp.dat, family=binomial(",link,")))")
							eval(parse(text=command_str)) ;
							model.fit <- data.frame(anova(null.fit,res_LGM))
							model.fit <- model.fit[c(2,1),c(2,1,4,3)]
							model.fit <- cbind(model.fit,c(pchisq(model.fit[1,3],model.fit[1,4],lower.tail=F),NA))
							rownames(model.fit) <- c('Proposed model','Null model')
							colnames(model.fit) <- c('Deviance','DF(Deviacne)','Chisq','DF(Chisq)','P-value')
							R2HTML::HTML(Digits(model.fit),file="./test.html",align="left",digits=15)
							anova.msg <- "<li> Note : 'Null model' means the model including only intercept and 'Proposed model' means the model including all explanatory variables including interaction effect."
							R2HTML::HTML(anova.msg,file="./test.html")

							R2HTML::HTML(R2HTML::as.title(paste0("Variable Effect with Type ",ss," SS")),HR=4,file="./test.html")
							warn.desc <- ifelse(ss=='I',"<li> Note : In type I test, 'Null model' means the model including only intercept. Terms are added sequentially (first to last).",
									ifelse(ss=='II',"<li> Note : In type II test, 'Proposed model' means the model including all explanatory variables except interaction effect. The other rows are the testing results for each main effect after the other main effect.",
									"<li> Note : In type III test, 'Proposed model' means the model including all explanatory variables including interaction effect. The other rows are the testing results for each effect after the other effect."))
							
							if(ss=='I'){
								AT <- try(anova(res_LGM,test="Chisq"),s=T)
								rowNames <- rownames(AT)
								if(class(AT)[1]!='try-error'){
									AT <- data.frame(AT)
									AT <- AT[,c(4,3,2,1,5)]
									rownames(AT) <- rowNames; rownames(AT)[1] <- 'Null model'
									colnames(AT) <- c('Deviance','DF(Deviacne)','Chisq','DF(Chisq)','P-value')
									R2HTML::HTML(Digits(AT),file="./test.html",align="left",digits=15)
									R2HTML::HTML(warn.desc,file="./test.html")
								} else {
									warn.msg6 <- "<li> Error : Fail to fit the Model."
									R2HTML::HTML(warn.msg6,file="./test.html")
								}
							}

							if(ss %in% c('II','III')){
								# II type : ignore interaction term
								# III type : calculate SS including interaction term
								RN <- vars
								if(ss=='II' & length(grep(':',RN))>0) {
									RN <- RN[-grep(':',RN)]
									warn.AT2 <- '<li> Warning : Test for interaction effects is not provided in type II test. Use type III test for interaction effects.'
								}
								# warn.AT3 <- '<li> Note : If there is indeed no interaction, then type II is statistically more powerful than type III.'

								AT.form.full <- paste(dep_var,'~',paste(RN,collapse=' + '))
								command_str<- paste0("full.fit<- try(glm(formula(AT.form.full), data=temp.dat, family=binomial(",link,")))")

								eval(parse(text=command_str))
								options(contrasts=c("contr.sum", "contr.poly"))
								res <- try(car::Anova(full.fit,type='III',test="LR"))
								options(contrasts=c('contr.treatment','contr.poly'))

								if(class(res)[1]!='try-error'){
									dev <- deviance(full.fit)
									dev.df <- df.residual(full.fit)
									devs <- dev+res[,1]
									dev.dfs <- dev.df + res[,2]
									AT <- data.frame(c(dev,devs),c(dev.df,dev.dfs),c(NA,res[,1]),c(NA,res[,2]),c(NA,res[,3]))
									colnames(AT) <- c('Deviance','DF(Deviacne)','Chisq','DF(Chisq)','P-value')
									rownames(AT) <- c('Proposed model',RN)

									R2HTML::HTML(Digits(AT),file="./test.html",align="left",digits=15)
									R2HTML::HTML(warn.desc,file="./test.html")
									if(exists('warn.AT2')) R2HTML::HTML(warn.AT2,file="./test.html")
								} else {
									warn.msg7 <- "<li> Error : Fail to fit the Model."
									R2HTML::HTML(warn.msg7,file="./test.html")
								}
							}
						}
					}

					## Goodness of fit
					R2HTML::HTML(R2HTML::as.title("Model Fitness Measurements"), HR=3, file="./test.html") ;
					fitness_print	<- data.frame(numeric(0), numeric(0)) ;
					fitness_print	<- rbind(fitness_print, c(sum(residuals(res_LGM, type="deviance")^2), res_LGM$df.residual)) ;
					fitness_print	<- rbind(fitness_print, c(sum(residuals(res_LGM, type="pearson")^2), res_LGM$df.residual)) ;
					LL <- logLik(res_LGM)
					fitness_print	<- rbind(fitness_print, c(-2*LL[1], attr(LL,'df'))) ;
					fitness_print	<- rbind(fitness_print, c(AIC(res_LGM), NA)) ;
					fitness_print	<- rbind(fitness_print, c(BIC(res_LGM), NA)) ;
					names(fitness_print)	<- c("Value", "DF") ;
					row.names(fitness_print)	<- c("Deviance", "Pearson's chi-square", "-2*log-likelihood", "AIC", "BIC") ;
					R2HTML::HTML(Digits(fitness_print), file="./test.html", digits=15, align="left",caption="<div style='text-align:left'> <li> A model with a smaller value is better.") ;

					R2HTML::HTML(R2HTML::as.title("Pseudo R-squared Measures"),HR=4,file="./test.html")
					temp.dat <<- temp.dat
					A6 <- as.data.frame(t(pscl::pR2(res_LGM)[-c(1:2)]))
					rm(temp.dat,envir=parent.env(environment()))
					R2HTML::HTML(Digits(A6),file="./test.html",align="left",row.names=F,digits=15)

					if(GOF){
						R2HTML::HTML(R2HTML::as.title("Goodness of Fit Test"),HR=3,file="./test.html")
						if(length(vars)!=0){
							## log-liklehood ratio
							R2HTML::HTML(R2HTML::as.title("Analysis-of-Deviance Table (Likelihood Ratio Test)"),HR=4,file="./test.html")
							if(!noint){
								command_str<- paste0("null.fit<- try(glm(formula(gsub('~.+','~ 1',form.0)), data=temp.dat, family=binomial(",link,")))")
								eval(parse(text=command_str)) ;
								model.fit <- data.frame(anova(null.fit,res_LGM))
								model.fit <- model.fit[c(2,1),c(2,1,4,3)]
								model.fit <- cbind(model.fit,c(pchisq(model.fit[1,3],model.fit[1,4],lower.tail=F),NA))
								rownames(model.fit) <- c('Proposed model','Null model')
								colnames(model.fit) <- c('Deviance','DF(Deviacne)','Chisq','DF(Chisq)','P-value')
								R2HTML::HTML(Digits(model.fit),file="./test.html",align="left",digits=15)
								anova.msg <- "<li> Note : 'Null model' means the model including only intercept and 'Proposed model' means the model including all explanatory variables including interaction effect."
								R2HTML::HTML(anova.msg,file="./test.html")
							} else {
								R2HTML::HTML('<li> Warning : Likelihood Ratio Test is not supported for the model without intercept.',file="./test.html")
							}
							##Hosmer-Lemeshow Test
							GOF2 <- ResourceSelection::hoslem.test(res_LGM$y,fitted(res_LGM))
							A4 <- data.frame(GOF2$statistic,GOF2$parameter,GOF2$p.value)
							colnames(A4) <- c("X-squared","df","p-value")
							rownames(A4) <- ""
							R2HTML::HTML(R2HTML::as.title("Hosmer-Lemeshow test"),HR=4,file="./test.html")
							R2HTML::HTML(Digits(A4),file="./test.html", align="left",row.names=F,digits=15)
						} else {
							warn.msg12 <- "<li> Warning : Goodness of fit test is not supported for intercept only model."
							R2HTML::HTML(warn.msg12,file="./test.html")
						}
					}

					# classification table
					# classtab=TRUE;cutpoint=.5
					if(classtab){
						R2HTML::HTML(R2HTML::as.title(paste0("Classification Table - cut point : ",Cutoff)),HR=3,file="./test.html")		
						R2HTML::HTML(R2HTML::as.title("Training dataset"),HR=4,file="./test.html")		
						Observed <- temp.dat[,dep_var]
						level <- levels(temp.dat[,dep_var])
						Predicted <- as.factor(ifelse(fitted(res_LGM)<Cutoff,level[1],level[2]))
						Predicted <- factor(Predicted,levels=level)
						A7 <- table(Observed,Predicted)
						rownames(A7) <- level
						Total <- colSums(A7)
						A7 <- rbind(A7,Total)
						Total <- rowSums(A7)
						A7 <- cbind(A7,Total)						
						colnames(A7) <- paste("Fitted_",colnames(A7),sep="")
						rownames(A7) <- paste("Observed_",rownames(A7),sep="")
						R2HTML::HTML(A7,file="./test.html",align="left",digits=4)
						
						R2HTML::HTML(R2HTML::as.title("Diagnostic Accuracy Measure"),HR=5,file="./test.html")		
						aucs <- unlist(ROCR::performance(ROCR::prediction(res_LGM$fitted,res_LGM$y),'auc')@y.values)
						AM.1 <- data.frame(Measure=c('True positive (TP)','True negative (TN)','False positive (FP)','False negative (FN)'),
											Description=c(rep('',4)),
											Value=c(A7[2,2],A7[1,1],A7[1,2],A7[2,1]))
						AM.2 <- data.frame(Measure=c('AUC','Accuracy','Balanced accuracy','Sensitivity (Sens)','Specificity (Spec)','False positive rate (FPR)','False negative rate (FNR)'),
											Description=c('','(TP+TN)/(TP+TN+FP+FN)','(TP/(TP+FN)+TN/(TN+FP))/2','TP/(TP+FN)','TN/(TN+FP)','FP/(TN+FP)','FN/(TP+FN)'),
											Value=c(aucs,(AM.1[1,3]+AM.1[2,3])/sum(AM.1[,3]),(AM.1[1,3]/(AM.1[1,3]+AM.1[4,3])+AM.1[2,3]/(AM.1[2,3]+AM.1[3,3]))/2,AM.1[1,3]/(AM.1[1,3]+AM.1[4,3]),AM.1[2,3]/(AM.1[2,3]+AM.1[3,3]),AM.1[3,3]/(AM.1[2,3]+AM.1[3,3]),AM.1[4,3]/(AM.1[1,3]+AM.1[4,3])))
						AM.3 <- data.frame(Measure=c('Positive predictive value (PPV)','Negative predictive value (NPV)','False discovery rate (FDR)','False omission rate (FOR)'),
											Description=c('TP/(TP+FP)','TN/(TN+FN)','FP/(TP+FP)','FN/(TN+FN)'),
											Value=c(AM.1[1,3]/(AM.1[1,3]+AM.1[3,3]),AM.1[2,3]/(AM.1[2,3]+AM.1[4,3]),AM.1[3,3]/(AM.1[1,3]+AM.1[3,3]),AM.1[4,3]/(AM.1[2,3]+AM.1[4,3])))
						AM.4 <- data.frame(Measure=c('Positive likelihood ratio (LR+)','Negative lieklihood ratio (LR-)'),
											Description=c('Sens/FPR','FNR/Spec'),
											Value=c(AM.2[4,3]/AM.2[6,3],AM.2[7,3]/AM.2[5,3]))
						AM.5 <- data.frame(Measure=c('Diagnostic odds ratio (DOR)','F1 score'),
											Description=c('LR+/LR-','2/(1/Sens+1/PPV)'),
											Value=c(AM.4[1,3]/AM.4[2,3],2/(1/AM.2[4,3]+1/AM.3[1,3])))
						AM <- rbind(AM.1,AM.2,AM.3,AM.4,AM.5)
						AM[,3] <- Digits(AM[,3])
						R2HTML::HTML(AM,file="./test.html",align="left",digits=15,row.names=F)

						if(exists('test.dataset')){
							R2HTML::HTML(R2HTML::as.title("Testing dataset"),HR=4,file="./test.html")		
							Observed <- test.dataset[,dep_var]
							level <- levels(temp.dat[,dep_var])

							pp <- predict(res_LGM,newdata=test.dataset,type='response')
							Predicted <- as.factor(ifelse(pp<Cutoff,level[1],level[2]))
							Predicted <- factor(Predicted,levels=level)
							A7 <- table(Observed,Predicted)
							rownames(A7) <- level
							Total <- colSums(A7)
							A7 <- rbind(A7,Total)
							Total <- rowSums(A7)
							A7 <- cbind(A7,Total)						
							colnames(A7) <- paste("Predicted_",colnames(A7),sep="")
							rownames(A7) <- paste("Observed_",rownames(A7),sep="")
							R2HTML::HTML(A7,file="./test.html",align="left",digits=4)
							
							R2HTML::HTML(R2HTML::as.title("Measures of Diagnostic Accuracy"),HR=5,file="./test.html")		
							test.res <- na.omit(cbind(predict(res_LGM,newdata=test.dataset,type='response'),test.dataset[,dep_var]))
							aucs <- unlist(ROCR::performance(ROCR::prediction(test.res[,1],test.res[,2]),'auc')@y.values)
							AM.1 <- data.frame(Measure=c('True positive (TP)','True negative (TN)','False positive (FP)','False negative (FN)'),
												Description=c(rep('',4)),
												Value=c(A7[2,2],A7[1,1],A7[1,2],A7[2,1]))
							AM.2 <- data.frame(Measure=c('AUC','Accuracy','Balanced accuracy','Sensitivity (Sens)','Specificity (Spec)','False positive rate (FPR)','False negative rate (FNR)'),
												Description=c('','(TP+TN)/(TP+TN+FP+FN)','(TP/(TP+FN)+TN/(TN+FP))/2','TP/(TP+FN)','TN/(TN+FP)','FP/(TN+FP)','FN/(TP+FN)'),
												Value=c(aucs,(AM.1[1,3]+AM.1[2,3])/sum(AM.1[,3]),(AM.1[1,3]/(AM.1[1,3]+AM.1[4,3])+AM.1[2,3]/(AM.1[2,3]+AM.1[3,3]))/2,AM.1[1,3]/(AM.1[1,3]+AM.1[4,3]),AM.1[2,3]/(AM.1[2,3]+AM.1[3,3]),AM.1[3,3]/(AM.1[2,3]+AM.1[3,3]),AM.1[4,3]/(AM.1[1,3]+AM.1[4,3])))
							AM.3 <- data.frame(Measure=c('Positive predictive value (PPV)','Negative predictive value (NPV)','False discovery rate (FDR)','False omission rate (FOR)'),
												Description=c('TP/(TP+FP)','TN/(TN+FN)','FP/(TP+FP)','FN/(TN+FN)'),
												Value=c(AM.1[1,3]/(AM.1[1,3]+AM.1[3,3]),AM.1[2,3]/(AM.1[2,3]+AM.1[4,3]),AM.1[3,3]/(AM.1[1,3]+AM.1[3,3]),AM.1[4,3]/(AM.1[2,3]+AM.1[4,3])))
							AM.4 <- data.frame(Measure=c('Positive likelihood ratio (LR+)','Negative lieklihood ratio (LR-)'),
												Description=c('Sens/FPR','FNR/Spec'),
												Value=c(AM.2[4,3]/AM.2[6,3],AM.2[7,3]/AM.2[5,3]))
							AM.5 <- data.frame(Measure=c('Diagnostic odds ratio (DOR)','F1 score'),
												Description=c('LR+/LR-','2/(1/Sens+1/PPV)'),
												Value=c(AM.4[1,3]/AM.4[2,3],2/(1/AM.2[4,3]+1/AM.3[1,3])))
							AM <- rbind(AM.1,AM.2,AM.3,AM.4,AM.5)
							AM[,3] <- Digits(AM[,3])
							R2HTML::HTML(AM,file="./test.html",align="left",digits=15,row.names=F)
						}
					}

					#### Figures - jhan : ggplot style로 수정
					if(ROCcurve) {
						R2HTML::HTML(R2HTML::as.title("ROC Curve"),HR=3,file="./test.html",append=TRUE)
						R2HTML::HTML(R2HTML::as.title("Training dataset"),HR=4,file="./test.html",append=TRUE)
						REx_ANA_PLOT()
						train.aucs <- unlist(ROCR::performance(ROCR::prediction(res_LGM$fitted,res_LGM$y),'auc')@y.values)
						train.ROCdat <- fortify(ROCR::performance(ROCR::prediction(res_LGM$fitted,res_LGM$y),'tpr','fpr'))
						print(ggplot(data=train.ROCdat,
								aes(x=False.positive.rate, y=True.positive.rate)) + geom_line(size=1) +
							geom_abline(intercept=0, slope=1, color="grey") +   annotate("text", x = 0.9, y = 0.05, label = paste0("AUC : ", Digits(train.aucs))) +
							labs(x="False Positive Rate", y="True Positive Rate") +
							theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()))
						REx_ANA_PLOT_OFF("")
						
						if(exists('test.dataset')){
							R2HTML::HTML(R2HTML::as.title("Testing dataset"),HR=4,file="./test.html",append=TRUE)
							REx_ANA_PLOT()
							test.res <- na.omit(cbind(predict(res_LGM,newdata=test.dataset,type='response'),test.dataset[,dep_var]))
							test.aucs <- unlist(ROCR::performance(ROCR::prediction(test.res[,1],test.res[,2]),'auc')@y.values)
							test.ROCdat <- fortify(ROCR::performance(ROCR::prediction(test.res[,1],test.res[,2]),'tpr','fpr'))
							print(ggplot(data=test.ROCdat,
									aes(x=False.positive.rate, y=True.positive.rate)) + geom_line(size=1) +
								geom_abline(intercept=0, slope=1, color="grey") +   annotate("text", x = 0.9, y = 0.05, label = paste0("AUC : ", Digits(test.aucs))) +
								labs(x="False Positive Rate", y="True Positive Rate") +
								theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()))
							REx_ANA_PLOT_OFF("")
							
							## Combined
							R2HTML::HTML(R2HTML::as.title("Training & Testing dataset"),HR=4,file="./test.html",append=TRUE)
							REx_ANA_PLOT(h=500,w=575)
							combined.dat <- rbind(train.ROCdat,test.ROCdat)
							combined.dat$dataset <- factor(c(rep('Training',nrow(train.ROCdat)),rep('Testing',nrow(test.ROCdat))),levels=c('Training','Testing'))
							print(ggplot(data=combined.dat,
									aes(x=False.positive.rate, y=True.positive.rate,group=dataset,col=dataset)) + geom_line(size=1) +
								geom_abline(intercept=0, slope=1, color="grey") +   annotate("text", x = 0.75, y = 0.05, label = paste0("AUC : ", Digits(train.aucs),"(Training), ",Digits(test.aucs),"(Testing)")) +
								labs(x="False Positive Rate", y="True Positive Rate") +
								theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()))
							REx_ANA_PLOT_OFF("")
						}
					}

					#### validataion ####
					if(Valid.method=='Partition'){
						R2HTML::HTML(R2HTML::as.title(VM),HR=3,file="./test.html")
					} else {
						if(Cross.method=='LOOCV'){
							R2HTML::HTML(R2HTML::as.title("Leave-One-Out Cross Validation"),HR=3,file="./test.html")
						} else {
							if(k>nrow(temp.dat)) {
								warn.valid <- paste0('<li> Warning : The number of folds cannot be greater than the number of non-missing observations. (# of non-missing observations: ',nrow(temp.dat),', # of folds specified by user: ',k,', readjusted # of folds: ',nrow(temp.dat),')')  
								k <- nrow(temp.dat)
								R2HTML::HTML(R2HTML::as.title("Leave-One-Out Cross Validation"),HR=3,file="./test.html")
							} else {
								R2HTML::HTML(R2HTML::as.title(paste0(k,"-Fold Cross Validation")),HR=3,file="./test.html")
							}
						}
					}
					if(length(vars)!=0){
						if(Valid.method=='Partition'){
							level <- levels(temp.dat[,dep_var])
							if(exists('test.dataset')){
								test.p <- predict(res_LGM,newdata=test.dataset,type='response')
								test.Predicted <- as.factor(ifelse(test.p<Cutoff,level[1],level[2]))
								test.Observed <- test.dataset[,dep_var]
								test.res <- test.Predicted==test.Observed
								test.trueness <- sum(test.res,na.rm=T)
								test.n <- sum(!is.na(test.res))
								test.Accuracy <- test.trueness/test.n
								test.Kappa <- psych::cohen.kappa(x=cbind(test.Predicted,test.Observed))$kappa
								test.AUC <- unlist(ROCR::performance(ROCR::prediction(test.p,test.Observed),'auc')@y.values)
							} else {
								train.Perc <- 100
							}
							train.Observed <- temp.dat[,dep_var]
							train.Predicted <- factor(ifelse(fitted(res_LGM)<Cutoff,level[1],level[2]),levels=level)
							train.res <- train.Predicted==train.Observed
							train.trueness <- sum(train.res,na.rm=T)
							train.n <- sum(!is.na(train.res))
							train.Accuracy <- train.trueness/train.n
							train.Perc <- ifelse(exists('test.dataset'),train.n/(train.n+test.n)*100,100)
							train.Kappa <- psych::cohen.kappa(x=cbind(train.Predicted,train.Observed))$kappa
							train.AUC <- unlist(ROCR::performance(ROCR::prediction(res_LGM$fitted,res_LGM$y),'auc')@y.values)

							vali <- data.frame(N.observed=train.n,Percent=train.Perc,N.trueness=train.trueness,Accuracy=train.Accuracy,Kappa=train.Kappa,AUC=train.AUC)
							rownames(vali) <- 'Training'
							if(exists('test.dataset')){
								vali <- rbind(vali,data.frame(N.observed=test.n,Percent=(100-train.Perc),N.trueness=test.trueness,Accuracy=test.Accuracy,Kappa=test.Kappa,AUC=test.AUC))
								rownames(vali)[2] <- 'Testing'
							}
							colnames(vali)[1] <- 'N.non-missing<br>observations'
							R2HTML::HTML(Digits(vali),file="./test.html",align="left",digits=15)
							R2HTML::HTML("<li> Note : Be careful of overfitting in interpreting the result of the training dataset",file="./test.html")
							if(exists('warn.DP4')) R2HTML::HTML(warn.DP4,file="./test.html")
							if(exists('warn.DP7')) R2HTML::HTML(warn.DP7,file="./test.html")
							if(exists('warn.DP8')) R2HTML::HTML(warn.DP8,file="./test.html")
							R2HTML::HTML("<li> Note : For more measures of diagnostic accuracy, please use classification table option.",file="./test.html")
						} else {
							levs <- levels(temp.dat[,dep_var])
							levels(temp.dat[,dep_var]) <- make.names(levs)
							ctrl <- trainControl(method = ifelse(k==nrow(temp.dat)|Cross.method=='LOOCV',"LOOCV","cv"), number = k, savePredictions = TRUE,returnResamp='all',classProbs=TRUE)
							newdat <- temp.dat[,c(dep_var,indep_numeric_var,indep_cat_var),drop=F]
							command_str	<- paste0("mod_fit <- train(res_LGM$formula,  data=newdat, method='glm',family=binomial(",link,"),trControl = ctrl, tuneLength = 5)")
							eval(parse(text=command_str)) ;
													
							cv.res <- mod_fit$pred
							cv.pred <- ifelse(cv.res[,make.names(levs)[levs!=dep_ref]]>Cutoff,make.names(levs)[levs!=dep_ref],make.names(levs)[levs==dep_ref])
							cv.obs <- as.character(cv.res[,'obs'])
							if(k==nrow(temp.dat)|Cross.method=='LOOCV'){
								Accuracy <- sum(cv.pred==cv.obs)/length(cv.pred)
								Kappa <- psych::cohen.kappa(x=cbind(cv.pred,cv.obs))$kappa
								AUC <- unlist(ROCR::performance(ROCR::prediction(cv.res[,make.names(levs)[levs!=dep_ref]],cv.obs),'auc')@y.values)
								res.acc.II <- data.frame(Accuracy,Kappa,AUC)
							} else {
								get.acc <- function(fold){
									# fold : fold name
									idx <- which(cv.res[,'Resample']==fold)
									Accuracy <- sum(cv.pred[idx]==cv.obs[idx])/length(idx)
									CT <- table(data.frame(Obs=factor(cv.obs[idx],levels=levs),Pred=factor(cv.pred[idx],levels=levs)))
									Kappa <- psych::cohen.kappa(x=CT)$kappa
									AUC <- unlist(ROCR::performance(ROCR::prediction(cv.res[idx,make.names(levs)[levs!=dep_ref]],cv.obs[idx]),'auc')@y.values)
									return(data.frame(Accuracy,Kappa,AUC,fold))
								}
								res.acc <- lapply(unique(cv.res[,'Resample']),get.acc)
								res.acc.I <- do.call(rbind,res.acc)
								res.acc.II <- data.frame(Accuracy = mean(res.acc.I[,'Accuracy']),AccuracySD = sd(res.acc.I[,'Accuracy']),Kappa = mean(res.acc.I[,'Kappa']),KappaSD = sd(res.acc.I[,'Kappa']),AUC=mean(res.acc.I[,'AUC']),AUCSD=sd(res.acc.I[,'AUC']))
								colnames(res.acc.II)[c(2,4,6)] <- c('SD(Accuracy)','SD(Kappa)','SD(AUC)')
							}
							levels(temp.dat[,dep_var]) <- levs
							R2HTML::HTML(Digits(res.acc.II),file="./test.html",align="left",digits=15,row.names=F)
							if(exists('warn.valid')) R2HTML::HTML(warn.valid,file="./test.html")
						} 
					} else {
						warn.msg10 <- "<li> Warning : Validation is not supported for intercept only model."
						R2HTML::HTML(warn.msg10,file="./test.html")
					}

					## Variable selection
					if(select) stepAIC.wj(object=res_LGM_1,temp.dat,dep_var,type='binom',noint,direct,keep_var,hr=3,vars,link=link,odds=odds,CI=CI,confint.level=confint.level,VIP=VIP,VIF=VIF)
					if(exists('warn.res1')|exists('warn.res2')){
						R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file="./test.html")
						if(exists('warn.res1')) R2HTML::HTML(warn.res1,file="./test.html")
						if(exists('warn.res2')) R2HTML::HTML(warn.res2,file="./test.html")
					}
				}
			}
			# Used R packages
			R2HTML::HTML(R2HTML::as.title("Used R Packages"),HR=2,file="./test.html")
			pkg.list <- list(list("Generalized linear regression for binomial data","glm","stats"),
					 list("Confidence interval for regression coefficients","confint","stats"),
					 list("Variance inflation factor (VIF)","vif","rms"),
					 list("ANOVA table",c("anova","Anova"),c("stats","car")),
					 list("Model fitness measurements","residuals, logLik, AIC, BIC","stats"),
					 list("Pseudo R-squared measures","pR2","pscl"),
					 list("Goodness of fit test",c("anova","hoslem.test"),c("stats","ResourceSelection")),
					 list("ROC curve","prediction, performance","ROCR"),
					 list("Cross validation","trainControl, train","caret"),
					 list("Predicted value","predict","stats"),
					 list("Confidence & Prediction interval for fitted & predicted value","predict","stats"),
					 list("Unstandardized residual","residuals","stats"),
					 list("Standardized residual","stdres","MASS"),
					 list("Studentized residual","studres","MASS"),
					 list("Cook's distance","cooks.distance","stats"),
					 list("Linear predictor","predict","stats"),
					 list("Diagonals of hat matrix","hatvalues","stats"))
			R2HTML::HTML(used.pkg(pkg.list), file="./test.html")

			# Analysis end time
			R2HTML::HTMLhr(file="./test.html")
			R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". Rex : Generalized Linear Regression for Binomial Data",sep=""),file="./test.html")
			R2HTML::HTMLhr(file="./test.html")
		})
		if(exists('O')){
			if (length(Dep_var) > 1) stop("Multiple dependent variables with outputs are now allowed")
			FIN.RESULT <- list(html=html.output,Output=O)
		} else {
			FIN.RESULT[which(Dep_var%in%dep_var)] <- paste(html.output,collapse="\n")
		}
	}

	return(FIN.RESULT)
}
