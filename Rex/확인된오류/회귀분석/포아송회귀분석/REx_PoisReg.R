# 포아송회귀분석
REx_PoisReg <- function(dataset, res_var, quan_var=NULL, qual_var=NULL, vars=NULL, offset=NULL, link=c("log", "identity", "sqrt"), noint=FALSE, over_disp=FALSE, 
			Valid.method=c('Partition','Cross'),Part.method=c('all','percent','variable'),train.perc=70,Part.var=NULL, Cross.method=c('LOOCV','KFOLD'),k=10,Pred.var=NULL,
			CI=TRUE,confint.level=0.95,exp_estim=FALSE,VIF=FALSE,ANOVA=TRUE,ss_type="III",GOF=FALSE,VIP=FALSE,ODtest=FALSE,
			Predict_train=FALSE,Predict_CI_train=FALSE,confint.level_train=0.95,Resid=FALSE,stdResid=FALSE,studResid=FALSE,cook_distance=FALSE,linear_pred=FALSE,hat_value=FALSE,
			Predict_test=FALSE,Predict_CI_test=FALSE,confint.level_test=0.95,
			Predict_CI_pred=FALSE,confint.level_pred=0.95,Part_index=FALSE,
			Select=FALSE,direct=c("forward","backward","both"),keep_var=NULL){
			
	####### Packages #######
	load.pkg(c("R2HTML", "car", "AICcmodavg", "rms","caret","MASS"))

	###################
	html.output <- capture.output({
		R2HTML::HTML(R2HTML::as.title("Poisson Regression"), HR=1, file="./test.html", append=FALSE) ;

		## Warnings
		# Dependent variable type
		if(!is.numeric(dataset[,res_var])) {
			warn.msg1 <- '<li> Error : Dependent variable should be non-negative integer. Analysis has been stopped.'
		} else if (any(dataset[,res_var] < 0,na.rm=T) | !all(check.integer(dataset[,res_var]),na.rm=T)) {
			warn.msg1 <- "<li> Error : Dependent variable should be non-negative integer. Analysis has been stopped."
		}

		# explanatory variable type
		if(is.null(vars)) {
			Vars <- c()
		} else {
			Vars <-  unique(unlist(strsplit(vars,":")))
		}
		qual_var <- qual_var[qual_var%in%Vars]; 
		if(length(qual_var)==0) qual_var <- NULL
		quan_var <- quan_var[quan_var%in%Vars]
		if(length(quan_var)==0) quan_var <- NULL

		if(!is.null(quan_var)) {
			is.nom <- sapply(quan_var, function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg2	<- paste0("<li> Warning : The type of variable '", paste(quan_var[is.nom],collapse=', '), "' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")
				ec		<- quan_var[is.nom]
				for(jj in ec) vars <- vars[-grep(jj,vars)]
				if(length(vars)==0)	vars <- NULL
			}
		}
		if(!is.null(qual_var)) {
			is.num <- sapply(qual_var, function(i) is.numeric(dataset[,i]))
			if(any(is.num)) warn.msg3 <- paste0("<li> Warning : The type of variable '", paste(qual_var[is.num],collapse=', '), "' is numeric but selected as the qualitative variable. It was coerced into character.")
			for(i in qual_var) dataset[,i] <- factor(dataset[,i])
		}
		if(!is.null(offset)){
			is.nom <- sapply(offset, function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.offset1	<- paste0("<li> Warning : The type of variable '", paste(offset[is.nom],collapse=', '), "' is not numeric but selected as offset variable. It was excluded from the analysis.")
				offset <- offset[!is.nom]
				if(length(offset)==0)	offset <- NULL
			}

			if(!is.null(offset)){
				is.nonposi <- sapply(offset, function(i) any(dataset[,i]<=0,na.rm=T))
				if(any(is.nonposi)){
					warn.offset2	<- paste0("<li> Warning : The variable '", paste(offset[is.nonposi],collapse=', '), "' contain non-positive observations (Offset variable should be positive). It was excluded from the analysis.")
					offset <- offset[!is.nonposi]
					if(length(offset)==0)	offset <- NULL
				}
			}
		}

		# no intercept & no explanatory variable : error
		if(noint&is.null(vars)) warn.msg4 <- '<li> Error : With no intercept, at least 1 independent variable should be selected. Analysis has been stopped.'

		# Warning in VS
		if(Select) {
			if(is.null(vars)) {
				Select <- FALSE
				warn.VS1 <- '<li> Warning : Variable selection is not supported for intercept-only model.'
			}
			if(Select & over_disp){
				Select <- FALSE
				warn.VS2 <- '<li> Warning : Variable selection is not supported for overdispersion model.'
			}
			keep_var <- keep_var[keep_var%in%vars]
			if(length(keep_var)==0) keep_var <- NULL
		}
		getOutput <- FALSE

		## Data processing
		original.dataset <- dataset
		
		# predict variable
		if(!is.null(Pred.var)){
			pred.level <- na.omit(unique(dataset[,Pred.var]))
			# No training & test dataset
			if(!2%in%pred.level) {
				warn.DP1 <- paste0("<li> Error : '2' is not observed in the split variable for prediction, '",Pred.var,"'. That is, no observations are assigned to the training & test dataset. Analysis has been stopped.")
			} else {
				if(!1%in%pred.level) {	
					# Only 2 is obaserved
					warn.DP3 <- paste0("<li> Warning : '1' is not observed in the split variable for prediction, '",Pred.var,"'. That is, no observations are assigned to the prediction dataset. Prediction is not supported in this analysis.")
					Predict_CI_pred <- Predict_PI_pred <- FALSE
				} else if(!all(pred.level%in%c(1,2))) {
					# Value other than 1 and 2
					warn.DP2 <- paste0("<li> Warning : Values other than '1' and '2' are observed in the split variable for prediction, '",Pred.var,"'. Observations with the value other than '1' and '2' are not included the analysis.")
					pred.dataset <- dataset[dataset[,Pred.var]==1,]	# prediction dataset
					dataset <- dataset[dataset[,Pred.var]==2,]	# training & test dataset
				}
			}
		}
		
		# training & test
		if(Valid.method=='Partition'){
			if(Part.method=='percent'){
				n.train <- round(nrow(dataset)*train.perc/100)
				
				if(train.perc==100 | n.train==nrow(dataset)){
					warn.DP4 <- "<li> Warning : No observations are assigned to the test dataset due to too high percent for the training dataset. Validation using test dataset is not supported in this analysis."
					Predict_test <- Predict_CI_test <- Predict_PI_test <- FALSE
				} else if(n.train==0){
					warn.DP5 <- "<li> Error : No observations are assigned to the training dataset due to too low percent for the training dataset. To secure sufficient number of observations for the training dataset, increase the percent. Analysis has been stopped."
				} else {
					train <- sample(seq(nrow(dataset)),n.train)

					# test set
					test.dataset <- dataset[-train,,drop=F]
					# training set
					dataset <- dataset[train,,drop=F]
					
					if(nrow(original.dataset)!=(nrow(dataset)+nrow(test.dataset)))	warn.DP9 <- "<li> Warning : Observations with missing dependent variable were not assigned to training dataset or test dataset."
					
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
					warn.DP7 <- paste0("<li> Warning : '2' is not observed in the split variable for validation, '",Part.var,"'. That is, no observations are assigned to the test dataset. Validation using test dataset is not supported in this analysis.")
					Predict_prob_tet <- Predict_g_test <- FALSE
				} else {
					# Value other than 1 and 2
					if(!all(part.level%in%c(1,2))) {
						warn.DP8 <- paste0("<li> Warning : Values other than '1' and '2' are observed in the split variable for validation, '",Part.var,"'. Observations with the value other than '1' and '2' are not included the analysis.")
					}
					test.dataset <- dataset[dataset[,Part.var]==2,]	# test dataset
					dataset <- dataset[dataset[,Part.var]==1,]	# training dataset
				}
			}
		}
		#### Then, original.dataset : original dataset, dataset : training dataset, test.dataset : test dataset, pred.dataset : prediction dataset.
		
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
			raw.dat <- dataset
			temp.dat <- dataset[complete.cases(dataset[,c(Vars,res_var,offset),drop=F]),,drop=F]

			if(noint){
				form.0 <- paste(res_var,'~',ifelse(is.null(offset),paste(paste(vars,collapse=' + '),'-1'),paste(paste(vars,collapse=' + '),paste0("+ offset(", paste(offset, collapse=" + "), ")"),'-1')))
			} else {
				form.0 <- paste(res_var,'~',paste(ifelse(is.null(vars),ifelse(is.null(offset),1,paste0("offset(", paste(offset, collapse=" + "), ")")),ifelse(is.null(offset),paste(vars,collapse=' + '),paste(paste(vars,collapse=' + '),paste0("+ offset(", paste(offset, collapse=" + "), ")"))))))
			}
			form.1 <- as.formula(form.0)

			link_actual	<- link
			dist	<- ifelse(over_disp, quasipoisson, poisson) ;

			fit_glm	<- res_Pois <- try(glm(form.1, data=temp.dat, family=dist(link=link_actual)))

			if('try-error'%in%class(fit_glm)){
				R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file="./test.html")
				R2HTML::HTML("<li> Warning : Poisson regression was failed to fit. Analysis has been stopped.",file="./test.html")
			} else if(!fit_glm$conv) {
				R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file="./test.html")
				warn.msg.conv <- "<li> Warning : Poisson regression did not converge. Analysis has been stopped."
				R2HTML::HTML(warn.msg.conv,file="./test.html")
			} else {
				## confidence interval
				getIntervals <- function(Pois_obj,newdat,con.level){
					# CI for linear predictor
					linpred <- predict(Pois_obj,newdata=newdat,se.fit=T)
					lower <- linpred$fit-qnorm((1-con.level)/2,lower.tail=F)*linpred$se.fit
					upper <- linpred$fit+qnorm((1-con.level)/2,lower.tail=F)*linpred$se.fit

					# Get CI for fitted value by transforming the CI for linear predictor
					LinkFt <- function(x,link){
						tf <- switch(link, log=exp(x), identity=x, sqrt=x^2)
						return(tf)
					}
					t.lower <- LinkFt(lower,link); t.lower[t.lower<0] <- 0
					t.upper <- LinkFt(upper,link)
					t.CI <- cbind(t.lower,t.upper)
					return(t.CI)
				}
		
				## 다음의 값들은 엑셀 시트에 저장
				# Training dataset
				if(Predict_train) {
					temp.0 <- data.frame(Fitted_Pois=predict(res_Pois,newdata=temp.dat,type='response'))
					temp <- data.frame(Fitted_Pois=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
					mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
					mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
					O <- data.frame(Fitted_train_Pois=mer.1[,3])
				}
				if(Predict_CI_train){
					temp.0 <- as.data.frame(getIntervals(Pois_obj=res_Pois,newdat=temp.dat,con.level=confint.level_train))
					temp <- data.frame(temp=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
					mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
					mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),c((ncol(mer)-1),ncol(mer))]
					colnames(mer.1) <- paste0('Fitted_',round(confint.level_train*100,0),'CI_',c('Lower','Upper'),'_train_Pois')
				
					if(exists('O')) {
						O <- cbind(O,mer.1)
					} else {
						O <- mer.1
					}
				}
				if(Resid) {
					temp.0 <- data.frame(Resid_Pois=resid(res_Pois))
					temp <- data.frame(Resid_Pois=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
					mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
					mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
					if(exists('O')) {
						O <- cbind(O,unstdResid_train_Pois=mer.1[,3])
					} else {
						O <- data.frame(unstdResid_train_Pois=mer.1[,3])
					}
				}
				if(stdResid) {
					temp.0 <- data.frame(stdResid_Pois=stdres(res_Pois))
					is.singular <- "try-error"%in%class(try(temp.0[,1],silent=T))
					if(is.singular){
						warn.res1 <- "<li> Error : Cannnot calculate standardized residuals because of the (near) linear dependences in explatory variables. Please check multicollinearity problem."
					} else {
						temp <- data.frame(stdResid_Pois=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
						mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
						mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
						if(exists('O')) {
							O <- cbind(O,stdResid_train_Pois=mer.1[,3])
						} else {
							O <- data.frame(stdResid_train_Pois=mer.1[,3])
						}
					}
				}
				if(studResid) {
					temp.0 <- data.frame(studResid_Pois=studres(res_Pois))
					is.singular <- "try-error"%in%class(try(temp.0[,1],silent=T))
					if(is.singular){
						warn.res2 <- "<li> Error : Cannnot calculate studentized residuals because of the (near) linear dependences in explatory variables. Please check multicollinearity problem."
					} else {
						temp <- data.frame(studResid_Pois=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
						mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
						mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
						if(exists('O')) {
							O <- cbind(O,studResid_train_Pois=mer.1[,3])
						} else {
							O <- data.frame(studResid_train_Pois=mer.1[,3])
						}
					}
				}
				if (cook_distance) {
					temp.0 <- data.frame(CookDist_Pois=cooks.distance(res_Pois))
					temp <- data.frame(CookDist_Pois=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
					mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
					mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
					if(exists('O')) {
						O <- cbind(O,CookDist_train_Pois=mer.1[,3])
					} else {
						O <- data.frame(CookDist_train_Pois=mer.1[,3])
					}
				}
				if (linear_pred) {
					temp.0 <- data.frame(LinearPred_Pois=res_Pois$linear.predictors)
					temp <- data.frame(LinearPred_Pois=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
					mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
					mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
					if(exists('O')) {
						O <- cbind(O,LinearPred_train_Pois=mer.1[,3])
					} else {
						O <- data.frame(LinearPred_train_Pois=mer.1[,3])
					}
				}
				if (hat_value) {
					temp.0 <- data.frame(HatValue_Pois=hatvalues(res_Pois))
					temp <- data.frame(HatValue_Pois=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
					mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
					mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
					if(exists('O')) {
						O <- cbind(O,HatValue_train_Pois=mer.1[,3])
					} else {
						O <- data.frame(HatValue_train_Pois=mer.1[,3])
					}
				}
				
				# test dataset
				if(exists('test.dataset') & Predict_test) {
					temp.0 <- data.frame(Fitted_Pois=predict(res_Pois,newdata=test.dataset,type='response'))
					temp <- data.frame(Fitted_Pois=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
					mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
					mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
					if(exists('O')) {
						O <- cbind(O,Predicted_test_Pois=mer.1[,3])
					} else {
						O <- data.frame(Predicted_test_Pois=mer.1[,3])
					}
				}
				if(exists('test.dataset') & Predict_CI_test){
					temp.0 <- as.data.frame(getIntervals(Pois_obj=res_Pois,newdat=test.dataset,con.level=confint.level_test))
					temp <- data.frame(temp=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
					mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
					mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),c((ncol(mer)-1),ncol(mer))]
					colnames(mer.1) <- paste0('Predicted_',round(confint.level_test*100,0),'CI_',c('Lower','Upper'),'_test_Pois')
				
					if(exists('O')) {
						O <- cbind(O,mer.1)
					} else {
						O <- mer.1
					}
				}
				
				# prediction dataset
				if(exists('pred.dataset')){
					# Probability
					temp.0 <- data.frame(Predicted_Pois=predict(res_Pois,newdata=pred.dataset,type='response'))
					temp <- data.frame(Predicted_Pois=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
					mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
					mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
					if(exists('O')) {
						O <- cbind(O,Predicted_pred_Pois=mer.1[,3])
					} else {
						O <- data.frame(Predicted_pred_Pois=mer.1[,3])
					}
					if(Predict_CI_pred){
						temp.0 <- as.data.frame(getIntervals(Pois_obj=res_Pois,newdat=pred.dataset,con.level=confint.level_pred))
						temp <- data.frame(temp=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
						mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
						mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),c((ncol(mer)-1),ncol(mer))]
						colnames(mer.1) <- paste0('Predicted_',round(confint.level_pred*100,0),'CI_',c('Lower','Upper'),'_pred_Pois')
					
						if(exists('O')) {
							O <- cbind(O,mer.1)
						} else {
							O <- mer.1
						}
					}
				}
				
				if(Part_index){
					temp <- rep('None',nrow(original.dataset))
					temp[rownames(original.dataset)%in%rownames(dataset)] <- 'Training'
					if(exists('test.dataset')) temp[rownames(original.dataset)%in%rownames(test.dataset)] <- 'Test'
					if(exists('pred.dataset')) temp[rownames(original.dataset)%in%rownames(pred.dataset)] <- 'Prediction'
					if(exists('O')) {
						O <- cbind(O,Partition_idx_Pois=temp,stringsAsFactors=F)
					} else {
						O <- data.frame(Partition_idx_Pois=temp,stringsAsFactors=F)
					}
				}

				if(exists('O'))	O[is.na(O)] <- ''

				### Default output
				# Data Structure
				R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file="./test.html")
				total.var <- ncol(raw.dat)
				used.var <- ifelse(is.null(vars),0,length(c(unique(unlist(strsplit(vars,":"))))))+length(c(Part.var,Pred.var))+length(offset)
				none.n <- nrow(original.dataset)-nrow(dataset)-ifelse(exists("test.dataset"),nrow(test.dataset),0)-ifelse(exists("pred.dataset"),nrow(pred.dataset),0)
				total.n <- paste0(nrow(original.dataset)," (Training: ",nrow(dataset),ifelse(exists("test.dataset"),paste0(", Test: ",nrow(test.dataset)),""),ifelse(exists("pred.dataset"),paste0(", Prediction: ",nrow(pred.dataset)),""),ifelse(none.n!=0,paste0(", None: ",none.n),""),")")
				
				DS <- matrix(c('Number of observations','Number of total variables','Number of used variables',total.n,total.var,used.var+1),ncol=2)
				R2HTML::HTML(DS,file="./test.html",align="left")
				if(exists('warn.DP4')) R2HTML::HTML(warn.DP4,file="./test.html")
				if(exists('warn.DP7')) R2HTML::HTML(warn.DP7,file="./test.html")
				if(exists('warn.DP8')) R2HTML::HTML(warn.DP8,file="./test.html")
				if(exists('warn.DP9')) R2HTML::HTML(warn.DP9,file="./test.html")
				if(exists('warn.DP2')) R2HTML::HTML(warn.DP2,file="./test.html")
				if(exists('warn.DP3')) R2HTML::HTML(warn.DP3,file="./test.html")

				## Varibale list
				R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file="./test.html")
				if(is.null(vars)) {
					Vars <- c()
				} else {
					Vars <-  unique(unlist(strsplit(vars,":")))
				}
				qual <- c(qual_var[qual_var%in%Vars],Part.var,Pred.var)
				quan <- c(quan_var[quan_var%in%Vars],res_var,offset)
				varlist <- matrix(c('Quantitative variable',paste(quan,collapse=', ')),ncol=2,byrow=T)
				if(length(qual)!=0) varlist <- rbind(varlist,c('Qualitative variable',paste(qual,collapse=', ')))
				R2HTML::HTML(varlist,file="./test.html",align="left")
				if(exists('warn.msg2')) R2HTML::HTML(warn.msg2, file="./test.html", row.names=FALSE)
				if(exists('warn.msg3')) R2HTML::HTML(warn.msg3, file="./test.html", row.names=FALSE)
				if(exists('warn.offset1')) R2HTML::HTML(warn.offset1, file="./test.html", row.names=FALSE)
				if(exists('warn.offset2')) R2HTML::HTML(warn.offset2, file="./test.html", row.names=FALSE)

				# Analysis Description
				R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file="./test.html")
				link_print	<- switch(link, log = "Log", identity = "Identity", sqrt = "Square root") ;
				AD <- matrix(c('Dependent variable',res_var),ncol=2,byrow=T)
				if(!is.null(vars)) AD <- rbind(AD,c('Explanatory variable',paste(vars,collapse=', ')))
				if(!is.null(offset)) AD <- rbind(AD,c('Offset variable',paste(offset,collapse=', ')))
				AD <- rbind(AD,matrix(c('Over-dispersion',over_disp,'Distribution',ifelse(over_disp,'Quasipoisson','Poisson'),'Intercept included',!noint,'Link function',link_print,'Variable selection',Select),ncol=2,byrow=T))
				if(Select) {
					direction <- switch(direct,forward="Forward selection",backward="Backward elimination",both="Stepwise regression")
					AD <- rbind(AD,matrix(c('Variable selection method',direction),ncol=2,byrow=T))
					if(!is.null(keep_var)) AD <- rbind(AD,c('Fixed variable for variable selection',paste(keep_var,collapse=', ')))
				}
				if(Valid.method=='Partition'){
					VM <- ifelse(!exists('test.dataset'),"Internal validation using training dataset","Internal validation using training/test dataset")
					if(exists('test.dataset')){
						if(Part.method=='percent'){
							PM <- paste0("Randomly split by percent (training: ",train.perc,"%, test: ",100-train.perc,"%)")
						} else if(Part.method=='variable'){
							PM <- paste0("Split by variable '",Part.var,"' (1: training, 2: test)")
						}
					}
				} else {
					VM <- ifelse(Cross.method=='LOOCV',"Leave-One-Out cross validation",paste0(k,"-fold cross validation"))
				}
				AD <- rbind(AD,matrix(c('Validation method',VM),ncol=2,byrow=T))
				if(exists('PM')) AD <- rbind(AD,matrix(c('Data splitting method for validataion',PM),ncol=2,byrow=T))
				if(exists('pred.dataset')){
					AD <- rbind(AD,matrix(c('Prediction using new dataset','TRUE','Data splitting method for prediction',paste0("Split by variable '",Pred.var,"' (1: prediction, 2: training/test)")),ncol=2,byrow=T))
				} else {
					AD <- rbind(AD,matrix(c('Prediction using new dataset','FALSE'),ncol=2,byrow=T))
				}
				R2HTML::HTML(AD,file="./test.html",align="left")
				if(exists('warn.VS1')) R2HTML::HTML(warn.VS1,file="./test.html")
 				if(exists('warn.VS2')) R2HTML::HTML(warn.VS2,file="./test.html")
				if(exists('warn.DP4')) R2HTML::HTML(warn.DP4,file="./test.html")
				if(exists('warn.DP7')) R2HTML::HTML(warn.DP7,file="./test.html")
				if(exists('warn.DP3')) R2HTML::HTML(warn.DP3,file="./test.html")
				


				## Analysis Results
				R2HTML::HTML(R2HTML::as.title("Results of Poisson Regression"),HR=2,file="./test.html")

				# Coefficient Estimate
				R2HTML::HTML(R2HTML::as.title("Coefficient Estimates"),HR=3,file="./test.html")
				VS1 <- capture.output(summary(fit_glm))
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

				CE <- as.data.frame(summary(res_Pois)$coef)
				colnames(CE) <- c('Estimate','SE','Z-value','P-value')
				if(CI){
					tmp <- merge(CE,confint.default(res_Pois, level=confint.level),by="row.names",all=TRUE,sort=F)
					rownames(tmp) <- tmp[,1]
					CE <- tmp[,-1]
					colnames(CE)[(ncol(CE)-1):ncol(CE)] <- paste0(c('Lower','Upper'),' bound of<br>',confint.level*100,'% CI for<br>Estimate')
				}

				if(exp_estim) {
					ORs <- exp(CE[,1,drop=F])
					colnames(ORs) <- 'exp(Estimate)'
					CE <- cbind(CE,ORs)
					CE <- CE[,c(1,ncol(CE),2:(ncol(CE)-1))]
					if(CI){
						tmp <- merge(CE,exp(confint.default(res_Pois, level=confint.level)),by="row.names",all=TRUE,sort=F)
						rownames(tmp) <- tmp[,1]
						CE <- tmp[,-1]
						colnames(CE)[(ncol(CE)-1):ncol(CE)] <- paste0(c('Lower','Upper'),' bound of<br>',confint.level*100,'% CI for<br>exp(Estimate)')
					}
				}

				if(VIF)	{
					if(is.null(vars)){
						warn.vif <- '<li> Warning : VIF is not supported for intercept-only model.'
					} else {
						VIF.1 <- rms::vif(res_Pois)
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
						vip <- varImp(res_Pois)
						if(!noint) {
						  vip <- rbind(NA,vip)
						  rownames(vip)[1] <- "(Intercept)"
						}

						tmp <- merge(CE,vip,by='row.names',all=TRUE)
						rownames(tmp) <- tmp[,1]
						colnames(tmp)[ncol(tmp)] <- "VIP"
						CE <- tmp[,-1]
					} else {
						warn.vip <- "<li> Warning : Variable importance table is not supported for intercept only model."
					}
				}
				R2HTML::HTML(Digits(CE),file="./test.html",align="left",digits=15)
				if(exists('warn.vif')) R2HTML::HTML(warn.vif, file="./test.html", row.names=FALSE)
				if(exists('warn.vip')) R2HTML::HTML(warn.vip, file="./test.html", row.names=FALSE)
				R2HTML::HTML(IF,file="./test.html", align="left",digits=15)
				
				#### Anova table, Goodness of fit, R-squared and classfication table
				# ANOVA=TRUE;GOF=TRUE;ng=10;R2=TRUE;classtab=TRUE;cutpoint=.5
				if(ANOVA){
					R2HTML::HTML(R2HTML::as.title("Analysis-of-Deviance Table"),HR=3,file="./test.html")
					if(length(vars)==0){
						warn.AT0 <- "<li> Warning : ANOVA table is not supported for intercept only model."
						R2HTML::HTML(warn.AT0,file="./test.html")
					} else if(noint){
						warn.AT1 <- "<li> Warning : ANOVA table is not supported for the model without intercept."
						R2HTML::HTML(warn.AT1,file="./test.html")
					} else {
						R2HTML::HTML(R2HTML::as.title("Model Effect (Goodness of Fit Test)"),HR=4,file="./test.html")
						form.null <- paste0(gsub('~.+','~ 1',form.0),ifelse(class(offset)=='character',paste0('+offset(',offset,')'),''))
						command_str<- paste0("null.fit<- try(glm(formula(",form.null,"), data=temp.dat, family=",ifelse(over_disp, 'quasipoisson', 'poisson'),"(",link,")))")
						eval(parse(text=command_str)) ;
						model.fit <- data.frame(anova(null.fit,res_Pois))
						model.fit <- model.fit[c(2,1),c(2,1,4,3)]
						model.fit <- cbind(model.fit,c(pchisq(model.fit[1,3],model.fit[1,4],lower.tail=F),NA))
						rownames(model.fit) <- c('Proposed model','Null model')
						colnames(model.fit) <- c('Deviance','DF(Deviacne)','Chisq','DF(Chisq)','P-value')
						R2HTML::HTML(Digits(model.fit),file="./test.html",align="left",digits=15)
						anova.msg <- "<li> Note : 'Null model' means the model including only intercept and 'Proposed model' means the model including all explanatory variables including interaction effect."
						R2HTML::HTML(anova.msg,file="./test.html")

						R2HTML::HTML(R2HTML::as.title(paste0("Variable Effect with Type ",ss_type," SS")),HR=4,file="./test.html")
						warn.desc <- ifelse(ss_type=='I',"<li> Note : In type I test, 'Null model' means the model including only intercept. Terms are added sequentially (first to last).",
								ifelse(ss_type=='II',"<li> Note : In type II test, 'Proposed model' means the model including all explanatory variables except interaction effect. The other rows are the testing results for each main effect after the other main effect.",
								"<li> Note : In type III test, 'Proposed model' means the model including all explanatory variables including interaction effect. The other rows are the testing results for each effect after the other effect."))
						
						if(ss_type=='I'){
							AT <- try(anova(res_Pois,test="Chisq"),s=T)
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

						if(ss_type %in% c('II','III')){
							# II type : ignore interaction term
							# III type : calculate SS including interaction term
							RN <- vars
							if(ss_type=='II' & length(grep(':',RN))>0) {
								RN <- RN[-grep(':',RN)]
								warn.AT2 <- '<li> Warning : Test for interaction effects is not provided in type II test. Use type III test for interaction effects.'
							}
							# warn.AT3 <- '<li> Note : If there is indeed no interaction, then type II is statistically more powerful than type III.'

							AT.form.full <- paste(res_var,'~',paste(RN,collapse=' + '),ifelse(class(offset)=='character',paste0('+ offset(',offset,')'),''))
							command_str<- paste0("full.fit<- try(glm(formula(AT.form.full), data=temp.dat, family=",ifelse(over_disp, 'quasipoisson', 'poisson'),"(",link,")))")
							eval(parse(text=command_str))
							options(contrasts=c("contr.sum", "contr.poly"))
							res <- try(car::Anova(full.fit,type='III',test="LR"))

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
								warn.AT4 <- "<li> Error : Fail to fit the Model."
								R2HTML::HTML(warn.AT4,file="./test.html")
							}
						}
					}
				}

				## Goodness of fit
				R2HTML::HTML(R2HTML::as.title("Model Fitness Measurements"), HR=3, file="./test.html") ;
				fitness_print	<- data.frame(numeric(0), numeric(0)) ;
				fitness_print	<- rbind(fitness_print, c(sum(residuals(fit_glm, type="deviance")^2), fit_glm$df.residual)) ;
				fitness_print	<- rbind(fitness_print, c(sum(residuals(fit_glm, type="pearson")^2), fit_glm$df.residual)) ;
				if(!over_disp){
					LL <- logLik(fit_glm)
					fitness_print	<- rbind(fitness_print, c(-2*LL[1], attr(LL,'df'))) ;
					fitness_print	<- rbind(fitness_print, c(AIC(fit_glm), NA)) ;
					fitness_print	<- rbind(fitness_print, c(BIC(fit_glm), NA)) ;
					row.names(fitness_print)	<- c("Deviance", "Pearson's chi-square", "-2*log-likelihood", "AIC", "BIC") ;
				} else {
					row.names(fitness_print)	<- c("Deviance", "Pearson's chi-square") ;
				}
				names(fitness_print)	<- c("Value", "DF") ;
				R2HTML::HTML(Digits(fitness_print), file="./test.html", digits=15, align="left",caption="<div style='text-align:left'> <li> A model with a smaller value is better.") ;

				if(GOF){
					R2HTML::HTML(R2HTML::as.title("Goodness of Fit Test"),HR=3,file="./test.html")
					if(length(vars)!=0){
						## log-liklehood ratio
						R2HTML::HTML(R2HTML::as.title("Analysis-of-Deviance Table (Likelihood Ratio Test)"),HR=4,file="./test.html")
						if(!noint){
							## log-liklehood ratio test
							form.null <- paste0(gsub('~.+','~ 1',form.0),ifelse(class(offset)=='character',paste0('+offset(',offset,')'),''))
							command_str<- paste0("null.fit<- try(glm(formula(",form.null,"), data=temp.dat, family=",ifelse(over_disp, 'quasipoisson', 'poisson'),"(",link,")))")
							eval(parse(text=command_str)) ;
							model.fit <- data.frame(anova(null.fit,res_Pois))
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
					} else {
						R2HTML::HTML("<li> Warning : Goodness of fit test is not supported for intercept only model.",file="./test.html")
					}
				}

				## Overdispersion Test (ref: http://127.0.0.1:15678/library/AER/html/dispersiontest.html)
				if(ODtest) {
					R2HTML::HTML(R2HTML::as.title("Test for Overdispersion"),HR=3,file="./test.html",append=TRUE)
					load.pkg('AER')
					ODres <- AER::dispersiontest(glm(form.1, data=temp.dat, family=poisson(link=link_actual)))
					ODres1 <- data.frame(ODres$estimate, ODres$statistic, ODres$p.value)
					colnames(ODres1) <- c('Dispersion','Z-value','P-value')
					R2HTML::HTML(Digits(ODres1),file="./test.html",align="left",digits=15,row.names=F)
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
					# Function that returns Root Mean Squared Error
					rmse <- function(error){
						sqrt(mean(error^2,na.rm=T))
					}
					# Function that returns MAE
					mae <- function(error){
						mean(abs(error),na.rm=T)
					}
					# Function that returns R-squared (correlation)
					rsq <- function(pred,obs){
						cor(pred,obs,use='complete.obs')
					}
					
					if(Valid.method=='Partition'){
						if(exists('test.dataset')){
							test.Predicted <- predict(res_Pois,newdata=test.dataset,type='response')
							test.Observed <- test.dataset[,res_var]
							test.Resid <- test.Predicted-test.Observed
							
							test.RMSE <- rmse(test.Resid)
							test.MAE <- mae(test.Resid)
							test.Rsq <- rsq(test.Predicted,test.Observed)
							test.n <- sum(!is.na(test.Resid))
						} else {
							train.Perc <- 100
						}
						train.Observed <- temp.dat[,res_var]
						train.Predicted <- predict(res_Pois,type='response')
						train.Resid <- train.Predicted-train.Observed
						
						train.RMSE <- rmse(train.Resid)
						train.MAE <- mae(train.Resid)
						train.Rsq <- rsq(train.Predicted,train.Observed)
						train.n <- sum(!is.na(train.Resid))
						
						train.Perc <- ifelse(exists('test.dataset'),train.n/(train.n+test.n)*100,100)

						vali <- data.frame(N.observed=train.n,Percent=train.Perc,RMSE=train.RMSE,Rsquared=train.Rsq)
						rownames(vali) <- 'Training'
						if(exists('test.dataset')){
							vali <- rbind(vali,data.frame(N.observed=test.n,Percent=(100-train.Perc),RMSE=test.RMSE,Rsquared=test.Rsq))
							rownames(vali)[2] <- 'Test'
						}
						colnames(vali)[1] <- 'N.non-missing<br>observations'
						R2HTML::HTML(Digits(vali),file="./test.html",align="left",digits=15)
						if(exists('warn.DP4')) R2HTML::HTML(warn.DP4,file="./test.html")
						if(exists('warn.DP7')) R2HTML::HTML(warn.DP7,file="./test.html")
						if(exists('warn.DP8')) R2HTML::HTML(warn.DP8,file="./test.html")
					} else {
						ctrl <- caret::trainControl(method = ifelse(k==nrow(temp.dat)|Cross.method=='LOOCV',"LOOCV","cv"), number = k)
						newdat <- temp.dat[,c(res_var,quan_var,qual_var,offset),drop=F]
						command_str	<- paste0("mod_fit <- train(res_Pois$formula,  data=newdat, method='glm',family=",ifelse(over_disp, 'quasipoisson', 'poisson'),"(",link,"),trControl = ctrl, tuneLength = 5)")
						eval(parse(text=command_str)) ;
						A8 <- mod_fit$results
						if(k==nrow(temp.dat)|Cross.method=='LOOCV'){
							res.acc.II <- A8[,c('RMSE','Rsquared')]
						} else {
							res.acc.II <- A8[,c('RMSE','RMSESD','Rsquared','RsquaredSD')]
							colnames(res.acc.II)[c(2,4)] <- c('SD(RMSE)','SD(Rsquared)')
						}
						R2HTML::HTML(Digits(res.acc.II),file="./test.html",align="left",digits=15,row.names=F)
						if(exists('warn.valid')) R2HTML::HTML(warn.valid,file="./test.html")
					} 
				} else {
					warn.msg10 <- "<li> Warning : Validation is not supported for intercept only model."
					R2HTML::HTML(warn.msg10,file="./test.html")
				}
				
				## Variable selection
				if(Select) stepAIC.wj(object=fit_glm,dataset=temp.dat,dep_var=res_var,type='poisson',noint=noint,direct=direct,keep_var=keep_var,hr=3,vars=vars,link=link_actual,offset=offset,dist=dist,CI=CI,confint.level=CI_level,VIF=VIF,exp_estim=exp_estim)
			}
		}

		# Used R packages
		R2HTML::HTML(R2HTML::as.title("Used R Packages"),HR=2,file="./test.html")
		pkg.list <- list(list("Poisson regression","glm","stats"),
				 list("Confidence interval for regression coefficients","confint","stats"),
				 list("Variance inflation factor","vif","rms"),
				 list("ANOVA table",c("anova","Anova"),c("stats","car")),
				 list("Model fitness measurements","residuals, logLik, AIC, BIC","stats"),
				 list("Goodness of fit test",'anova','stats'),
				 list("Over-dispersion test","dispersiontest","AER"),
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
		
		R2HTML::HTMLhr(file="./test.html") ;
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". Rex : Poisson Regression",sep=""),file="./test.html") ;
		R2HTML::HTMLhr(file="./test.html") ;
	})

	if(exists('O')) {
		return(list(html=html.output,Output=O))
	} else {
		return(html.output)
	}
}
