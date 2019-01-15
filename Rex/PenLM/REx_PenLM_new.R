# penalized linear regression

REx_PenLM <- function(dataset,dep_var,indep_cat_var=NULL,indep_numeric_var=NULL,
  vars=NULL,noint=FALSE,standardize=TRUE,
  Penalty=c('Ridge','Lasso','EN'),alpha=0.5,
  TuningPara=c('Grid','Custom'),Grid.Num=NULL,Custom.List=NULL,
  Cross.method=c('KFOLD','LOOCV'),k=10,
  AccMS=c('MSE','MAE'),
  Part.method=c('all','percent','variable'),
  train.perc=70,Part.var=NULL,
  Profile=FALSE,Best_model_print=TRUE,CI=TRUE,confint.level=0.95,
  VIF=FALSE,ANOVA=TRUE,ss=c('I','II','III'),GOF=FALSE,Plot=FALSE,
  Best_model_save=FALSE,Predict_train=FALSE,Predict_CI_train=FALSE,
  Predict_PI_train=FALSE,confint.level_train=0.95,Resid=FALSE,
  stdResid=FALSE,studResid=FALSE,cook_distance=FALSE,hat_value=FALSE,
  Predict_test=FALSE,Predict_CI_test=FALSE,Predict_PI_test=FALSE,confint.level_test=0.95,Part_index=FALSE){

  load.pkg(c("R2HTML", "rms", "MASS","glmnet","plotmo"))
  
  #
  # "./test.html" => stdout() 
  #
  
  html.output <- capture.output({
    # Title
    R2HTML::HTML(R2HTML::as.title("Penalized Linear Regression"),
                 HR=1,file="./test.html",append=FALSE)

    ## Warnings
    # Response variable type
    if(!is.numeric(dataset[,dep_var])) warn.msg1 <- '<li> Error : Dependent variable should be numeric. Analysis has been stopped.'

    # explanatory variable type
    if(is.null(vars)) {
      Vars <- c()
    } else {
      Vars <-  vars
    }
    indep_cat_var <- indep_cat_var[indep_cat_var%in%Vars]
    if(length(indep_cat_var)==0) indep_cat_var <- NULL
    indep_numeric_var <- indep_numeric_var[indep_numeric_var%in%Vars]
    if(length(indep_numeric_var)==0) indep_numeric_var <- NULL

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
      Vars <-  vars
    }
    indep_cat_var <- indep_cat_var[indep_cat_var%in%Vars]
    if(length(indep_cat_var)==0) indep_cat_var <- NULL
    indep_numeric_var <- indep_numeric_var[indep_numeric_var%in%Vars]
    if(length(indep_numeric_var)==0) indep_numeric_var <- NULL
    
    # no intercept & no explanatory variable : error
    if(noint & is.null(vars)) 
      warn.msg4 <- '<li> Error : With no intercept, at least 1 independent variable should be selected. Analysis has been stopped.'

    ## Data processing
    original.dataset <- dataset
    
    # training & test
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
          
        if(nrow(original.dataset)!=(nrow(dataset)+nrow(test.dataset)))  warn.DP9 <- "<li> Warning : Observations with missing dependent variable were not assigned to training dataset or test dataset."
          
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
        test.dataset <- dataset[dataset[,Part.var]==2,]  # test dataset
        dataset <- dataset[dataset[,Part.var]==1,]  # training dataset
      }
    }
    #### Then, original.dataset : original dataset, dataset : training dataset, test.dataset : test dataset.

    # independent variables should be a matrix with 2 or more columns 
    temp.form <- paste0(dep_var,' ~ ', paste(vars, collapse=' + '))
    if(noint) temp.form <- paste0(temp.form,-1)
    mm <- model.matrix(as.formula(temp.form),data=dataset)
    if(ncol(mm)<2) warn.msg5 <- '<li> Error : Independent variables should be 2 or more.Analysis has been stopped.'

    # warnings
    if(exists('warn.msg1')|exists('warn.msg4')|exists('warn.msg5')|exists('warn.DP1')|exists('warn.DP6')){
      R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file="./test.html")
      if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file="./test.html")
      if(exists('warn.msg4')){
        if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file="./test.html")
        R2HTML::HTML(warn.msg4,file="./test.html")
      }
      if(exists('warn.msg5')) R2HTML::HTML(warn.msg5,file="./test.html")
      if(exists('warn.DP1')) R2HTML::HTML(warn.DP1,file="./test.html")
      if(exists('warn.DP6')) R2HTML::HTML(warn.DP6,file="./test.html")
    } else {
      var_info <- c(dep_var,indep_numeric_var,indep_cat_var)
      raw.dat <- dataset
      temp.dat <- dataset[complete.cases(dataset[,var_info,drop=F]),,drop=F]

      ## Model formula
      form.0 <- paste0(dep_var,' ~ ', paste(vars, collapse=' + '))
      if(noint) form.0 <- paste(form.0, -1)
      form.1 <- as.formula(form.0)
      
      ## Penalty & alpha
      if(Penalty=='Ridge') alpha=0
      if(Penalty=='Lasso') alpha=1

      ###
      ## first Model fitting : res_PLM_0
      ## inference => res_PLM
      #res_PLM <- suppressWarnings(try(cv_glmnet_wrapper(form.1, data=dataset, 
      #  family="gaussian", alpha=alpha, nlambda=Grid.Num, lambda=Custom.List, 
      #  standardize=standardize, intercept=intercept)))
      if(TuningPara=='Custom'){
        Grid.Num <- NULL
      } else {
        Custom.List <- NULL
      }
      res_PLM_0 <- cv_glmnet_wrapper(form.1, data=temp.dat, lambda=Custom.List, 
        type.measure=ifelse(AccMS=='MSE',"mse",'mae'), nfolds=k, 
        family="gaussian", alpha=alpha, nlambda=Grid.Num, 
        standardize=standardize, intercept=!noint, opt.lambda="lambda.min")
      res_PLM <- res_PLM_0  

      ## p>n case
      if(length(vars)>nrow(temp.dat) & Best_model_print & Best_model_save){
        warn.largerp <- '<li> Warning : Since the number of parameters of the final model is larger than the number of observations, results of linear regression are not provided.'
		    Best_model_print <- Best_model_save <- FALSE
      }

      if(class(res_PLM_0)=='try-error'){
        R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file="./test.html")
        warn.msg.conv <- "<li> Warning : Penalized Linear regression was not fitted. Analysis has been stopped."
        R2HTML::HTML(warn.msg.conv,file="./test.html")
      } else {
        # Data Structure
        R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file="./test.html")
        total.var <- ncol(original.dataset)
        used.var <- ifelse(is.null(vars),0,length(c(unique(unlist(strsplit(vars,":"))),Part.var)))
        none.n <- nrow(original.dataset)-nrow(dataset)-ifelse(exists("test.dataset"),nrow(test.dataset),0)
        total.n <- paste0(nrow(original.dataset)," (Number of missings: ",sum(!complete.cases(original.dataset[,c(var_info,Part.var),drop=F])),")")
        DS <- matrix(c('Number of observations','Number of total variables','Number of used variables',total.n,total.var,used.var+1),ncol=2)
        R2HTML::HTML(DS,file="./test.html",align="left")
        if(exists('warn.DP7')) R2HTML::HTML(warn.DP7,file="./test.html")
        if(exists('warn.DP8')) R2HTML::HTML(warn.DP8,file="./test.html")
        if(exists('warn.DP9')) R2HTML::HTML(warn.DP9,file="./test.html")

        # Varibale list
        R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file="./test.html")
        if(is.null(vars)) {
          Vars <- c()
        } else {
          Vars <-  vars
        }
        qual <- c(indep_cat_var[indep_cat_var%in%Vars],Part.var)
        quan <- c(indep_numeric_var[indep_numeric_var%in%Vars],dep_var)
        varlist <- matrix(c('Quantitative variable',paste(quan,collapse=', ')),ncol=2,byrow=T)
        if(length(qual)!=0) varlist <- rbind(varlist,c('Qualitative variable',paste(qual,collapse=', ')))
        R2HTML::HTML(varlist,file="./test.html",align="left")
        if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file="./test.html")
        if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file="./test.html")

        ## Analysis Description
        R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file="./test.html")
        AD <- matrix(c('Dependent variable',dep_var),ncol=2,byrow=T)
        if(!is.null(vars)) AD <- rbind(AD,c('Explanatory variable',paste(vars,collapse=', ')))
        AD <- rbind(AD,matrix(c('Intercept included',!noint,'Standardization for explanatory variables',standardize,'Penalty method',switch(Penalty,Lasso='Lasso',Ridge='Ridge',EN='Elastic net')),ncol=2,byrow=T))
        if(Penalty=='EN') AD <- rbind(AD,c('- Alpha for Elastic net',alpha))
        AD <- rbind(AD,matrix(c('Candidate tuning parameter',switch(TuningPara,Custom='User-defined list',Grid=paste0('Grid search (Number of grids : ',Grid.Num,')')),
                                '- Cross validation method',switch(Cross.method,KFOLD=paste0(k,'-fold cross validation'),LOOCV='Leave-One-Out cross validation'),
                                '- Accuracy measure',switch(AccMS,MSE='Mean squared error',MAE='Mean absolute error')),ncol=2,byrow=T))
        # validation method
        VM <- ifelse(!exists('test.dataset'),"Validation using training dataset","Validation using training/test dataset")
        if(exists('test.dataset')){
          if(Part.method=='percent'){
            PM <- paste0("Randomly split by percent (training: ",train.perc,"%, test: ",100-train.perc,"%)")
          } else if(Part.method=='variable'){
            PM <- paste0("Split by variable '",Part.var,"' (1: training, 2: test)")
          }
        }
        AD <- rbind(AD,matrix(c('Validation method',VM),ncol=2,byrow=T))
        if(exists('PM')) AD <- rbind(AD,matrix(c('- Data splitting method for validataion',PM),ncol=2,byrow=T))
        R2HTML::HTML(AD,file="./test.html",align="left")
        if(exists('warn.DP7')) R2HTML::HTML(warn.DP7,file="./test.html")
        if(exists('warn.DP3')) R2HTML::HTML(warn.DP3,file="./test.html")

        ##### Main Results
        R2HTML::HTML(R2HTML::as.title("Results of Penalized Linear Regression"),HR=2,file="./test.html")
        
        if(Profile){
          # solution path according to tuning parameter lambda
          R2HTML::HTML(R2HTML::as.title("Solution Path for the Best model"),HR=3,file="./test.html")
          REx_ANA_PLOT()
          plotmo::plot_glmnet(res_PLM_0$glmnet, xvar="norm", col=rainbow(length(vars),v=0.85),label=TRUE)
          REx_ANA_PLOT_OFF('')

          REx_ANA_PLOT()
          plotmo::plot_glmnet(res_PLM_0$glmnet, xvar="lambda", col=rainbow(length(vars),v=0.85), label=TRUE)
          REx_ANA_PLOT_OFF('')
          
          REx_ANA_PLOT()
          plotmo::plot_glmnet(res_PLM_0$glmnet, xvar="dev", col=rainbow(length(vars),v=0.85), label=TRUE)
          REx_ANA_PLOT_OFF('')
          
          R2HTML::HTML(R2HTML::as.title("Cross-Validation Plot"),HR=3,file="./test.html")
          REx_ANA_PLOT()
          glmnet::plot.cv.glmnet(res_PLM_0)
          REx_ANA_PLOT_OFF('')
        } 
        
        #### validataion ####
        R2HTML::HTML(R2HTML::as.title(VM),HR=3,file="./test.html")
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
          
          if(exists('test.dataset')){
            #
            # should check the standarization !!!
            #
            test.newx <- model.matrix(form.1,data=test.dataset)
            test.Predicted <- predict_glmnet_wrapper(res_PLM_0$glmnet, 
                                                     newx=test.newx, s=res_PLM_0$lambda.min, type='response')
            #R2HTML::HTML(R2HTML::as.title(test.Predicted),HR=3,file="./test.html")
            test.Observed <- test.dataset[,dep_var]
            test.Resid <- test.Predicted-test.Observed
            
            test.RMSE <- rmse(test.Resid)
            test.MAE <- mae(test.Resid)
            test.Rsq <- rsq(test.Predicted,test.Observed)
            test.n <- sum(!is.na(test.Resid))
          } else {
            train.Perc <- 100
          }
          train.Observed <- temp.dat[,dep_var]
          
          train.newx <- model.matrix(form.1,data=temp.dat)
          train.Predicted <- predict_glmnet_wrapper(res_PLM_0$glmnet, 
                                                    newx=train.newx, s=res_PLM_0$lambda.min, type='response')
          # train.Predicted : n by nlambda
          train.Resid <- train.Predicted-train.Observed
          train.RMSE <- rmse(train.Resid)
          train.MAE <- mae(train.Resid)
          train.Rsq <- rsq(train.Predicted,train.Observed)
          train.n <- sum(!is.na(train.Resid))
          
          train.Perc <- ifelse(exists('test.dataset'),train.n/(train.n+test.n)*100,100)
          
          vali <- data.frame(N.observed=train.n,Percent=train.Perc,RMSE=train.RMSE,
                             MAE=train.MAE,Rsquared=train.Rsq)
          rownames(vali) <- 'Training'
          if(exists('test.dataset')){
            vali <- rbind(vali,data.frame(N.observed=test.n,Percent=(100-train.Perc),
                                          RMSE=test.RMSE,MAE=test.MAE,Rsquared=test.Rsq))
            rownames(vali)[2] <- 'Test'
          }
          colnames(vali)[1] <- 'N.non-missing<br>observations'
          R2HTML::HTML(Digits(vali),file="./test.html",align="left",digits=15)
          if(exists('warn.DP7')) R2HTML::HTML(warn.DP7,file="./test.html")
          if(exists('warn.DP8')) R2HTML::HTML(warn.DP8,file="./test.html")
        } else {
          warn.msg10 <- "<li> Warning : Validation is not supported for intercept only model."
          R2HTML::HTML(warn.msg10,file="./test.html")
        }

        #### Regression Results
        ## res_PLM: final model
        #!!! R2HTML::HTML(R2HTML::as.title(Best_model_print),HR=3,file="./test.html")
        
        if(exists('warn.largerp')){
          R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file="./test.html")
          R2HTML::HTML(warn.largerp,file="./test.html")
        }
        
        if(Best_model_print){
          # Coefficient Estimate
          R2HTML::HTML(R2HTML::as.title("Results of Linear Regression for the Best Model"),HR=2,file="./test.html")

          ###
          ypos <- match(dep_var, colnames(temp.dat))
          coef <- predict_glmnet_wrapper(res_PLM_0$glmnet.fit, newx=temp.dat[, -ypos], 
                    s=res_PLM_0$lambda.min, type="coefficients")
          
          vars_ols <- coef@Dimnames[[1]][which(coef!=0)]
          vars_ols <- vars_ols[!vars_ols%in%'(Intercept)']
          if(length(vars_ols)==0) vars_ols <- NULL

          ## Model formula
          form.0.ols <- paste0(dep_var,' ~ ', ifelse(is.null(vars_ols),1,paste(vars_ols, collapse=' + ')))
          
          if(noint) form.0.ols <- paste(form.0.ols, -1)
          form.1.ols <- as.formula(form.0.ols)
          res_PLM <- lm(form.1.ols, data=temp.dat)

          R2HTML::HTML(R2HTML::as.title("Coefficient Estimates"),HR=3,file=stdout())
          CE <- as.data.frame(summary(res_PLM)$coef)
          colnames(CE) <- c('Estimate','SE','T-value','P-value')
          if(CI){
            tmp <- merge(CE,confint(res_PLM, level=confint.level),by="row.names",all=TRUE,sort=F)
            rownames(tmp) <- tmp[,1]
            CE <- tmp[,-1]
            colnames(CE)[(ncol(CE)-1):ncol(CE)] <- paste0(c('Lower','Upper'),' bound of<br>',confint.level*100,'% CI')
          }
          if(VIF)	{
            if(is.null(vars_ols)){
              warn.VIF <- '<li> Warning : VIF is not supported for intercept-only model.'
            } else {
              VIF.1 <- rms::vif(res_PLM)
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

          R2HTML::HTML(Digits(CE),file="./test.html",align="left",digits=15)
          if(exists('warn.VIF')) R2HTML::HTML(warn.VIF,file="./test.html")
          
          # Anova table and R-squared
          if(ANOVA){
            R2HTML::HTML(R2HTML::as.title("Analysis-of-Variance Table"),HR=3,file="./test.html")
            if(length(vars_ols)==0){
              warn.AT1 <- "<li> Warning : ANOVA table is not supported for intercept only model."
              R2HTML::HTML(warn.AT1,file="./test.html")
            } else if(noint){
              warn.AT2 <- "<li> Warning : ANOVA table is not supported for the model without intercept."
              R2HTML::HTML(warn.AT2,file="./test.html")
            } else {
              R2HTML::HTML(R2HTML::as.title("Model Effect (Goodness of Fit Test)"),HR=4,file="./test.html")
              Anova <- as.matrix(anova(lm(formula(paste(dep_var,'~1',sep='')),x=TRUE,data=temp.dat),res_PLM))
              A1 <- rbind(Anova[2,4:3], Anova[2:1,c(2,1)])
              MS <- c(A1[1:2,1]/A1[1:2,2],NA)
              A2 <- data.frame(A1,MS=MS,Fvalue=c(Anova[2,5],NA,NA),Pvalue=c(Anova[2,6],NA,NA),R2=c(summary(res_PLM)$r.squared,NA,NA),adj.R2=c(summary(res_PLM)$adj.r.squared,NA,NA))
              colnames(A2) <- c('SS','DF','MS','F-value','P-value','R2','adj.R2')
              rownames(A2) <- c('Regression','Residual','Total')
              R2HTML::HTML(Digits(A2),file="./test.html",align="left",digits=15)
              
              R2HTML::HTML(R2HTML::as.title(paste0("Variable Effect with Type ",ss," SS")),HR=4,file="./test.html")
              warn.desc <- ifelse(ss=='I',"<li> Note : In type I test, terms are added sequentially (first to last).",
                                  ifelse(ss=='II',"<li> Note : In type II test, each row is the testing result for each main effect after the other main effect.",
                                         "<li> Note : In type III test, each row is the testing result for each effect after the other effect."))
              
              if(ss=='I'){
                AT <- try(anova(res_PLM,test="Chisq"),s=T)
                rowNames <- rownames(AT)
                if(class(AT)[1]!='try-error'){
                  AT <- data.frame(AT)
                  AT <- AT[,c(2,1,3,4,5)]
                  rownames(AT) <- rowNames
                  colnames(AT) <- c('SS','DF','MS','F-value','P-value')
                  R2HTML::HTML(Digits(AT),file="./test.html",align="left",digits=15)
                  R2HTML::HTML(warn.desc,file="./test.html")
                } else {
                  warn.AT3 <- "<li> Error : Fail to fit the Model."
                  R2HTML::HTML(warn.AT3,file="./test.html")
                }
              }
              #R2HTML::HTML(var_info_ols,file="./test.html",align="left",digits=15)
              if(ss %in% c('II','III')){
                # II type : ignore interaction term
                # III type : calculate SS including interaction term
                RN <- vars_ols
                if(ss=='II' & length(grep(':',RN))>0) {
                  RN <- RN[-grep(':',RN)]
                  warn.AT4 <- '<li> Warning : Test for interaction effects is not provided in type II test. Use type III test for interaction effects.'
                }
                # warn.AT5 <- '<li> Note : If there is indeed no interaction, then type II is statistically more powerful than type III.'
                
                AT.form.full <- paste(dep_var,'~',paste(RN,collapse=' + '))
                full.fit <- lm(formula(AT.form.full),data=temp.dat)
                
                options(contrasts=c("contr.sum", "contr.poly"))
                res <- try(car::Anova(full.fit,type='III'))
                
                if(class(res)[1]!='try-error'){
                  res <- data.frame(res)
                  MS <- res[,1]/res[,2]
                  AT <- data.frame(res[,c(1:2)],MS,res[,c(3,4)])
                  colnames(AT) <- c('SS','DF','MS','F-value','P-value')
                  
                  R2HTML::HTML(Digits(AT[-nrow(AT),]),file="./test.html",align="left",digits=15)
                  R2HTML::HTML(warn.desc,file="./test.html")
                  if(exists('warn.AT4')) R2HTML::HTML(warn.AT4,file="./test.html")
                } else {
                  warn.AT6 <- "<li> Error : Fail to fit the Model."
                  R2HTML::HTML(warn.AT6,file="./test.html")
                }
              }
            }
          }
          
          R2HTML::HTML(R2HTML::as.title("Model Fitness Measurements"),HR=3,file="./test.html")
          fitness_print	<- data.frame(numeric(0), numeric(0)) ;
          fitness_print	<- rbind(fitness_print, c(sum(residuals(res_PLM, type="deviance")^2), res_PLM$df.residual)) ;
          fitness_print	<- rbind(fitness_print, c(sum(residuals(res_PLM, type="pearson")^2), res_PLM$df.residual)) ;
          LL <- logLik(res_PLM)
          fitness_print	<- rbind(fitness_print, c(-2*LL[1], attr(LL,'df'))) ;
          fitness_print	<- rbind(fitness_print, c(AIC(res_PLM), NA)) ;
          fitness_print	<- rbind(fitness_print, c(BIC(res_PLM), NA)) ;
          names(fitness_print)	<- c("Value", "DF") ;
          row.names(fitness_print)	<- c("Deviance", "Pearson's chi-square", "-2*log-likelihood", "AIC", "BIC") ;
          R2HTML::HTML(Digits(fitness_print), file="./test.html", 
                       align="left",digits=15,caption="<div style='text-align:left'> <li> A model with a smaller value is better.") ;
          
          if(GOF){
            R2HTML::HTML(R2HTML::as.title("Goodness of Fit Test (Likelihood Ratio Test)"),HR=3,file="./test.html")
            if(length(vars_ols)!=0){
              ## log-liklehood ratio
              if(!noint){
                null.model <- lm(formula(paste(dep_var,"~1")),data=temp.dat)
                GOF1 <- anova(null.model, res_PLM, test ="Chisq")
                A3 <- as.data.frame(GOF1)[,c(2,1,4,3,5)]
                rownames(A3) <- c("Null Model","Proposed Model")
                colnames(A3) <- c('RSS','DF(RSS)','Chisq','DF(Chisq)','P-value')
                R2HTML::HTML(Digits(A3),file="./test.html",align="left",digits=15)
              } else {
                R2HTML::HTML('<li> Warning : Likelihood Ratio Test is not supported for the model without intercept.',file="./test.html")
              }
            } else {
              warn.GOF <- "<li> Warning : Goodness of fit test is not supported for intercept only model."
              R2HTML::HTML(warn.GOF,file="./test.html")
            }
          }
          
          # plot
          if(Plot) {
            R2HTML::HTML(R2HTML::as.title("Graphs for Regression Diagnostics"),
                         HR=3,file="./test.html",append=TRUE)
            load.pkg("ggfortify")
            REx_ANA_PLOT(800,1200)
            print(autoplot(res_PLM, which=1:6) + theme_bw() + 
                    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), plot.title = element_text(hjust = 0.5)))
            REx_ANA_PLOT_OFF("")
          }
        } ###Best_model_end

        ### Save results in the Excel spread sheet 
        if(Best_model_save){
          if(Predict_train) {
            temp.0 <- data.frame(Fitted_PLM=predict(res_PLM, newdata=temp.dat))
            temp <- data.frame(Fitted_PLM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
            mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
            mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
            O <- data.frame(Fitted_train_PLM=mer.1[,3])
          }
          if(Predict_CI_train){
            temp.0 <- data.frame(predict(res_PLM,newdata=temp.dat,interval='confidence',level=confint.level_train))
            temp <- data.frame(temp=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
            mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
            mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),c((ncol(mer)-1),ncol(mer))]
            colnames(mer.1) <- paste0('Fitted_',round(confint.level_train*100,0),'CI_',c('Lower','Upper'),'_train_PLM')
            if(exists('O')) {
              O <- cbind(O,mer.1)
            } else {
              O <- mer.1
            }
          }
          if(Predict_PI_train){
            temp.0 <- data.frame(predict(res_PLM,newdata=temp.dat,interval='prediction',level=confint.level_train))
            temp <- data.frame(temp=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
            mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
            mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),c((ncol(mer)-1),ncol(mer))]
            colnames(mer.1) <- paste0('Fitted_',round(confint.level_train*100,0),'PI_',c('Lower','Upper'),'_train_PLM')
            if(exists('O')) {
              O <- cbind(O,mer.1)
            } else {
              O <- mer.1
            }
          }
          if(Resid) {
            temp.0 <- data.frame(Resid_PLM=resid(res_PLM))
            temp <- data.frame(Resid_PLM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
            mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
            mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
            if(exists('O')) {
              O <- cbind(O,unstdResid_train_PLM=mer.1[,3])
            } else {
              O <- data.frame(unstdResid_train_PLM=mer.1[,3])
            }
          }
          if(stdResid) {
            temp.0 <- data.frame(stdResid_PLM=MASS::stdres(res_PLM))
            is.singular <- "try-error"%in%class(try(temp.0[,1],silent=T))
            if(is.singular){
              warn.res1 <- "<li> Error : Cannnot calculate standardized residuals because of the (near) linear dependences in explatory variables. Please check multicollinearity problem."
            } else {
              temp <- data.frame(stdResid_PLM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
              mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
              mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
              if(exists('O')) {
                O <- cbind(O,stdResid_train_PLM=mer.1[,3])
              } else {
                O <- data.frame(stdResid_train_PLM=mer.1[,3])
              }
            }
          }
          if(studResid) {
            temp.0 <- data.frame(studResid_PLM=MASS::studres(res_PLM))
            is.singular <- "try-error"%in%class(try(temp.0[,1],silent=T))
            if(is.singular){
              warn.res2 <- "<li> Error : Cannnot calculate standardized residuals because of the (near) linear dependences in explatory variables. Please check multicollinearity problem."
            } else {
              temp <- data.frame(studResid_PLM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
              mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
              mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
              if(exists('O')) {
                O <- cbind(O,studResid_train_PLM=mer.1[,3])
              } else {
                O <- data.frame(studResid_train_PLM=mer.1[,3])
              }
            }
          }
          if(cook_distance){
            temp.0 <- data.frame(CookDist_PLM=cooks.distance(res_PLM))
            temp <- data.frame(CookDist_PLM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
            mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
            mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
            if(exists('O')) {
              O <- cbind(O,CookDist_train_PLM=mer.1[,3])
            } else {
              O <- data.frame(CookDist_train_PLM=mer.1[,3])
            }
          }
          if(hat_value){
            temp.0 <- data.frame(HatValue_PLM=hatvalues(res_PLM))
            temp <- data.frame(HatValue_PLM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
            mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
            mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
            if(exists('O')) {
              O <- cbind(O,HatValue_train_PLM=mer.1[,3])
            } else {
              O <- data.frame(HatValue_train_PLM=mer.1[,3])
            }
          }
        
          # test dataset
          if(exists('test.dataset') & Predict_test) {
            temp.0 <- data.frame(Fitted_PLM=predict(res_PLM,newdata=test.dataset))
            temp <- data.frame(Fitted_PLM=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
            mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
            mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),]
            if(exists('O')) {
              O <- cbind(O,Predicted_test_PLM=mer.1[,3])
            } else {
              O <- data.frame(Predicted_test_PLM=mer.1[,3])
            }
          }
          if(exists('test.dataset')&Predict_CI_test){
            temp.0 <- data.frame(predict(res_PLM,newdata=test.dataset,interval='confidence',level=confint.level_test))
            temp <- data.frame(temp=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
            mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
            mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),c((ncol(mer)-1),ncol(mer))]
            colnames(mer.1) <- paste0('Predicted_',round(confint.level_test*100,0),'CI_',c('Lower','Upper'),'_test_PLM')
            if(exists('O')) {
              O <- cbind(O,mer.1)
            } else {
              O <- mer.1
            }
          }
          if(exists('test.dataset') & Predict_PI_test){
            temp.0 <- data.frame(predict(res_PLM,newdata=test.dataset,interval='prediction',level=confint.level_test))
            temp <- data.frame(temp=rep(NA,nrow(original.dataset)),row.names=rownames(original.dataset))
            mer <- merge(temp,temp.0,by='row.names',sort=F,all=T)
            mer.1 <- mer[order(as.numeric(as.character(mer$Row.names))),c((ncol(mer)-1),ncol(mer))]
            colnames(mer.1) <- paste0('Predicted',round(confint.level_test*100,0),'PI_',c('Lower','Upper'),'_test_PLM')
            if(exists('O')) {
              O <- cbind(O,mer.1)
            } else {
              O <- mer.1
            }
          }
        
          if(Part_index){
            temp <- rep('None',nrow(original.dataset))
            temp[rownames(original.dataset)%in%rownames(dataset)] <- 'Training'
            if(exists('test.dataset')) temp[rownames(original.dataset)%in%rownames(test.dataset)] <- 'Test'
            if(exists('O')) {
              O <- cbind(O,Partition_idx_PLM=temp,stringsAsFactors=F)
            } else {
              O <- data.frame(Partition_idx_PLM=temp,stringsAsFactors=F)
            }
          }
        
          if(exists('O'))	O[is.na(O)] <- ''
          if(exists('warn.res1')|exists('warn.res2')){
            R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file="./test.html")
            if(exists('warn.res1')) R2HTML::HTML(warn.res1,file="./test.html")
            if(exists('warn.res2')) R2HTML::HTML(warn.res2,file="./test.html")
          }
		}
      }
    }

    # Used R packages
    R2HTML::HTML(R2HTML::as.title("Used R Packages"),HR=2,file="./test.html")
    pkg.list <- list(list("Penalized linear regression","cv.glmnet","glmnet"),
                     list("Solution path","plot_glmnet","plotmo"),
                     list("Cross-validation plot","plot.cv.glmnet","glmnet"),
                     list("Linear regression","lm","stats"),
                     list("Confidence interval for regression coefficients","confint","stats"),
                     list("Variance inflation factor (VIF)","vif","rms"),
                     list("ANOVA table",c("anova","Anova"),c("stats","car")),
                     list("Model fitness measurements","residuals, logLik, AIC, BIC","stats"),
                     list("Goodness of fit test","anova","stats"),
                     list("Fitted value","fitted","stats"),
                     list("Predicted value","predict","stats"),
                     list("Confidence & Prediction interval for fitted & predicted value","predict","stats"),
                     list("Unstandardized residual","residuals","stats"),
                     list("Standardized residual","stdres","MASS"),
                     list("Studentized residual","studres","MASS"),
                     list("Cook's distance","cooks.distance","stats"),
                     list("Diagonals of hat matrix","hatvalues","stats"))
    R2HTML::HTML(used.pkg(pkg.list), file="./test.html")

    # Analysis end time
    R2HTML::HTMLhr(file="./test.html")
    R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),
                       ". Rex : Penalized Linear Regression",sep=""),file="./test.html")
    R2HTML::HTMLhr(file="./test.html")
  })

  if(exists('O')){
    return(list(html=html.output,Output=O))
  } else {
    return(html.output)
  }
}

library("R2HTML")
library("glmnet")
REx_ANA_PLOT <- function(w=500,h=500) {
  ## save plot as temp file
  REx.plot.tempfn <<- paste(tempdir(), "\\REx_temp.png", sep="")
  png(filename=REx.plot.tempfn, width=w, height=h)
}

REx_ANA_PLOT_OFF <- function(caption) {
  dev.off()
  load.pkg("markdown")
  ## read temp file as a binary string
  img <- paste(markdown:::.b64EncodeFile(REx.plot.tempfn))
  #R2HTML::HTML(paste("<p align=left><img src='", img, "' /><br /><font class=caption>", caption, "</font></p>", sep=""),file=stdout())
  R2HTML::HTML(paste("<p align=left><img src='", img, "' /><br /><font class=caption>", 
                     caption, "</font></p>", sep=""),file="./test.html")  
}
load.pkg <- function(pkgs) {
  for (i in pkgs) {
    ret <- inst.pkg(i)
    if (ret != 0) {
      if (ret == 2) stop(paste0("package [", i, "] installation error!"))
      else if (ret == 1) {
        if (class(require(i, quietly=TRUE, character.only=TRUE)) == 'try-error')
          stop(paste0("package [", i, "] is not loaded!"))
      }
    }
  }
}
# Function for formatting numbers
Digits <- function(x) {
  x1 <- x
  if(!is.vector(x)){
    # all NA
    allNA <- apply(x,2,function(gg) all(is.na(gg)))
    x <- x[,!allNA,drop=F]
  }
  
  if(length(x)!=0){
    x[abs(x)>=0.0001&!is.na(x)] <- round(x[abs(x)>0.0001&!is.na(x)],4)
    xx = x[abs(x)<0.0001&!is.na(x)&x!=0]
    if(length(xx)>0){
      x[abs(x)<0.0001&!is.na(x)&x!=0] <- paste0(gsub('e','x10<sup>',format(x[abs(x)<0.0001&!is.na(x)&x!=0],scientific=T,digits=4)),'</sup>')
    }
    x[x==0] <- '0'
    if(!is.vector(x)) {
      x1[,!allNA] <- x
    } else {
      x1 <- x
    }
  }
  x1[is.na(x1)] <- ""
  return(x1)
}

inst.pkg <- function(pkg, repos="http://healthstat.snu.ac.kr/CRAN") {
  if (!suppressWarnings(require(pkg, quietly=TRUE, character.only=TRUE))) {
    ret <- try(local.install.packages(pkg, repos=repos), T)
    if (class(ret) == 'try-error') 2
    else 1
  }
  0
}

# Package List
indiv.pkg.info <- function(pkg.info) {
  funs <- strsplit(pkg.info[[2]],",\\s*")[[1]]
  paste0("<li> ",pkg.info[[1]]," : ",
         paste0(paste0(paste0("function ",
                              paste(sapply(funs,function(v)paste0("'<a href=\"https://www.rdocumentation.org/packages/", pkg.info[[3]], "/topics/", v, "\" target='_new'>", v, "</a>'")),collapse=", "),
                              " of R package '"),"<a href=\"https://www.rdocumentation.org/packages/", pkg.info[[3]], "\" target='_new'>", pkg.info[[3]],"</a>'"),collapse=", "))
}
used.pkg <- function(pkg.list) paste(c(sapply(pkg.list,indiv.pkg.info),"<li> All results other than those mentioned above were written with basic functions of R."),collapse=" <br> ")

# Check integer
check.integer <- function(x) x == round(x)
