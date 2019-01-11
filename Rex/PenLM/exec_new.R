setwd('C://Users//rewki//Documents//WJ_codes//Rex//PenLM')
source('./REx_glmnet.R', encoding = 'UTF-8')
source('./REx_PenLM_new.R')
source('./glmnet.wrapper.1.0.R')
set.seed(1)
n <- 100; p <- 6
x <- matrix(rnorm(n*p), n, p)
epsilon <- rnorm(n,0,0.01)
beta <- c(1,3,0,0,-2,-2)
y = x %*% beta + epsilon
xname <- paste0("x",1:p)
temp <- data.frame(y=y,x=x)
colnames(temp)[2:(p+1)] <- xname

#ex0 <- REx_PenLM(temp,dep_var='y',indep_cat_var=NULL,indep_numeric_var=xname,vars=xname,
#  Penalty='Lasso',TuningPara='Custom',Grid.Num=100, Cross.method='LOOCV',k=10,
#  AccMS='MSE',Part.method='all',Profile=TRUE, Best_model_print=TRUE,ss="III")

ex0 <- REx_PenLM(temp,dep_var='y',
                 indep_cat_var=NULL,
                 indep_numeric_var=xname,vars=xname,
                 Penalty='Lasso', 
                 TuningPara='Custom',Grid.Num=100,
                 Cross.method='KFOLD',k=10,
                 AccMS='MSE',Part.method='percent', train.perc=70, Profile=TRUE,
                 Best_model_print=TRUE, ss="III")

#REx_PenLM <- function(dataset,dep_var,indep_cat_var=NULL,indep_numeric_var=NULL,
#                      vars=NULL,intercept=FALSE,standardize=TRUE,
#                      Penalty=c('Ridge','Lasso','EN'),alpha=0.5,
#                      TuningPara=c('Grid','Custom'),Grid.Num=100,Custom.List=NULL,
#                      Cross.method=c('KFOLD','LOOCV'),k=10,AccMS=c('MSE','MAE'),
#                      Part.method=c('all','percent','variable'),train.perc=70,Part.var=NULL,
#                      Profile=FALSE,Best_model_print=TRUE,CI=TRUE,confint.level=0.95,
#                      VIF=FALSE,ANOVA=TRUE,ss=c('I','II','III'),GOF=FALSE,Plot=FALSE,
#                      Best_model_save=FALSE,Predict_train=FALSE,Predict_CI_train=FALSE,
#                      Predict_PI_train=FALSE,confint.level_train=0.95,Resid=FALSE,
#                      stdResid=FALSE,studResid=FALSE,cook_distance=FALSE,hat_value=FALSE,
#                      Predict_test=FALSE,Predict_CI_test=FALSE,Predict_PI_test=FALSE,confint.level_test=0.95,Part_index=FALSE){
  
#fit <- glmnet_wrapper(y~., data=temp, family = "gaussian", 
#  alpha = 1, nlambda = 100, lambda = NULL, standardize = TRUE, intercept = TRUE) 

#fit2 <- cv_glmnet_wrapper(y~., data=temp, type.measure = "mse", nfolds = 10, 
#  family = "gaussian", alpha = 1, nlambda = 100, standardize = TRUE, intercept = TRUE, opt.lambda="lambda.min")

#predict.cv.glmnet(fit2, newx=, s = c("lambda.1se", "lambda.min"), ...) 

#ex1 <- REx_glmnet(temp, res_var='y', quan_var=xname, 
#         vars=xname, family = "gaussian", standardize = TRUE, intercept=TRUE, 
#         alpha=1, Valid.method = "Noselection", #c("Partition", "Cross", "Noselection"), 
#         Part.method = "percent", #c("percent", "variable"),
#         train.perc = 70, Part.var = NULL,
#         Cross.method= "LOOCV", #c("LOOCV", "KFOLD"),
#         nfolds = 10,
#         TuningPara = "Grid", 
#         Grid.num = 100, Custom.List = NULL,
#         type.measure="mse", opt.lambda="lambda.min", 
#         solutionpath_plot=TRUE, Predict=FALSE, Coef=FALSE)

#ex2 <- REx_glmnet(temp, res_var='y', quan_var=xname,          vars=xname, family = "gaussian", standardize = TRUE, intercept=TRUE,          alpha=1/2, Valid.method = "Cross", #"Partition", #c("Partition", "Cross", "Noselection"), 
#         Part.method = "percent", #c("percent", "variable"),         train.perc = 80, Part.var = NULL,
#         Cross.method= "KFOLD", #c("LOOCV", "KFOLD"),         nfolds = 5,
#         TuningPara = "Grid",          Grid.num = 100, Custom.List = c(3,2.5,2,1,1.5,0),
#         type.measure="mse", opt.lambda="lambda.min",          solutionpath_plot=T, Predict=T, Coef=T)


#aa <- glmnet.wrapper(y~x1+x2, data=temp, family = "gaussian", alpha=1, nlambda = 100, lambda = NULL, standardize = TRUE, intercept = TRUE) 
#write.table(ex1$html, "a.html")
