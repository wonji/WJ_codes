#setwd('C://Users//rewki//Documents//WJ_codes//Rex//PenLM')
setwd("~/WJ_codes/Rex/PenLM/")
#source('./REx_glmnet.R', encoding = 'UTF-8')
source('./REx_PenLM_new.R')
source('../확인된오류/prior_ft.txt')

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
                 TuningPara='Grid',Grid.Num=100,
                 Cross.method='KFOLD',k=10,
                 AccMS='MSE',Part.method='percent', train.perc=70, Profile=TRUE,
                 Best_model_print=TRUE, ss="III")

ex0 <- REx_PenLM(temp,dep_var='y',
                 indep_cat_var=NULL,
                 indep_numeric_var=xname,vars=xname,noint=T,
                 Penalty='Lasso', 
                 TuningPara='Grid',Grid.Num=100,
                 Cross.method='LOOCV',
                 AccMS='MSE',Part.method='percent', train.perc=70, Profile=TRUE,
                 Best_model_print=TRUE, ss="III")


ex1 <- REx_PenLM(temp,dep_var='y',
                 indep_cat_var=NULL,
                 indep_numeric_var=xname,vars=xname,
                 Penalty='EN',alpha=0.6, 
                 TuningPara='Grid',Grid.Num=200,
                 Cross.method='KFOLD',k=10,
                 AccMS='MSE',Part.method='percent', train.perc=70, Profile=TRUE,
                 Best_model_print=TRUE, ss="III")

ex2 <- REx_PenLM(temp,dep_var='y',
                 indep_cat_var=NULL,
                 indep_numeric_var=c('x3','x4'),vars=c('x3','x4'),noint=T,
                 Penalty='Lasso', 
                 TuningPara='Grid',Grid.Num=100,
                 Cross.method='LOOCV',
                 AccMS='MSE',Part.method='percent', train.perc=70, Profile=TRUE,
                 Best_model_print=TRUE, ss="III")


dataset <- temp;
dep_var='y';
indep_cat_var=NULL;
indep_numeric_var=xname;vars=xname;noint=F;
Penalty='Lasso'; 
TuningPara='Grid';Grid.Num=100;
Cross.method='LOOCV';
AccMS='MSE';Part.method='percent'; train.perc=70; Profile=TRUE;
Best_model_print=TRUE; ss="III"
