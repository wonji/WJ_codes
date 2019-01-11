#library("R2HTML")
#library("car")
#library("lme4")
#library("lattice")
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
REx_glmnet <- function(dataset, res_var, quan_var=NULL, vars=NULL, 
  family = c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian"),
  standardize = TRUE, intercept=FALSE, weights=NULL, alpha=1, 
  Valid.method = c("Partition", "Cross", "Noselection"), 
  Part.method = c("percent", "variable"),
  train.perc = 70, Part.var=NULL,
  Cross.method=c("LOOCV", "KFOLD"),
  nfolds = 10, 
  TuningPara = c("Grid", "Custom"), 
  Grid.num = 100, Custom.List = NULL,
  type.measure="mse", opt.lambda="lambda.min", 
  solutionpath_plot=FALSE, Predict=FALSE, Coef=FALSE){

  ###################
  #### Required packages : R2HTML, glmnet ####
  ###################
  #load.pkg(c("R2HTML", "glmnet"))
  n <- nrow(dataset)

  html.output <- capture.output({
    # Title
    R2HTML::HTML(R2HTML::as.title("Penalized Regression Analysis"),
      HR=1, file="./test.html", append=FALSE)
    
    warn.msg <- NULL
    
    ## Warnings
    if(!is.numeric(dataset[, res_var])) 
      warn.msg <- c(warn.msg, '\a Error : Response variable should be numeric. 
        Analysis has been stopped.')
    
    if(is.null(vars)){
      Vars <- c()
    }else{
      Vars <- unique(unlist(strsplit(vars,":")))
    }

    #qual_var <- qual_var[qual_var %in% Vars]
    quan_var <- c(quan_var[quan_var %in% Vars],res_var)
    
    if(!is.null(quan_var)){
      is.nom <- sapply(quan_var,function(i) !is.numeric(dataset[,i]))
      if(any(is.nom)){
        warn.msg <- c(warn.msg, paste0("\a Warning : The type of variable '",
          paste(quan_var[is.nom],collapse=', '), 
          "' is not numeric but selected as 
          the quantitative variable. It was excluded from the analysis."))
        
        # explanatory variables
        vars <- vars[!vars %in% ec]
        if(length(vars)==0) vars <- NULL
      }
    }
    #if(!is.null(qual_var)) {
    #  is.num <- sapply(qual_var,function(i) is.numeric(dataset[,i]))
    #  if(any(is.num))
    #    warn.msg5 <- paste0("\a Warning : The type of variable '", paste(qual_var[is.num],collapse=', '),"' is continuous but selected as the qualitative variable. It was coreced into character.")
    #  for(i in qual_var) dataset[,i] <- factor(dataset[,i])
    #}
    if(Valid.method == "Partition"){ # partition-validation
      if(Part.method == "percent"){
        if(train.perc < 0 | train.perc > 100){
          warn.msg <- c(warn.msg, '\a Error : Partition percent should be input appropriately!')
        }
        trainid <- sample(1:n, floor(train.perc/100 *n))
        validid <- setdiff(1:n, trainid)
      }else if(Part.method =="variable"){
        trainid <- (1:n)[var == 1]
        validid <- setdiff(1:n, trainid)
      }
    }else if(Valid.method == "Noselection"){
      trainid <- 1:n
      validid <- NULL
    }

    if(TuningPara=="Grid"){
      if(Grid.num <=2)
        warn.msg <- c(warn.msg, "number of grids should be greater than 2")
      nlambda = Grid.num
      lambda = NULL
    }else if(TuningPara == "Custom"){
      if(any(!is.numeric(Custom.List)) | is.null(Custom.List))
        stop("Custom option need the list of candidate lambdas(decreasing order).")
      lambda = Custom.List
      nlambda = NULL
    }

    if(!is.null(warn.msg)){
      R2HTML::HTML(R2HTML::as.title("Warnings"), HR=2, file="./test.html")
      R2HTML::HTML(warn.msg, file="./test.html")
    }else{
      ## Model formula
      fe <- paste(vars, collapse=' + ')
      if(fe=="") 
        fe <- NULL
      
      form.0 <- paste0(res_var,' ~ ', paste(fe, collapse=' + '))
      
      if(!intercept) 
        form.0 <- paste(form.0, -1)
      
      form.1 <- as.formula(form.0)
      
      ## Data processing
      # id var
      #dataset[,id_var] <- factor(dataset[,id_var])
      if(is.null(weights)){
        weights <- rep(1, NROW(dataset))
      }
      
      ## Model Fitting
      if(Valid.method =="Partition" | Valid.method == "Noselection"){
        fit <- suppressWarnings(try(glmnet_wrapper(form.1, data=dataset, 
                 family=family, trainid=trainid, validid=validid, weights, 
                 alpha=alpha, nlambda=nlambda, lambda=lambda, 
                 standardize=standardize, intercept=intercept)))
      }else if(Valid.method =="Cross"){
        if(Cross.method == "LOOCV"){
          foldid <- 1:n
          fit <- suppressWarnings(try(cv_glmnet_wrapper(form.1, data=dataset, 
                   weights, lambda=lambda, 
                   type.measure=type.measure, nfolds=nfolds, foldid,
                   family=family, alpha=alpha, nlambda=nlambda, 
                   standardize=standardize, intercept=intercept)))
        }else if(Cross.method =="KFOLD"){
          fit <- suppressWarnings(try(cv_glmnet_wrapper(form.1, data=dataset, 
                   weights, lambda=lambda, 
                   type.measure=type.measure, nfolds=nfolds, 
                   family=family, alpha=alpha, nlambda=nlambda, standardize=standardize, 
                   intercept=intercept, opt.lambda=opt.lambda)))
        }
      }

      ## Data Structure
      R2HTML::HTML(R2HTML::as.title("Data Structure"), HR=2, file="./test.html")
      total.var <- NCOL(dataset)
      used.var <- ifelse(is.null(vars), 0, length(unique(unlist(strsplit(vars,":")))))
      DS <- matrix(c('Number of observations','Number of total variables',
        'Number of used variables',nrow(dataset),total.var,used.var+1), ncol=2)
      R2HTML::HTML(DS,file="./test.html", innerBorder = 1, align="left")
      
      ## Varibale list
      R2HTML::HTML(R2HTML::as.title("Variable List"), HR=2, file="./test.html")
      if(is.null(vars)){
        Vars <- c()
      }else{
        Vars <- unique(unlist(strsplit(vars,":")))
      }
      #qual <- qual_var[qual_var%in%Vars]
      quan <- c(quan_var[quan_var%in%Vars],res_var)
      varlist <- matrix(c('Quantitative variable', paste(quan,collapse=', ')), ncol=2, byrow=T)
      #if(length(qual)!=0) varlist <- rbind(varlist,c('Qualitative variable',paste(qual,collapse=', ')))
      
      R2HTML::HTML(varlist, file="./test.html", innerBorder = 1, align="left")

      ## Analysis Description
      R2HTML::HTML(R2HTML::as.title("Analysis Description"), HR=2, file="./test.html")
      AD <- matrix(c('Response variables', res_var,
                     'Explanatory variable', paste(vars,collapse=', '),
                     'Intercept included', intercept),ncol=2,byrow=T)

      if(Valid.method=="Partition"){
        varmethod <- "Partition"
        varmethod2 <- ifelse(Part.method=="percent", 
                             paste0(train.perc, "\\%"), 
                             paste("partition by ", variable))
      }else if(Valid.method=="Cross"){
        varmethod <- "Cross-Validation"
        varmethod2 <- ifelse(Cross.method=="LOOCV", "leave-one-out", 
                             paste0(nfolds, "-fold CV"))
      }else if(Valid.method=="Noselection"){
        varmethod <- "No model selection"
        varmethod2 <- ""
      }
      AD <- rbind(AD, matrix(c('Distribution of response variable', family,
                               'Link function', "identiy", 
                               'Partition', varmethod,
                               'method', varmethod2), ncol=2, byrow=T))
      R2HTML::HTML(AD, file="./test.html", innerBorder = 1, align="left")
      
      ### Model Fitting
      R2HTML::HTML(R2HTML::as.title(paste0("Penalized Regression Model Fit by ", varmethod)),
        HR=2, file="./test.html")

      ## Check convergence
      CC <- capture.output(summary(fit))
      is.conv <- grep("Model failed to converge", CC)
      
      if(length(is.conv)==0){
        ## Time variable
        # Estimated effect
        R2HTML::HTML(R2HTML::as.title("Estimated effects"), HR=3, file="./test.html")
        xx <- as.matrix(coef(fit, s=fit$lambda))
        colnames(xx) <- paste0("lam",1:length(fit$lambda))
        CE <- data.frame(xx[-2,])
        #CE <- cbind(CE,pt(abs(CE[,3]),df=1,lower.tail=F))
        #colnames(CE) <- c('Estimate','SE','T-value','P-value')
        #CE[,ncol(CE)] <- format(CE[,ncol(CE)],scientific=T,digits=4)
        
        R2HTML::HTML(CE,file="./test.html", innerBorder=1, align="left", digits=3)
        
        # Measure of GOF according to lambda
        R2HTML::HTML(R2HTML::as.title("Measures of Goodness of Fit"), HR=3, file="./test.html")
        if(Valid.method!="Cross"){
          R2HTML::HTML(R2HTML::as.title(deparse(fit$call)), HR=7, file="./test.html")
          GOF <- as.data.frame(cbind(Df = fit$df, `%Dev` = signif(fit$dev.ratio, fit$digits), 
                                  Lambda = signif(fit$lambda, fit$digits)))
          R2HTML::HTML(GOF, file="./test.html", innerBorder = 1, align="left", digits=4,row.names=F)
        }else if(Valid.method =="Cross"){
          fit2 <- fit$glmnet.fit
          R2HTML::HTML(R2HTML::as.title(deparse(fit2$call)), HR=7, file="./test.html")
          GOF <- as.data.frame(cbind(Df = fit2$df, `%Dev` = signif(fit2$dev.ratio, 4), 
                                     Lambda = signif(fit2$lambda, 4)))
          R2HTML::HTML(GOF, file="./test.html", innerBorder = 1, align="left", digits=4,row.names=F)
        }
        # Solution Path Plot
        if(solutionpath_plot & length(fit$lambda) > 1){
          if(Valid.method!="Cross"){
            R2HTML::HTML(R2HTML::as.title("Solution Path Plot"),HR=3,file="./test.html")
            REx_ANA_PLOT()
            par(mfrow=c(2,2))
            glmnet::plot.glmnet(fit, xvar="norm", label=TRUE)
            glmnet::plot.glmnet(fit, xvar="lambda", label=TRUE)
            glmnet::plot.glmnet(fit, xvar="dev", label=TRUE)
            REx_ANA_PLOT_OFF('')
          }else{
            R2HTML::HTML(R2HTML::as.title("Cross-Validation Plot"),HR=3,file="./test.html")
            REx_ANA_PLOT()
            glmnet::plot.cv.glmnet(fit)
            REx_ANA_PLOT_OFF('')
          }
        }
        #[1] "a0"        "beta"      "df"        "dim"       "lambda"    "dev.ratio" "nulldev"   "npasses"   "jerr"     
        #[10] "offset"    "call"      "nobs"
        if(Predict) {
          R2HTML::HTML(names(fit),file="./test.html", innerBorder=1, align="left", digits=3)
          #fx <- predict.glmnet(fit, x, s=fit$lambda, type="response")
          prediction <- data.frame(prediction=fit$fx)
        }
        if(Coef) {
          xx <- as.matrix(coef(fit, s=fit$lambda))
          colnames(xx) <- paste0("lam", 1:length(fit$lambda))
          coefficient = data.frame(xx[-2,])
        }
      } else {
        conv.msg <- paste('\a', CC[is.conv])
        R2HTML::HTML(conv.msg, file="./test.html")
      }
    }
    R2HTML::HTMLhr(file="./test.html")
    R2HTML::HTML(paste0("Analysis is finished at ", 
      Sys.time(), ". REx : Penalized Regression"), "./test.html")
    R2HTML::HTMLhr(file="./test.html")
  })
  R2HTML::HTMLhr(Predict, file="./test.html")
  if(Predict & Coef){
    print("############# output 'O'") 
    return(list(html=html.output,Output=list(prediction, coefficient)))
  }else if(Predict){
    print("html.output")
    return(list(html=html.output,Output=list(prediction)))
  }else if(Coef){
    print("html.output")
    return(list(html=html.output,Output=list(coefficient)))
  }else{
    print("html.output")
    return(list(html=html.output))
  }
}

###### Examples
#birth <- read.csv("http://healthstat.snu.ac.kr/homepage/files//Intro/birth.csv",head=T)
#birth$time <- rep(1:25,each=20)
#birth$obs <- rep(1:20,25)

#ex1 <- REx_glmnet_regression(birth, res_var='bweight', time_var='time', id_var='obs', qual_var=c('lowbw','preterm','hyp','sex'), quan_var=c('gestwks','matage'), vars=c('lowbw','preterm','hyp','sex','gestwks','matage','gestwks:matage'),intercept=FALSE, time=c('1','2','3','4','5'), time_code=c(0,1,2,3,4), rand_int=TRUE, rand_time=TRUE, rand_corr=FALSE, var_method='ML', resid_plot=TRUE, profile_plot=TRUE, Predict=TRUE, Resid=TRUE);
#ex2 <- REx_glmnet_regression(birth, res_var='bweight', time_var='time', id_var='obs', qual_var=c('lowbw','preterm','hyp','sex'), quan_var=c('gestwks','matage'), vars=NULL,intercept=FALSE, time=c('1','2','3','4','5'), time_code=c(0,1,2,3,4), rand_int=TRUE, rand_time=TRUE, rand_corr=FALSE, var_method='ML', resid_plot=TRUE, profile_plot=TRUE, Predict=TRUE, Resid=TRUE);
#ex3 <- REx_glmnet_regression(birth, res_var='bweight', time_var='time', id_var='obs', qual_var=c('lowbw','preterm','hyp','sex'), quan_var=c('gestwks','matage'), vars=NULL,intercept=TRUE, time=c('1','2','3','4','5'), time_code=c(0,1,2,3,4), rand_int=TRUE, rand_time=TRUE, rand_corr=FALSE, var_method='ML', resid_plot=TRUE, profile_plot=TRUE, Predict=TRUE, Resid=TRUE);

#write.table(ex1$html, "a.html")


