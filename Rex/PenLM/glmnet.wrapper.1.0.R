###
### glmnet wrapper function provides formula
###

glmnet_wrapper <- function(formula, data,
  family = c("gaussian", "binomial"), 
  alpha = 1, nlambda = 100, lambda = NULL, standardize = TRUE, intercept = TRUE, ...) 
{
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric") #w <- as.vector(model.weights(mf))

  ###
  ###
  #print(is.empty.model(mt))
  if(is.empty.model(mt)){
    x <- NULL
    z <- list(coefficients = if (is.matrix(y))  matrix(, 0, 3) else numeric(), 
          residuals = y, fitted.values = 0*y, weights = w, rank = 0L, 
          df.residual = if (!is.null(w)) sum(w != 0) else if (is.matrix(y)) nrow(y) else length(y))
  }else {
    x <- model.matrix(mt, mf)
    #R2HTML::HTML("HERE", file="./test.html")
    z <- glmnet::glmnet(x, y, family = family, offset = NULL, 
           alpha = alpha, nlambda = nlambda, lambda = lambda, 
           standardize = standardize, intercept = intercept, ...)
    
    # model selection via validation data set
#    if(!is.null(validid)){
#      pred <- predict.glmnet(z, x[validid,], s=z$lambda, type="response")
#      R2HTML::HTML(c(dim(pred), length(z$lambda)), file="./test.html")
#      res <- pred-matrix(rep(y[validid], times=length(z$lambda)), nrow=length(validid))
#      mse <- apply(res^2, 2, mean)
#      R2HTML::HTML(mse, file="./test.html")
        
#      z <- glmnet::glmnet(x, y, family = family, weights = w, offset = NULL, 
#             alpha = alpha, lambda = z$lambda[which.min(mse)], standardize = standardize, 
#             intercept = intercept, ...)
#    }
  }
  z$mt <- mt
  z$mf <- mf
  z$fx <- predict.glmnet(z, x, s=z$lambda, type="response")
  z$digits = max(3, getOption("digits") - 3)
  z
}
predict_glmnet_wrapper <- function(object, newx, s=NULL, type = c("link", "response", 
                                                                 "coefficients", "nonzero", "class"))
{
  if(class(newx) !="matrix")
    newx <- as.matrix(newx)
  
  predict.glmnet(object, newx, s, type)
}
predict_cv.glmnet_wrapper <- function(object, newx, s = c("lambda.1se", "lambda.min"), ...) 
{
  if(class(newx) !="matrix")
    newx <- as.matrix(newx)
  
  predict.cv.glmnet(object, newx, s="lambda.min", ...) 
}
###
### cv.glmnet wrapper function provides formula
###
cv_glmnet_wrapper <- function(formula, data, 
  weights, offset=NULL, lambda = NULL, 
  type.measure = c("mse", "mae"), nfolds = 10, 
  family = c("gaussian", "binomial"), 
  alpha = 1, nlambda = 100, standardize = TRUE, intercept = TRUE, opt.lambda="lambda.min")
{
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")

  ###
  ###
  ###
  if(is.empty.model(mt)){
    x <- NULL
    z <- list(coefficients = if (is.matrix(y))  matrix(, 0, 3) else numeric(), 
              residuals = y, fitted.values = 0*y, weights = w, rank = 0L, 
              df.residual = if (!is.null(w)) sum(w != 0) else if (is.matrix(y)) nrow(y) else length(y))
    
    if (!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }
  }else {
    x <- model.matrix(mt, mf)
    z <- glmnet::cv.glmnet(x, y, lambda = lambda, 
           type.measure = type.measure, nfolds = nfolds, 
           family = "gaussian", alpha = alpha, nlambda = nlambda, 
           standardize = standardize, intercept = intercept) 
  }  
  z$mt <- mt
  z$mf <- mf
  z$fx <- predict.cv.glmnet(z, x, s=opt.lambda, type="response")
  z$digits = max(3, getOption("digits") - 3)
  z
}


#set.seed(1)
#n <- 100
#x1 <- rnorm(n,0,1)
#x2 <- rnorm(n,0,1)
#epsilon <- rnorm(n,0,0.001)
#y = 2*x1 + 3*x2 +epsilon
#temp <- data.frame(y=y,x1=x1,x2=x2)

#aa <- glmnet_reg_wrapper(y~x1+x2, data=temp, family = "gaussian", alpha=1, nlambda = 100, lambda = NULL, standardize = TRUE, intercept = TRUE) 
#aa$fx  







#