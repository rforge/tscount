tsglm <- function(ts, model=list(past_obs=NULL, past_mean=NULL, xreg=NULL, external=NULL), link=c("identity", "log"), distr=c("poisson", "nbinom"), ...){
  model_names <- c("past_obs", "past_mean", "xreg", "external")
  stopifnot( #Are the arguments valid?
    all(names(model) %in% model_names)
  )
  model <- model[model_names]
  names(model) <- model_names
  if(is.null(model$xreg)) model$xreg <- matrix(0, nrow=length(ts), ncol=0) else model$xreg <- as.matrix(model$xreg)
  if(length(model$external)==0) model$external <- rep(FALSE, ncol(model$xreg)) else model$external <- as.logical(model$external) #the default value for model$external is FALSE (i.e. an internal covariate effect)
  if(length(model$external)==1) model$external <-  rep(model$external, ncol(model$xreg)) else model$external <- as.logical(model$external) #if only one value for model$external is provided, this is used for all covariates
  if(any(is.na(ts)) || any(is.na(model$xreg))) stop("Cannot make estimation with missing values in time series or regressor")
  stopifnot( #Are the arguments valid?
    model$past_obs%%1==0,
    model$past_mean%%1==0,
    length(ts)==nrow(model$xreg),    
    length(model$external)==ncol(model$xreg)
  )  
  link <- match.arg(link)  
  distr <- match.arg(distr)
  if(link=="log" && distr=="nbinom") stop("Negative binomial distribution is currently only available for the identity link")
  cl <- match.call()
  #Estimating the mean structure:
  if(link=="identity") meanfit <- ingarch.fit(ts=ts, model=model, ...)
  if(link=="log") meanfit <- loglin.fit(ts=ts, model=model, ...)
  if(length(meanfit$coefficients)==0){ #if no final estimation is done, then the function returns a list with less elements and without the class 'tsglm'
    result <- list(init=meanfit$init, call=cl, n_obs=meanfit$n_obs, ts=meanfit$ts, model=meanfit$model, link=link)
    return(result)
  }
  #Estimating the distribution:
  disfit <- distr.fit(meanfit, distr=distr)
  info.matrix_corrected <- if(is.null(meanfit$outerscoreprod)) NULL else apply(as.numeric(1/meanfit$fitted.values + disfit$sigmasq)*meanfit$outerscoreprod, c(2,3), sum)
  if(distr=="poisson") loglik <- sum(dpois(ts, lambda=meanfit$fitted.values, log=TRUE))
  if(distr=="nbinom") loglik <- sum(dnbinom(ts, mu=meanfit$fitted.values, size=disfit$distrcoefs[["size"]], log=TRUE))
  result <- c(
    list(coefficients=meanfit$coefficients, init=meanfit$init, residuals=meanfit$residuals, fitted.values=meanfit$fitted.values, linear.predictors=meanfit$linear.predictors, logLik=loglik, score=meanfit$score, info.matrix=meanfit$info.matrix, info.matrix_corrected=info.matrix_corrected, call=cl, n_obs=meanfit$n_obs, ts=meanfit$ts, model=meanfit$model, link=link), #an extract of the object meanfit with an additional information matrix corrected for the fitted conditional distribution
    disfit
  )
  class(result) <- c("tsglm")
  return(result)
}
