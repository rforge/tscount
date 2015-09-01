predict.tsglm <- function(object, n.ahead=1, newobs=NULL, newxreg=NULL, level=0.95, global=FALSE, type=c("quantiles", "shortest", "onesided"), method=c("conddistr", "bootstrap"), B=1000, estim_error=c("ignore", "normapprox"), ...){
  tsglm.check(object)
  #Link and related functions:
  if(object$link=="identity"){
    trafo <- function(x) x #transformation function
    g_inv <- function(x) x  #inverse of link function
  }
  if(object$link=="log"){
    trafo <- function(x) if(!is.null(x)) log(x+1) else NULL #transformation function
    g_inv <- function(x) exp(x) #inverse of link function
  }
  newxreg <- if(is.null(newxreg)) matrix(0, nrow=n.ahead, ncol=ncol(object$xreg)) else as.matrix(newxreg)  #if no covariates are provided, these are set to zero
  stopifnot(n.ahead>0,
            n.ahead%%1==0,
            ncol(newxreg)==ncol(object$xreg)
  )
  n <- object$n_obs
  model <- object$model
  p <- length(model$past_obs)
  P <- seq(along=numeric(p)) #sequence 1:p if p>0 and NULL otherwise
  q <- length(model$past_mean)
  Q <- seq(along=numeric(q)) #sequence 1:p if p>0 and NULL otherwise
  r <- ncol(object$xreg)
  R <- seq(along=numeric(r)) 
  xreg <- rbind(object$xreg, newxreg)
  if(nrow(xreg) < n+n.ahead) xreg <- rbind(xreg, matrix(xreg[nrow(xreg), ], nrow=n+n.ahead-nrow(xreg), ncol=r, byrow=TRUE)) #if not enough future values of the covariates are given, then the values of the last available time point are used, which is usually NOT sensible!
  new <- rep(NA, n.ahead)
  new[seq(along=newobs)] <- newobs
  ts <- c(object$ts, new)
  if(is.ts(object$ts)) ts <- ts(ts, start=start(object$ts), frequency=frequency(object$ts))
  kappa <- c(rep(NA, n-object$n_eff), object$linear.predictors, rep(NA, n.ahead))
  if(is.ts(object$ts)) kappa <- ts(kappa, start=start(object$ts), frequency=frequency(object$ts))
  for(t in n+(1:n.ahead)){
    kappa[t] <- sum(coef(object)*c(1, trafo(ts[t-model$past_obs]), kappa[t-model$past_mean]-if(r>0){sum((as.numeric(model$external)*coef(object)[1+p+q+R])*t(xreg[t-model$past_mean,]))}else{0}, xreg[t,])) 
    if(is.na(ts[t])) ts[t] <- g_inv(kappa[t]) #unobserved future observations are replaced by their prediction (by the conditional mean)
  }
  if(is.ts(object$ts)){
    pred <- window(g_inv(kappa), start=tsp(object$ts)[2]+1/frequency(object$ts)) #use time series class if input time series has this class
  }else{
    pred <- g_inv(kappa)[n+(1:n.ahead)]
  }
  result <- list(pred=pred)

  #Prediction intervals:
  stopifnot(
    length(level)==1,
    !is.na(level),
    level<1,
    level>=0
  )
  if(level>0){ #do not compute prediction intervals for level==0
    warning_messages <- NULL
    level_local <- if(global){1-(1-level)/n.ahead}else{level} #Bonferroni adjustment of the coverage rate of each individual prediction interval such that a global coverage rate is achieved
    method <- match.arg(method, several.ok=TRUE)
    conddistr_possible <- n.ahead==1 || (length(newobs)>=n.ahead-1 && !any(is.na(newobs[1:(n.ahead-1)]))) #method="conddistr" is only possible if all predictions are 1-step ahead predictions such that the true previous observations are available (this is the case if only one observation is to be predicted or when the first n.ahead-1 observations are given in argument 'newobs')
    if(all(c("conddistr", "bootstrap") %in% method)) method <- if(conddistr_possible){"conddistr"}else{"bootstrap"} #automatic choice of the method, prefer method="conddistr" if possible
    type <- match.arg(type) 
    if(type=="quantiles") a <- (1-level_local)/2
    if(type=="onesided") a <- 1-level_local
    if(type=="shortest"){
      largestdensityinterval <- function(probs, level){ #find shortest interval with given probability, probs are assumed to be the probabilities corresponding to the values seq(along=probs)-1 
        values <- seq(along=probs)-1
        ord <- order(probs, decreasing=TRUE)
        result <- range(values[ord][1:which(cumsum(probs[ord])>=level)[1]])
        return(result)
      }
    }
    estim_error <- match.arg(estim_error)
    if(estim_error=="normapprox" && method!="bootstrap") stop("Accounting for the estimation uncertainty by employing the normal approximation\nof the parameter estimation is only available if argument 'method' is set to\n\"bootstrap\".")
    if(method=="conddistr"){
      if(!conddistr_possible) stop("Computation of prediction intervals with argument 'method' set to \"bootstrap\"\ndoes only work for 1-step-ahead predictions. If argument 'n.ahead' is larger\nthan 1, future observations have to be provided in argument 'newobs'.") 
      if(type %in% c("quantiles", "onesided")){  
        qdistr_mod <- function(p, meanvalue, distr=c("poisson", "nbinom"), distrcoefs, upper=FALSE){ #allows for alternative definition of the quantile
          result <- qdistr(p=p, meanvalue=meanvalue, distr=distr, distrcoefs=distrcoefs)
          if(upper) result <- result + (pdistr(q=result, meanvalue=meanvalue, distr=distr, distrcoefs=distrcoefs)==p) #alternative definition of the quantile
          return(result)
        }
        lower <- if(type=="onesided"){integer(n.ahead)}else{qdistr_mod(a, meanvalue=pred, distr=object$distr, distrcoefs=object$distrcoefs)}
        upper <- qdistr_mod(1-a, meanvalue=pred, distr=object$distr, distrcoefs=object$distrcoefs, upper=TRUE)
        predint <- cbind(lower=lower, upper=upper)
      }
      if(type=="shortest"){
        cutoff <- qdistr(p=1-(1-level_local)/10, meanvalue=max(pred), distr=object$distr, distrcoefs=object$distrcoefs) #very large quantile which is almost certainly larger than the upper bound of the shortest prediction interval
        predint <- t(sapply(pred, function(predval) largestdensityinterval(probs=ddistr(x=0:cutoff, meanvalue=predval, distr=object$distr, distrcoefs=object$distrcoefs), level=level_local)))
      }
      predmed <- qdistr(p=0.5, meanvalue=pred, distr=object$distr, distrcoefs=object$distrcoefs)
      B <- NULL
    }
    if(method=="bootstrap"){
      if(!is.null(newobs)) stop("Computation of prediction intervals with argument 'method' set to \"bootstrap\"\ncan currently not take into account the future observations given in argument\n'newobs'. Either set argument 'method' to \"conddistr\" or do not compute\nprediction intervals at all by setting argument 'level' to zero.")
      stopifnot(
        length(B)==1,
        B%%1==0,
        B>=10
      )
      if(estim_error=="normapprox"){
        param_complete <- simparam(object=object, n=B)
        param <- param_complete$param
        n_invalid <- param_complete$n_invalid
        if(n_invalid > 0){
          bootstrap_message <- paste("In", n_invalid, "cases the bootstrap samples of the regression parameter generated\nfrom the normal approximation was not valid and the respective sample has been\ngenerated again. Note that this might affect the validity of the final result.")
          warning_messages <- c(warning_messages, bootstrap_message)
          warning(bootstrap_message)
        }        
      }else{
        param <- t(replicate(B, coef(object)))        
      }
      simfunc <- function(paramvec, object, xreg, n.ahead){
        object$coefficients <- paramvec
        result <- tsglm.sim(n=n.ahead, xreg=xreg[-(1:object$n_obs), , drop=FALSE], fit=object, n_start=0)$ts
        return(result)
      }
      futureobs <- matrix(apply(param, 1, simfunc, object=object, xreg=xreg, n.ahead=n.ahead), nrow=n.ahead, ncol=nrow(param))          
      if(type %in% c("quantiles", "onesided")){
        quantiles <- t(apply(futureobs, 1, quantile, probs=c(a, 1-a), type=1))
        lower <- if(type=="onesided"){integer(n.ahead)}else{quantiles[, 1]} 
        upper <- quantiles[, 2]
        predint <- cbind(lower=lower, upper=upper)
      } 
      if(type=="shortest") predint <- t(apply(futureobs, 1, function(x) largestdensityinterval(tabulate(x+1)/length(x), level=level_local)))
      predmed <- apply(futureobs, 1, median)
    }
    colnames(predint) <- c("lower", "upper")
    if(is.ts(object$ts)){ #use time series class if input time series has this class
        predint <- ts(predint, start=tsp(object$ts)[2]+1/frequency(object$ts), frequency=frequency(object$ts))
        predmed <- ts(predmed, start=tsp(object$ts)[2]+1/frequency(object$ts), frequency=frequency(object$ts))
    }
    result <- c(result, list(interval=predint, level=level, global=global, type=type, method=method, B=B, estim_error=estim_error, warning_messages=warning_messages, median=predmed))
  }
  return(result)
}
