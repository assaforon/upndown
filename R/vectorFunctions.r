#' Key Probability Vectors of Up-and-Down Designs
#'
#' Asymptotic and progression
#'
#' @details
#' Up-and-Down designs (UDDs) generate random walk behavior, whose theoretical properties can be summarized via a transition probability matrix (TPM). Given the number of doses \eqn{M}, and the value of the cdf \eqn{F} at each dose (i.e., the positive-response probabilities), the specific UDD rules uniquely determine the TPM.
#' 
#' 
#' @note When using the k-in-a-row design, set `matfun = kmatMarg`, not `kmatFull`.
#' 
#' @param cdf vector of positive-response probabilities at the doses
#' @param matfun The function to calculate the TPM. Depends on the specific design; see \code{\link{bcdmat}}. Must


pivec<-function(cdf, matfun, ...)
{
  m=length(cdf)
  imat=matfun(cdf, ...)
  
  vout=cumprod(c(1,imat[cbind(1:(m-1),2:m)]/imat[cbind(2:m,1:(m-1))]))
  vout/sum(vout)
}

######

advancevec<-function(startdose=NULL, cdf, n, matfun,...)
{
  require(expm)
  m=length(cdf)
  ## Starting vector
  vec0=rep(1/m,m)
  if (!is.null(startdose) && startdose %in% 1:m) {
    vec0=rep(0,m)
    vec0[startdose]=1
  }
  ## custom probability vector
  if (length(startdose)==m) vec0=startdose/sum(startdose)
  
  vec0 %*% (matfun(cdf = cdf,...) %^% (n-1))
}

###

cumulpi<-function(startdose=NULL, cdf, n, matfun, average = TRUE, exclude = 0,...)
{
  require(expm)
  m=length(cdf)
  ## Starting vector
  vec0=rep(1/m,m)
  if (!is.null(startdose) && startdose %in% 1:m) {
    vec0=rep(0,m)
    vec0[startdose]=1
  }
  ## custom probability vector
  if (length(startdose)==m) vec0=startdose/sum(startdose)
  
  progmat = matfun(cdf, ...)
  
  ovec=rep(0,m)
  if (is.null(exclude)) exclude=0
  for (a in exclude:(n-1)) ovec=ovec+vec0 %*% (progmat %^% a)
  if(average) ovec=ovec/(n-exclude)
  
  ovec
}





# For `pivec`,  an $M$-length vector with stationary/asymptotic visit frequencies.

