

#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	 
#' 
#' @example 
#' 
#' @export
#'
#' @return 
#'  
#' @param x1,y1 numeric vectors: sequence of administered doses and responses for group 1. Must be same length, and with `y1` coded `TRUE/FALSE` or 0/1.
#' @param x2,y2 numeric vectors: sequence of administered doses and responses for group 2. Must be same length, and with `y2` coded `TRUE/FALSE` or 0/1.
#' @param target The target response rate for which target dose estimate is requested. Must be a single number in \eqn{(0,1).}
#' @param conf The desired confidence level for the confidence interval (CI). Default \eqn{95\%}. Note the difference with single-group estimates where the recommended level is \eqn{90\%}. See more under "Value".

armDiffCI <- function(x1, y1, x2, y2, target, conf = 0.95, ...)
{
  est1 = udest(x=x1, y=y1, target=target, conf=conf, ...)
  est2 = udest(x=x2, y=y2, target=target, conf=conf, ...)
  
  wid1 = c(est1$point - est1[1,3], est1[1,4] - est1$point)
  wid2 = c(est2$point - est2[1,3], est2[1,4] - est2$point)
  pointdiff = est1$point - est2$point
  
  data.frame(point = pointdiff, lower = pointdiff - sqrt(wid1[1]^2 + wid2[2]^2), 
             upper = pointdiff + sqrt(wid1[2]^2 + wid2[1]^2) )
  
}
