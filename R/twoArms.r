#' Utility for significance testing between 2 Up-and-Down Arms
#' 
#' A t-test like approach to comparing two Up-and-Down arms, returning the point estimate and confidence interval of the difference.
#' 
#' If the two arms were run independently and used the same Up-and-Down design, then the CI for the difference can be found by *"Pythagorean addition"* (root of sum-of-squares) of the opposite parts of the two CIs (right with left and vice versa).
#' 
#' Note that if the arms are part of a multiple-arm experiment, then multiple hypothesis testing should be accounted for separately, preferably via False Discovery Rate adjustment.

#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	 
#' 
#' 
#' @export
#'
#' @return 
#'  - If `full = FALSE` (default), a data frame with the `point` difference (Arm 2 minus Arm 1),`lowerYY` and `upperYY` bounds, with `YY` being the confidence level.
#'  - If `full = TRUE`, a list containing that data frame, as well as each arm's analogous estimates.
#'  
#' @param x1,y1 numeric vectors: sequence of administered doses and responses for group 1. Must be same length, and with `y1` coded `TRUE/FALSE` or 0/1.
#' @param x2,y2 numeric vectors: sequence of administered doses and responses for group 2. Must be same length, and with `y2` coded `TRUE/FALSE` or 0/1.
#' @param target The target response rate for which target dose estimate is requested. Must be a single number in \eqn{(0,1).}
#' @param conf The desired confidence level for the confidence interval (CI). Default \eqn{95\%}. Note the difference with single-group estimates where the recommended level is \eqn{90\%}. See more under "Value".
#' @param full Should the function also return each arm's estimates (`TRUE`), or only the point+interval estimates of the difference (`FALSE`, default)?
#' @param ...	Other arguments passed on to \code{\link{udest}}. 


armDiffCI <- function(x1, y1, x2, y2, target, conf = 0.95, full = FALSE, ...)
{

# Validation is handled by the functions called below
  
  est1 = udest(x=x1, y=y1, target=target, conf=conf, ...)
  est2 = udest(x=x2, y=y2, target=target, conf=conf, ...)
  
  wid1 = c(est1$point - est1[1,3], est1[1,4] - est1$point)
  wid2 = c(est2$point - est2[1,3], est2[1,4] - est2$point)
  pointdiff = est2$point - est1$point
  
  tmp = data.frame(point = pointdiff, lower = pointdiff - sqrt(wid1[2]^2 + wid2[1]^2), 
              upper = pointdiff + sqrt(wid1[1]^2 + wid2[2]^2))
  names(tmp)[2:3] = paste(c('lower', 'upper'), round(100 * conf), 'conf', sep = '')

  if(full) return(list(arm1 = est1, arm2 = est2, difference = tmp))
  
  tmp
  
}
