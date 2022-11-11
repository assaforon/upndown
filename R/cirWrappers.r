#' Centered-Isotonic-Regression (CIR) Estimate for the Up-and-Down Target Dose
#' 
#' Centered Isotonic Regression (CIR) is an extension of isotonic regression (IR), substantially improving upon IR's estimation performance in the dose-response and dose-finding contexts (Oron and Flournoy 2017, Flournoy and Oron 2020). CIR is the recommended method for estimating up-and-down targets. 
#' 
#'  CIR and related methods are available in the `cir` package. The `udest()` function in the present package provides a convenient wrapper for `cir::quickInverse()`, with arguments already set to the appropriate values for estimating the target dose after an up-and-down experiment. The function also returns a confidence interval as default.
#' 
#' **WARNING!** You should not estimate target doses too far removed from the design's actual balance point (definitely no further than 0.1, e.g., estimating the 33rd percentile for a design whose balance point is the median). As Flournoy and Oron (2020) explain, observed response rates are biased away from the balance point. Even though `udest()` performs the rudimentary bias correction described in that article, practically speaking this correction's role is mostly to expand the confidence intervals in response to the bias. It cannot guarantee to provide reliable off-balance-point estimates.
#'

#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	  

#' 
#' @export
#'
#' @return A one-row data frame with 4 variables: `target`, `point` (the point estimate), `lowerXYconf, upperXYconf` (the confidence bounds, with `XY` standing for the percents, default `90`).
#'  
#' @param x numeric vector: sequence of administered doses, treatments, stimuli, etc.
#' @param y numeric vector: sequence of observed responses. Must be same length as `x`, and must be coded `TRUE/FALSE` or 0/1.
#' @param target The target response rate for which target dose estimate is requested. Must be a single number in \eqn{(0,1).}
#' @param balancePt In case the design's inherent balance point differs somewhat from `target`, specify it here to improve estimation accuracy. See Details for further explanation. Otherwise, this argument defaults to be equal to `target`.
#' @param conf The desired confidence level for the confidence interval. Default \eqn{90\%.} We do not recommend increasing to \eqn{95\%} unless you have \eqn{\sim 100} or more observations. 

#' @references 
#'  - Oron AP, Flournoy N.  Centered Isotonic Regression: Point and Interval Estimation for Dose-Response Studies. *Statistics in Biopharmaceutical Research* 2017; 9, 258-267. [Author's public version available on arxiv.org.](https://arxiv.org/pdf/1701.05964.pdf)
#'  - Flournoy N, Oron AP. Bias Induced by Adaptive Dose-Finding Designs. *Journal of Applied Statistics* 2020; 47, 2431-2442.  
#'  - Oron AP, Souter MJ, Flournoy N. Understanding Research Methods: Up-and-down Designs for Dose-finding. *Anesthesiology* 2022; 137:137–50. [See in particular the open-access Supplement.](https://cdn-links.lww.com/permalink/aln/c/aln_2022_05_25_oron_aln-d-21-01101_sdc1.pdf)
#'  
#'  
udest <- function(x, y, target, balancePt = target, conf = 0.9)
{
  require(cir)
  
if(length(target) > 1) stop("Experiment should have a single target.\n")
checkTarget(target)
if(length(x) != length(y)) stop('x, y, must have same length.\n')
checkDose(x)
checkResponse(y)

if(balancePt[1]!=target)
{
  if(length(balancPt) > 1) stop("Experiment can only have a   single balance point.\n")
  checkTarget(balancePt, tname='balancePt')
}

quickInverse(x=x, y=y, target=target, starget = balancePt,
             conf=conf, adaptiveShrink=TRUE,
              adaptiveCurve = (target != 0.5) )

}