#' Centered-Isotonic-Regression (CIR) Estimate for the Up-and-Down Target Dose
#' 
#' CIR is the recommended method for estimating up-and-down (UD) targets. This is a UD-adapted convenience wrapper for a function from the `cir` package.
#' 
#' Centered Isotonic Regression (CIR) is an extension of isotonic regression (IR), substantially improving upon IR's estimation performance in the dose-response and dose-finding contexts (Oron and Flournoy 2017, Flournoy and Oron 2020). CIR and related methods are available in the `cir` package. The `udest()` function in the present package provides a convenient wrapper for `cir::quickInverse()`, with arguments already set to the appropriate values for estimating the target dose after an up-and-down experiment.
#' 
#' WARNING! You should not estimate


#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	  


#' @inheritParams reversmean
#' 
#' @export
#' 
#' @param target The target response rate for which target dose estimate is requested. Must be a single number in $(0,1).$
#' @param balancePt In case the design's inherent balance point differs somewhat from `target`, specify it here to improve estimation accuracy. See Details for further explanation. Otherwise, this argument defaults to be equal to `target`
#' @param conf The desired confidence level for the confidence interval. Default $90\%,$ and we recommend not increasing to $95\%$ unless you have about $\sim 100$ or more observations. 


#' https://cdn-links.lww.com/permalink/aln/c/aln_2022_05_25_oron_aln-d-21-01101_sdc1.pdf

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