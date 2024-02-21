#' Control function for parameter and variance estimation
#'
#' This function allows customization of various aspects of the `fitter` argument functions in `fastCrrC`
#' and also allows for multicore processing.
#'
#' @param gtol Numeric: specified for `fitter`=`crrc`; is the threshold for stopping the algorithm based on
#'        a function of the gradient. The iteration stops when this function
#'        of the gradient is less than `gtol`. Default is 1e-6.
#'        algorithm. A value of 0 computes scores and variance at `init`,
#'        but performs no iterations. Default is 10.
#' @param eps Numeric: specified for `fitter`=`fastCrr`; the algorithm stops when the relative change in any
#'        coefficient is less than `eps`. Default is 1E-6.
#' @param maxiter Integer: the maximum number of iterations to achieve convergence. Default is 1000.
#' @param getBreslowJumps Logical: specified for `fitter`=`fastCrr`; whether to output jumps in Breslow estimator
#'        for the cumulative hazard. Default is TRUE.
#' @param standardize Logical: specified for `fitter`=`fastCrr`; whether to standardize the design matrix.
#'        Default is TRUE.
#' @param cindex Logical: whether to compute and output the concordance index for the model
#' @param multicore Logical: whether to use multicore processing. Default is FALSE.
#' @param mcores Integer: number of cores to use for multicore processing. Default is number of cores - 1.
#' @return A list of control parameters configured as specified.
#' @import parallel
#' @export


control <- function(gtol=1e-6,  maxiter = 1000,eps = 1E-6, getBreslowJumps = TRUE, standardize = TRUE,
                    cindex = FALSE, multicore = FALSE, mcores=detectCores()-1)
{
  clist                   <- list()
  clist$gtol              <- gtol
  clist$maxiter           <- maxiter
  clist$eps               <- eps
  clist$getBreslowJumps   <- getBreslowJumps
  clist$standardize       <- standardize
  clist$cindex            <- cindex
  clist$multicore         <- multicore
  clist$mcores            <- mcores

  return(clist)
}


