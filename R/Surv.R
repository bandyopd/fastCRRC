#' \code{Surv} function imported from \code{survival}
#'
#' This function is imported from the \code{survival} package. See
#' \code{\link[survival]{Surv}}.
#'
#' @importFrom survival Surv
#' @name export_Surv
#' @aliases Surv
#' @export Surv
NULL

#' \code{is.Surv} function imported from \code{survival}
#'
#' This function is imported from the \code{survival} package. See
#' \code{\link[survival]{is.Surv}}.
#'
#' @importFrom survival is.Surv
#' @name export_is.Surv
#' @aliases is.Surv
#' @export is.Surv
NULL

#' A new Surv()-like function to accommodate competing risk
#'
#' This convert the "status" vector to
#' 0 = right censor; 1 = event of interest; 2 = competing risks
#'
#' @param y is a numerical vector containing observed survival time
#' @param d is an indicator vector; could accommodate 0, 1, 2, ...
#' @param k is cause to indicate event of interest
#'
#' @return when d has more than 0 and 1,
#' this returns a matrix with 3 columns, time, status in 0-1-2, and the raw statsu
#'
#' @noRd
surv2 <- function(y, d, k) {
  if (length(unique(d)) <= 2) return(Surv(y, d))
  tmp <- Surv(y)
  tmp[,2] <- ifelse(d == k, 1, 2) * (d > 0)
  attr(tmp, "raw") <- d
  return(tmp)
}
