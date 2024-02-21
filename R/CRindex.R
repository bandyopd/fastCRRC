#' Concordance Index in the Presence of Competing Risks
#'
#' This function computes the concordance index (C-index) for prognostic models in the presence of competing risks.
#' The computation follows the methodology described by Wolbers, M. et al. (2014). The C-index quantifies
#' the agreement between predicted probabilities and observed outcomes, accounting for competing events
#' that preclude the event of interest.
#'
#' @param time Numeric vector of survival times
#' @param status Numeric vector indicating the event type. Status codes should correspond to different event types,
#' including the event of interest and competing events.
#' @param predicted Numeric vector of predicted probabilities or risk scores for the event of interest.
#' @param Cause_int Integer specifying the event of interest among competing risks.
#' Defaults to 1, indicating the first event type as the primary event of interest.
#' Adjust this parameter based on how event types are coded in `status`
#'
#' @return A numeric value representing the concordance index, where higher values indicate better model performance.
#'
#' @export
#'
#' @references
#' Wolbers, M., Blanche, P., Koller, M. T., Witteman, J. C., & Gerds, T. A. (2014). Concordance for prognostic models with competing risks.
#' Biostatistics (Oxford, England), 15(3), 526–539. https://doi.org/10.1093/biostatistics/kxt059


cindexCR <- function(time, status, predicted, Cause_int = 1) {
  time <- ifelse(time == 0, min(time[time > 0] / 2, 1e-5), time)
  reord <- order(time)
  Time_survival <- sort(time)
  Censoring <- ifelse(status == 0, 0, 1)[reord]
  Cause <- ifelse(status == 2, 2, 1)[reord]
  Prediction <- -log(predicted)[reord]
  Time <- max(Time_survival) + 1
  n <- length(Prediction)
  A <- Q <- Nt <- matrix(0, n, n)
  A <- upper.tri(A)
  if (any(duplicated(Time_survival)))  ## If there are ties in survival times
    A <- A[,rank(Time_survival, ties.method = "min")]
  B <- (1 - A) * ((Cause != Cause_int) * Censoring)[col(A)]
  tmp <- order(Prediction)
  if (any(duplicated(Prediction))) ## If there are ties in Prediction
    Q[tmp, tmp] <- lower.tri(A)[rank(Prediction, ties.method = "min"),]
  else ## if no ties
    Q[tmp, tmp] <- lower.tri(A)
  Nt[Time_survival < Time & Cause == Cause_int & Censoring == 1,] <- 1
  return(sum((A + B) * Q * Nt) / sum((A + B) * Nt))
}
