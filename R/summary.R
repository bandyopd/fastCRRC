is.fastCrrC <- function(x) inherits(x, "fastCrrC")


#' @exportS3Method coef fastCrrC
#' @importFrom stats model.matrix na.omit printCoefmat
coef.fastCrrC <- function(object, ...) {
  coef <- as.numeric(object$coefficient)
  names(coef) <- object$varNames
  coef
}

#' @exportS3Method vcov fastCrrC
vcov.fastCrrC <- function(object, ...) {
  vcov <- object$vcov
  colnames(vcov) <- rownames(vcov) <- object$varNames
  vcov
}

#' @exportS3Method print fastCrrC
print.fastCrrC <- function(x, digits = max(1L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

#' @exportS3Method summary fastCrrC
summary.fastCrrC <- function(object, ...) {
  if (!is.fastCrrC(object))
    stop("Must be fastCrrC class")
  ans <- object["call"]
  est.fastCrrC <- object$coefficient
  if (is.null(object$stderr)) se.fastCrrC <- rep(NaN, length(est.fastCrrC))
  else se.fastCrrC <- object$stderr
  z.fastCrrC <- as.numeric(est.fastCrrC) / as.numeric(se.fastCrrC)
  TAB <- data.frame(estimate = round(drop(est.fastCrrC), 4),
                    std.Error = round(drop(se.fastCrrC), 4),
                    z.value = round(z.fastCrrC, 3),
                    p.value = round(2 * pnorm(-abs(z.fastCrrC)), 4))
  rownames(TAB) <- object$varNames
  out <- list(call = object$call, coefficients = TAB)
  if (object$para$fitter == "fastCrr" & object$para$cindex == TRUE)
    out$cindex <- cindex(object)
  class(out) <- "summary.fastCrrC"
  out
}

#' @exportS3Method print summary.fastCrrC
print.summary.fastCrrC <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("fastCrrC Estimator")
  cat("\n")
  printCoefmat(as.matrix(x$coefficients), P.values = TRUE, has.Pvalue = TRUE)
  cat("\n")
  if (!is.null(x$cindex)) {
    cat("Concordance index: ", x$cindex)
    cat("\n")
  }
}

#' @exportS3Method confint fastCrrC
#' @importFrom stats qnorm
confint.fastCrrC <- function(object, parm, level = 0.95, ...) {
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm))
    parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]
  p <- (1 - level) / 2
  p <- c(p, 1 - p)
  prange <- qnorm(p)
  pct <- paste(format(100 * p, trim = TRUE, scientific = FALSE, digits = 3),"%")
  ci <- array(NA_real_, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  ses <- object$stderr
  parm <- match(parm, pnames)
  ci[] <- cf[parm] + ses[parm] %o% prange
  ci
}

#' Computes concordance
#' @noRd
cindex <- function(object) {
  tt <- object$cumHaz$time
  if (any(names(object$cumHaz) == "Haz"))
    HH <- object$cumHaz$Haz
  else
    HH <- object$cumHaz[,paste0("Haz", object$para$cause)]
  rm <- which(colnames(object$data$mm) == "(Intercept)")
  if (length(rm) > 0) X <- object$data$mm[,-rm]
  else X <- object$data$mm
  cindexCR(object$data$response[,1], object$data$response[,2],
           HH[findInterval(object$data$response[,1], tt, all.inside = TRUE)] * exp(X %*% object$coefficient))
}

