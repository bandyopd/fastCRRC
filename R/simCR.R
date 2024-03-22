#' Function needs for stable distribution
#' @references
#' Feller, W. (1991). An introduction to probability theory and its applications.
#' Volume 2 (Vol. 81). John Wiley & Sons. p. 583.
#'
#' @param a is the scale parameter for the stable distribution
#'
#' @noRd
g <- function(a) {
  iu <- complex(real = 0, imaginary = 1)
  return(abs(1 - iu * tan(pi * a / 2)) ^ (-1 / a))
}

#' Function to generate simulated competing risk data
#'
#' The function \code{simCR()} generates data from a hypothetical multicenter.
#' To generate correlated failure times, the subdistribution hazard with positive frailty is used to
#' maintain the PH assumption on sub-distribution.
#' The dataframe returns the cluster IDs, event status, event time, censoring time
#'
#' @param n is an integer specifying the number of clusters. Default at 100.
#' @param al is a numerical value generated from a positive stable distribution for cause 1 to control the within cluster correlation. Default is set at 0.3.
#' @param ga is a numerical value generated from a positive stable distribution for cause 2 to control the within cluster correlation. Default is set at 0.3.
#' @param rho is a numerical value specifying the degree of the within cluster dependence.
#' @param cen is a numerical value specifying the censoring rate. Only 0\%, 20\%, 40\%, 60\%, 65\%, and 80\% are currently allowed.
#' @param clsize an integer specifying the range of cluster sizes.
#' Each cluster size is randomly generated from a discrete uniform distribution with range [a, b].
#' where a and b are set to achieve the desired censoring rates of 0\%, 20\%, 40\%, 60\%, 65\%, and 80\%.
#' @param summary is an indicator indicating whether to print data summary
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom stabledist rstable
#' @import stats
#'
#' @example inst/examples/ex_simu.R
simCR <- function(n = 100, al = .3, ga = .3, rho = 0,
                  cen = 0.65, clsize = 2:50, summary = FALSE) {
  call <- match.call()
  b1 <- c(0.5, -0.5)
  b2 <- rep(1,length(b1))
  p <- length(b1)
  tau1 <- as.matrix(b1 / al)
  tau2 <- as.matrix(b2 / ga)
  dat <- data.frame(id = 1:n,
                    vi = rstable(n, alpha = al, beta = 1, gamma = g(al), delta = 0, pm = 1),
                    hi = rstable(n, alpha = ga, beta = 1, gamma = g(ga), delta = 0, pm = 1))
  m <- sample(clsize, n, TRUE)
  nn <- sum(m)
  dat <- dat[rep(dat$id, m),]
  rownames(dat) <- NULL
  ## Covariates at individual and cluster levels
  dat <- data.frame(dat, Z = mvrnorm(nrow(dat), rep(0, p), diag(1 - rho, p) + rho))
  dat$Z.2 <- unlist(tapply(dat$Z.2, dat$id, function(x) rep(mean(x), length(x))))
  Y <- as.matrix(subset(dat, select = paste0("Z.", 1:p)))
  dat$tZ <- as.numeric(exp(Y %*% tau1))
  dat$events <- 2 - with(dat, rbinom(sum(m), 1, 1 - exp(-vi * tZ)))
  tmp <- 1 - exp(- exp(Y %*% b1))
  uu <- runif(nn, 1 - tmp, 1)
  t2Z <- as.numeric(Y %*% tau2)
  b1Z <- as.numeric(Y %*% -b1)
  dat$Ti <- with(dat, ifelse(events == 2,
                             rexp(sum(m), hi * exp(t2Z)),
                             -log(1 - (-log(uu) * exp(b1Z)) ^ (1 / al))))
  if (!(cen %in% c(0, .2, .4, .6, .65, .8)))
    stop("Available censoring rates are 0%, 20%, 40%, 60%, 65%, and 80%.")
  if (cen == 0) dat$cen <- Inf
  if (cen == 0.2) dat$cen <- runif(nn, 0, 200)
  if (cen == 0.4) dat$cen <- runif(nn, 0, 0.7) #runif(nn, 0, .7)
  if (cen == 0.6) dat$cen <- runif(nn, 0, 5e-2)
  if (cen == 0.65) dat$cen <- runif(nn, 0, 3e-2)
  if (cen == 0.8) dat$cen <- runif(nn, 0, 2e-3)
  ## if(ics=="pois"){
  ##  dat$Yi <- pmin(dat$Ti, dat$cen)+rep(mui,m)
  ## }  else
  dat$time <- pmin(dat$Ti, dat$cen)
  dat$status <- 1 * (dat$Ti <= dat$cen)
  dat$status <- dat$status * dat$events
  if (summary) {
    cat("Call: \n")
    print(call)
    cat("\n")
    cat("Number of clusters generated:                   ", n, "\n")
    cat("Average number of observations per cluster:    ", round(nrow(dat) / n, 3), "\n")
    cat("Overall observed censoring rate:               ", round(mean(dat$status == 0), 3), "\n")
    cat("\n")
    cat("Before censoring:\n")
    cat("Proportion of event 1:                         ", round(mean(dat$events == 1), 3), "\n")
    cat("Proportion of event 2:                         ", round(mean(dat$events == 2), 3), "\n")
    cat("After censoring:\n")
    cat("Proportion of event 1:                         ", round(mean(dat$status == 1), 3), "\n")
    cat("Proportion of event 2:                         ", round(mean(dat$status == 2), 3), "\n")
    cat("Average proportion of event 1 per cluster:     ",
        round(mean(aggregate(I(status == 1) ~ id, data = dat, mean)[,2]), 3), "\n")
    cat("Average proportion of event 2 per cluster:     ",
        round(mean(aggregate(I(status == 2) ~ id, data = dat, mean)[,2]), 3), "\n")
    cat("\n")
  }
  data=subset(dat, select=c(id,time,status,Z.1,Z.2))
  return(data)
}
