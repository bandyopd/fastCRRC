#' Fast implementation of competing risks regression for clustered data
#'
#' Performs the marginal clustered competing risks regression based on the subdistribution hazards model.
#'
#' This is based on the subdistribution hazards model using the two-way linear scan approach to estimate model parameters and bootstrap to compute the variances
#'
#'
#' @param formula a formula expression, of the form `response ~ predictors`.
#' The `response` is a `Surv` object with right censoring.
#' In the case of competing risks, the competing risks indicator can be specified as the status indicator,
#' with 0 indicating right censoring, and 1, 2, ... indicating events other than right censoring.
#' @param data is an optional data.frame in which to interpret the variables occurring in the ‘formula’.
#' @param subset is an optional vector specifying a subset of observations to be used in the fitting process.
#' @param cause is the event of interest; default at 1 when not specified.
#' @param id is an optional vector used to identify the clsuters.
#' When missing, each individual row of `data` is presumed to represent a distinct subject.
#' The length of `id` should be the same as the number of observations.
#' @param fitter is an optional character string specifying the estimating procedures.
#' The available options are "fastCrr" and "crrc" corresponding to
#' the two-way linear scan with bootstrap for variance estimation and
#' the traditional marginal competing risk model using asymptotic variance estimation.
#' We recommend the two-way linear scan for n > 5000 or large sample or cluster and cluster sizes.
#' @param B is a numerical  value specifies the resampling number.
#' When B = 0, only the beta estimate will be displayed.
#' @param bMeth is an optional character string specifying the bootstrap method when `B > 0`.
#' Available methods are the clustered and two-step bootstrap.
#' @param cindex a logical value indicating whether the concordance index should be calculated.
#' Defaults to FALSE.
#' @param multicore a logical value indicating whether parallel processing should be used for computing the variance estimates.
#' Defaults to FALSE.
#' @param mcontrols is an optional list with parameters that pass to fastCrr or crrc, such as maxiter, tolerance, and more.
#' The arguments `gtol`, `maxiter`, and `eps` can be used when `fitter`=`crrc`.
#' The arguments `getBreslowJumps`, `maxiter`, `standardize`, `multicore`, and `mcores` can be only be set when `fitter`=`fastCrr`.
#'
#' @importFrom fastcmprsk fastCrr Crisk
#' @importFrom crrSC crrc
#' @importFrom mcreplicate mc_replicate
#' @import foreach parallel
#'
#' @example inst/examples/ex_fastCrrC.R
#' @export
fastCrrC <- function(formula, data, subset,
                     cause, id, fitter = c("fastCrr", "crrc"),
                     B = 0, bMeth = c("twoStep", "cluStep"),
                     mcontrols = control(gtol=1e-6,maxiter=1000, eps = 1E-6,
                                        getBreslowJumps = TRUE, standardize = TRUE,
                                        cindex = FALSE, multicore = FALSE, mcores=detectCores()-1)) {
  ## scall is to print; scall1 is to check; scall2 is to process
  scall <- scall1 <- scall2 <- match.call()
  fitter <- match.arg(fitter)
  bMeth <- match.arg(bMeth)
  if (missing(cause)) cause <- 1
  ## creating new calls to check and handling competing indicator
  resName <- all.vars(formula[[2]])
  tmp <- as.character(formula[2])
  res1 <- gsub(resName[2], paste0(resName[2], " > 0"), tmp)
  scall1[[2]] <- as.formula(paste0(res1, formula[1], formula[3]))
  res2 <- gsub("Surv", "surv2", tmp)
  res2 <- gsub(")", ", cause)", res2)
  scall2[[2]] <- as.formula(paste0(res2, formula[1], formula[3]))
  ## calling
  mnames <- c("", "formula", "data", "subset", "id")
  cnames <- names(scall2)
  cnames <- cnames[match(mnames, cnames, 0)]
  mcall1 <- scall1[cnames]
  mcall2 <- scall2[cnames]
  mcall1[[1]] <- mcall2[[1]] <- as.name("model.frame")
  m1 <- eval(mcall1, parent.frame())
  m2 <- eval(mcall2, parent.frame())
  ## Check
  if (!is.Surv(m1[[1]]) || ncol(unclass(m1[,1])) > 2 || attr(unclass(m1[,1]), "type") != "right")
    stop("fastCrrC only supports Surv object with right censoring.",  call. = FALSE)
  ## Extract and prepare DF
  mterms <- attr(m2, "terms")
  obj <- unclass(m2[, 1])
  id <- model.extract(m2, id)
  if (is.null(id)) id <- 1:nrow(obj)
  formula[[2]] <- NULL
  if (formula == ~1) stop("No covariates are detected.")
  df <- list(response = obj, id = id, mm = model.matrix(mterms, m2, NULL))
  ## Prepare output
  out <- NULL
  ## setting up the model controls
  model.control = mcontrols
    gtol              <- model.control$gtol
    maxiter           <- model.control$maxiter
    eps               <- model.control$eps
    getBreslowJumps   <- model.control$getBreslowJumps
    standardize       <- model.control$standardize
    cindex            <- model.control$cindex
    multicore         <- model.control$multicore
    mcores            <- model.control$mcores

  ## check if cindex is being used with fastCrr
    if (cindex == TRUE & fitter != "fastCrr") stop("cindex must be used with fitter = fastCrr")

  rm <- which(colnames(df$mm) == "(Intercept)")
  if (fitter == "fastCrr") {
    fit <- fastCrr(Crisk(df$response[,1], df$response[,2]) ~ df$mm - 1, variance = FALSE,
                   eps = 1E-6,max.iter = maxiter, getBreslowJumps = TRUE,standardize = TRUE)
    out$logLik <- fit$logLik
    ind <- findInterval(sort(df$response[,1]), fit$breslowJump$time, all.inside = TRUE)
    ## out$CIF <- data.frame(time = sort(df$response[,1]), cif = 1 - exp(-cumsum(fit$breslowJump$jump))[ind])
    ## attr(out$CIF$cif, "cause") <- cause
    out$cumHaz <- data.frame(time = sort(df$response[,1]), Haz = cumsum(fit$breslowJump$jump)[ind])
  }
  if (fitter == "crrc") {
    if (length(rm) > 0)
      fit <- crrc(df$response[,1], df$response[,2], cov1 = df$mm[,-rm], cluster = df$id)
    else
      fit <- crrc(df$response[,1], df$response[,2], cov1 = df$mm, cluster = df$id)
    out$logLik <- fit$loglik
  }
  out$coefficient <- fit$coef
  out$varNames <- colnames(df$mm)[-rm]
  ## Need bootstrap?
  if (B > 0) {
    performBoot = function(){
      df2 <- bootDt(df, bMeth)
      if (fitter == "fastCrr")
        return(fastCrr(Crisk(df2$response[,1], df2$response[,2]) ~ df2$mm - 1, variance = FALSE,
                       eps = 1E-6,max.iter = maxiter, getBreslowJumps = TRUE,standardize = TRUE)$coef)
      if (fitter == "crrc") {
        if (length(rm) > 0)
          return(crrc(df2$response[,1], df2$response[,2], cov1 = df2$mm[,-rm], cluster = df2$id,
                      gtol=1e-6,maxiter=maxiter)$coef)
        else
          return(crrc(df2$response[,1], df2$response[,2], cov1 = df2$mm, cluster = df2$id,
                      gtol=1e-6,maxiter=maxiter)$coef)
      }
    }
    if (multicore == TRUE)
      bs <- mc_replicate(B, performBoot(),mc.cores=mcores)
    else
      bs <- replicate(B, performBoot())
    out$vcov <- var(t(bs))
    out$stderr <- sqrt(diag(out$vcov))
  }
  out$call <- scall
  out$data <- df
  out$parameters <- list(cause = cause, fitter = fitter, bMeth = bMeth, control = control, multicore = multicore, cindex = cindex)
  out <- out[order(names(out))]
  class(out) <- "fastCrrC"
  return(out)
}


#' Cluster bootstrap methods
#' @noRd
bootDt <- function(df, bMeth){
  ## if there is no cluster
  if (length(df$id) == length(unique(df$id)))
    sampID <- tmp <- sample(df$id, replace = TRUE)
  else {
    dfcnt <- data.frame(id=as.integer(factor(df$id)),count= c(1:length(df$id)))
    sampClu <- sample(unique(df$id), length(unique(df$id)), TRUE)
    tmp <- split(dfcnt, df$id)[sampClu]
    if (bMeth == "twoStep")
      sampID <- do.call(rbind, lapply(tmp, function(d) d[sample(1:nrow(d), nrow(d), T),]))[,2]
    if (bMeth == "cluStep")
      sampID <- as.numeric(do.call(c,lapply(tmp, function(x) x[,2])))
  }
  df$response <- df$response[sampID,]
  df$mm <- df$mm[sampID,]
  df$id <- rep(1:length(tmp), sapply(tmp, length))
  return(df)
}
