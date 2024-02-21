#' Draw CIF
#'
#' @param x is an "\code{fastCrrC}" object
#' @param type types of plots to produce.
#' The available options are: one = plot CIF for the specified cause, all = plot CIF for all causes, stack = stack CIFs for all causes
#' \code{type = "all"} and \code{type = "stack"} will call \code{extend.fastCrrC()}
#' @param newdata is a data frame for an optional new data for CIF. Default to a vector of zero,
#' in which case the \code{plot()} function plots the baseline CIF.
#'
#' @importFrom ggplot2 ggplot aes geom_step geom_area xlab ylab ggtitle scale_fill_discrete
#' @import stats
#' @export
#'
#' @method plot fastCrrC
#' @example inst/examples/ex_plot.R
plot.fastCrrC <- function(x, type = c("one", "all", "stack"), newdata) {
  if (!is.fastCrrC(x)) stop("Must be fastCrrC class")
  if (x$para$fitter != "fastCrr")
    stop('CIF is only available when fitter = "fastCrr"')
  type <- match.arg(type)
  causes <- sort(unique(attr(x$data$response, "raw")))
  if (missing(newdata)) X <- matrix(0, 1, length(x$coefficient))
  else {
    tmp <- formula(x$call)
    tmp[[2]] <- NULL
    ## X <- model.matrix(as.formula(tmp), data = data.frame(x$data$mm))
    X <- model.matrix(as.formula(tmp), data = data.frame(newdata))
  }
  if (type == "one") {
    risk <- exp(sum(colMeans(X)[-1] * x$coefficient))
    if (any(names(x$cumHaz) == "Haz"))
      p <- ggplot(x$cumHaz, aes(time, 1 - exp(-Haz * risk))) + geom_step()
    else
      p <- ggplot(x$cumHaz, aes(time, eval(
        parse(text = paste0("1 - exp(-Haz", x$para$cause, " * risk)"))))) + geom_step()
  } else {
    if (length(unique(attr(x$data$response, "raw"))) > ncol(x$cumHaz)) {
      x <- extend.fastCrrC(x)
      assign(as.character(match.call()[[2]]), x, pos = 1)
    }
    risk <- exp(colSums(colMeans(X)[-1] * x$coefList))
    pstate <- subset(x$cumHaz, select = -time)
    pstate <- data.frame(1 - exp(-t(t(pstate) * risk)))
    if (type == "all") {
      pstate2 <- data.frame(time = rep(x$cumHaz$time, length(causes) - 1),
                            value = unlist(pstate),
                            event = factor(rep(paste0("cif", causes[-1]), each = nrow(x$cumHaz))))
      p <- ggplot(pstate2, aes(time, value, color = event)) + geom_step()
    }
    if (type == "stack") {
      pstate$s0 <- 1 - rowSums(pstate)
      pstate$s0 <- ifelse(pstate$s0 < 0, 0, pstate$s0)
      pstate <- pstate / rowSums(pstate)
      pstate2 <- data.frame(time = rep(x$cumHaz$time, length(causes)),
                            value = unlist(pstate),
                            event = factor(rep(c(paste0("cif", causes[-1]), "(S0)"), each = nrow(x$cumHaz))))
      p <- ggplot(pstate2, aes(time, value, fill = event)) + geom_area()
    }
  }
  p + ylab("Probability") + xlab("Time") + ggtitle("Stacked CIF plot") +
    scale_fill_discrete(name = "", labels = c("Overall survival", paste0("CIF cause ", causes[-1])))
}

#' Extend the "\code{fastCrrC}" object to accmodiate different causes
#'
#' This is an internal function and is called when a stacked CIF plot is called in the plot method.
#'
#' @noRd
extend.fastCrrC <- function(x) {
  if (!is.fastCrrC(x)) stop("Must be fastCrrC class")
  cause <- x$parameters$cause
  causes <- sort(unique(attr(x$data$response, "raw")))
  extra <- setdiff(causes, c(0, cause))
  extraFits <- lapply(extra, function(k) update(x, cause = k))
  cumHazList <- lapply(extraFits, function(y) y$cumHaz$Haz)
  cumHazList <- c(cumHazList, list(x$cumHaz$Haz))
  coefList <- lapply(extraFits, function(y) y$coefficient)
  coefList <- c(coefList, list(x$coefficient))
  names(cumHazList) <- names(coefList) <- c(extra, cause)
  cumHazList <- cumHazList[order(names(cumHazList))]
  coefList <- coefList[order(names(coefList))]
  pstate <- do.call(cbind, cumHazList)
  bstate <- do.call(cbind, coefList)
  pstate <- as.data.frame(pstate)
  colnames(pstate) <- paste0("Haz", causes[-1])
  colnames(bstate) <- paste0("cause", causes[-1])
  x$cumHaz <- cbind(time = x$cumHaz$time, pstate)
  x$coefList <- data.frame(bstate)
  x <- x[order(names(x))]
  class(x) <- "fastCrrC"
  invisible(x)
}


