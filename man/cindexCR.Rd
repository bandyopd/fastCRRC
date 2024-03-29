% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CRindex.R
\name{cindexCR}
\alias{cindexCR}
\title{Concordance Index in the Presence of Competing Risks}
\usage{
cindexCR(time, status, predicted, Cause_int = 1)
}
\arguments{
\item{time}{Numeric vector of survival times}

\item{status}{Numeric vector indicating the event type. Status codes should correspond to different event types,
including the event of interest and competing events.}

\item{predicted}{Numeric vector of predicted probabilities or risk scores for the event of interest.}

\item{Cause_int}{Integer specifying the event of interest among competing risks.
Defaults to 1, indicating the first event type as the primary event of interest.
Adjust this parameter based on how event types are coded in `status`}
}
\value{
A numeric value representing the concordance index, where higher values indicate better model performance.
}
\description{
This function computes the concordance index (C-index) for prognostic models in the presence of competing risks.
The computation follows the methodology described by Wolbers, M. et al. (2014). The C-index quantifies
the agreement between predicted probabilities and observed outcomes, accounting for competing events
that preclude the event of interest.
}
\references{
Wolbers, M., Blanche, P., Koller, M. T., Witteman, J. C., & Gerds, T. A. (2014). Concordance for prognostic models with competing risks.
Biostatistics (Oxford, England), 15(3), 526–539. https://doi.org/10.1093/biostatistics/kxt059
}
