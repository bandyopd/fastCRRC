
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Introduction to the `fastCRRC` R package

### Introduction

`fastCRRC` implements a fast and efficient clustered competing risks
(CCR) regression model for large correlated clustered data. This
regression model assesses the marginal effects of covariates on the
respective cumulative incidence functions while accommodating the
correlation induced due to clustering. Parameter estimation utilizes the
‘forward-backward scan’ algorithm to drastically reduce the computation
time for large clustered competing risk datasets. Variance estimation is
efficiently handled using a clustered bootstrap approach.

### CCR Model

Let $i=1,2,\ldots,n$ denote the clusters where each cluster has
$k=1, \ldots, M_i$ observations. For the $k$th member in cluster $i$,
let $T_{ik}$ and $C_{ik}$ denote the failure and right censoring times,
respectively. Then for right censored data, the corresponding observable
random variables become $X_{ik}= \min(T_{ik}, C_{ik})$ and
$\xi_{ik} = I(T_{ik} \leq C_{ik})\epsilon_{ik}$ where
$\epsilon_{ik} \in \{1, \ldots,l\}$ denotes the cause of failure,
$I(\cdot)$ is the indicator function, and
$Z_{ik} = (Z_{ik1},...,Z_{ikp})^{\top}$ is a $p \times 1$ vector of
time-independent covariates for individual $k$ in cluster $i$.

The marginal cumulative incidence function conditional on the covariates
is defined as
$F_1(t;Z_{ik}) = P(T_{ik} \leq t, \, \epsilon_{ik} = 1|Z_{ik})$ and the
marginal subdistribution hazard is modeled as

$$\lambda_1(t; Z_{ik}) = \lambda_{10}(t) \exp(\beta^\top_{0}Z_{ik}),$$

where $\lambda_{10}(\cdot)$ is a completely unspecified baseline
subdistribution hazard function and $\beta_0$ is a $p \times 1$ vector
of unknown regression parameters.

Vignettes is available
[here](https://html-preview.github.io/?url=https://github.com/edemprd/fastCRRC/blob/main/vignettes/fastCRRC.html)

### Installation

You can install the development version of `fastCRRC`from
[GitHub](https://github.com/edemprd/fastCRRC) with:

``` r
# install.packages("devtools")
devtools::install_github("edemprd/fastCRRC")
```

### Example

This is a basic example which shows you how to simulate data and fit the
CCR model using the package.

``` r
library(fastCRRC)
set.seed(1)
dat <- simCR(n = 10, al = .3, ga = .6, clsize = 2:10, summary = TRUE)
#> Call: 
#> simCR(n = 10, al = 0.3, ga = 0.6, clsize = 2:10, summary = TRUE)
#> 
#> Number of clusters generated:                    10 
#> Average number of observations per cluster:     7.9 
#> Overall observed censoring rate:                0.722 
#> 
#> Before censoring:
#> Proportion of event 1:                          0.772 
#> Proportion of event 2:                          0.228 
#> Before censoring:
#> Proportion of event 1:                          0.253 
#> Proportion of event 2:                          0.025 
#> Average proportion of event 1 per cluster:      0.235 
#> Average proportion of event 2 per cluster:      0.029

fit1 <- fastCrrC(Surv(time, status) ~ Z.1 + Z.2, data = dat, B = 100, fitter = "fastCrr")

summary(fit1)
#> Call:
#> fastCrrC(formula = Surv(time, status) ~ Z.1 + Z.2, data = dat, 
#>     fitter = "fastCrr", B = 100)
#> 
#> fastCrrC Estimator
#>     estimate std.Error z.value p.value
#> Z.1   0.4390    0.2709   1.621  0.1050
#> Z.2   0.0181    1.2406   0.015  0.9884
```

### References

- Zhou, B., Fine, J., Latouche, A., & Labopin, M. (2012). Competing
  risks regression for clustered data. Biostatistics, 13(3), 371-383.

- Kawaguchi, E. S., Shen, J. I., Suchard, M. A., & Li, G. (2021).
  Scalable algorithms for large competing risks data. Journal of
  Computational and Graphical Statistics, 30(3), 685-693.

<!--You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this.-->
