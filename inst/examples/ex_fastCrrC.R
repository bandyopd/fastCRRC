set.seed(1)
dat <- simCR(n = 10, al = .3, ga = .6, clsize = 2:10, summary = TRUE)

fit1 <- fastCrrC(Surv(time, status) ~ Z.1 + Z.2, data = dat, B = 10, fitter = "fastCrr")
fit2 <- fastCrrC(Surv(time, status) ~ Z.1 + Z.2, data = dat, B = 10, fitter = "crrc")

summary(fit1)
summary(fit2)

