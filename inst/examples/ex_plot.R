set.seed(1)
dat <- simCR(n = 100, al = .3, ga = .6, cen = .4, clsize = 2:5, summary = TRUE)
dat$status <- I(dat$status > 0)* sample(1:5, nrow(dat), TRUE)

foo <- fastCrrC(Surv(time, status) ~ Z.1 + Z.2, data = dat)
str(foo)
plot(foo)
plot(foo, newdata = data.frame(Z.1 = .5, Z.2 = 1))

## Plotting CIFs that involves other causes will automatically update the fit
plot(foo, type = "all")
plot(foo, type = "all", newdata = data.frame(Z.1 = .5, Z.2 = 1))
str(foo)

## Stacking CIF
plot(foo, type = "stack")
plot(foo, type = "stack", newdata = data.frame(Z.1 = .5, Z.2 = 1))
