library(dplyr)
library(deSolve)
library(bbmle)
source("../R/simulate.R")
source("../R/nllfun.R")

ss <- simulate_model() %>%
    filter(time/52 > 90)

fit_same_R0_theta <- mle2(nllfun_same_R0_theta,
                          start=list(log.R0=log(1.5),
                                     logit.theta=qlogis(0.1),
                                     log.delta=log(1/52),
                                     logit.rho12=qlogis(0.4),
                                     logit.rho21=qlogis(0.6)),
                          data=list(incidence1=ss$incidence1, incidence2=ss$incidence2),
                          method="Nelder-Mead",
                          skip.hessian = TRUE)

save("fit_same_R0_theta", file="fit_same_R0_theta.rda")
