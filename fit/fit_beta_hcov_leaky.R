library(dplyr)
library(deSolve)
library(bbmle)
source("../R/simulate.R")
source("../R/nllfun.R")

nrevssCDC_ILI <- read.csv("../data/nrevssCDC_ILI.csv")

nrevssCDC_ILI$OC43_ILI[nrevssCDC_ILI$OC43_ILI==0] <- min(nrevssCDC_ILI$OC43_ILI[nrevssCDC_ILI$OC43_ILI!=0])
nrevssCDC_ILI$HKU1_ILI[nrevssCDC_ILI$HKU1_ILI==0] <- min(nrevssCDC_ILI$HKU1_ILI[nrevssCDC_ILI$HKU1_ILI!=0])

fit_beta_hcov_leaky <- mle2(nllfun_same_R0_theta_phi,
                      start=list(log.R0=log(1.7),
                                 log.delta=log(1/40),
                                 logit.rho21=qlogis(0.7),
                                 logit.rho12=qlogis(0.5),
                                 logit.theta=qlogis(0.17),
                                 logit.phi=qlogis(0.5)),
                      data=list(
                          incidence1=nrevssCDC_ILI$OC43_ILI,
                          incidence2=nrevssCDC_ILI$HKU1_ILI,
                          tvec=seq(0, 29*52+26, by=1),
                          tmin=25-27/52
                      ),
                      method="Nelder-Mead",
                      skip.hessian = TRUE,
                      control=list(maxit=1e5))

save("fit_beta_hcov_leaky", file="fit_beta_hcov_leaky.rda")
