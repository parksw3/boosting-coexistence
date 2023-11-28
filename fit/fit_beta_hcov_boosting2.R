library(dplyr)
library(deSolve)
library(bbmle)
source("../R/simulate.R")
source("../R/nllfun.R")

load("fit_beta_hcov_boosting.rda")

nrevssCDC_ILI <- read.csv("../data/nrevssCDC_ILI.csv")

nrevssCDC_ILI$OC43_ILI[nrevssCDC_ILI$OC43_ILI==0] <- min(nrevssCDC_ILI$OC43_ILI[nrevssCDC_ILI$OC43_ILI!=0])
nrevssCDC_ILI$HKU1_ILI[nrevssCDC_ILI$HKU1_ILI==0] <- min(nrevssCDC_ILI$HKU1_ILI[nrevssCDC_ILI$HKU1_ILI!=0])

fit_beta_hcov_boosting <- mle2(nllfun_same_R0_theta_phi,
                      start=list(log.R0=coef(fit_beta_hcov_boosting)[["log.R0"]],
                                 log.delta=coef(fit_beta_hcov_boosting)[["log.delta"]],
                                 logit.rho21=0.9,
                                 logit.rho12=coef(fit_beta_hcov_boosting)[["logit.rho12"]],
                                 logit.theta=coef(fit_beta_hcov_boosting)[["logit.theta"]],
                                 logit.phi=coef(fit_beta_hcov_boosting)[["logit.phi"]]),
                      data=list(
                          incidence1=nrevssCDC_ILI$OC43_ILI,
                          incidence2=nrevssCDC_ILI$HKU1_ILI,
                          tvec=seq(0, 29*52+26, by=1),
                          tmin=25-27/52,
                          model=model_boosting
                      ),
                      method="BFGS",
                      skip.hessian = TRUE,
                      control=list(maxit=1e5))

save("fit_beta_hcov_boosting", file="fit_beta_hcov_boosting.rda")
