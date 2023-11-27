library(dplyr)
library(deSolve)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(ggpubr)
source("../R/simulate.R")
source("../R/nllfun.R")

load("../fit/fit_beta_hcov_leaky.rda")
load("../fit/fit_beta_hcov_boosting.rda")

nrevssCDC_ILI <- read.csv("../data/nrevssCDC_ILI.csv") %>%
    mutate(
        WEEKEND=as.Date(WEEKEND)
    )

nrevssCDC_ILI$OC43_ILI[nrevssCDC_ILI$OC43_ILI==0] <- min(nrevssCDC_ILI$OC43_ILI[nrevssCDC_ILI$OC43_ILI!=0])
nrevssCDC_ILI$HKU1_ILI[nrevssCDC_ILI$HKU1_ILI==0] <- min(nrevssCDC_ILI$HKU1_ILI[nrevssCDC_ILI$HKU1_ILI!=0])

ss_leaky <- with(as.list(coef(fit_beta_hcov_leaky)), {
    R0 <- exp(log.R0)
    theta <- plogis(logit.theta)
    phi <- plogis(logit.phi)*52
    delta <- exp(log.delta)
    rho12 <- plogis(logit.rho12)
    rho21 <- plogis(logit.rho21)

    tmp <- simulate_model(R01=R0, R02=R0,
                          theta1=theta, theta2=theta,
                          delta=delta,
                          rho12=rho12,
                          rho21=rho21,
                          model=model_leaky,
                          tvec=seq(0, 29*52+26, by=1)) %>%
        filter(time/52 > 25-27/52)

    scale <- calculate_scale(nrevssCDC_ILI$OC43_ILI, nrevssCDC_ILI$HKU1_ILI,
                             tmp$incidence1, tmp$incidence2)

    data.frame(
        weekend=nrevssCDC_ILI$WEEKEND,
        OC43=tmp$incidence1*exp(scale),
        HKU1=tmp$incidence2*exp(scale)
    )
})

ss_boosting <- with(as.list(coef(fit_beta_hcov_boosting)), {
    R0 <- exp(log.R0)
    theta <- plogis(logit.theta)
    phi <- plogis(logit.phi)*52
    delta <- exp(log.delta)
    rho12 <- plogis(logit.rho12)
    rho21 <- plogis(logit.rho21)

    tmp <- simulate_model(R01=R0, R02=R0,
                          theta1=theta, theta2=theta,
                          delta=delta,
                          rho12=rho12,
                          rho21=rho21,
                          model=model_boosting,
                          tvec=seq(0, 29*52+26, by=1)) %>%
        filter(time/52 > 25-27/52)

    scale <- calculate_scale(nrevssCDC_ILI$OC43_ILI, nrevssCDC_ILI$HKU1_ILI,
                             tmp$incidence1, tmp$incidence2)

    data.frame(
        weekend=nrevssCDC_ILI$WEEKEND,
        OC43=tmp$incidence1*exp(scale),
        HKU1=tmp$incidence2*exp(scale)
    )
})

lfit_leaky <- lm(log(c(nrevssCDC_ILI$OC43_ILI, nrevssCDC_ILI$HKU1_ILI))~1, offset=log(c(ss_leaky$OC43, ss_leaky$HKU1)))
lfit_boosting <- lm(log(c(nrevssCDC_ILI$OC43_ILI, nrevssCDC_ILI$HKU1_ILI))~1, offset=log(c(ss_boosting$OC43, ss_boosting$HKU1)))

AIC(lfit_leaky)
AIC(lfit_boosting)

g1 <- ggplot(nrevssCDC_ILI) +
    geom_point(aes(WEEKEND, OC43_ILI), shape=1) +
    geom_line(data=ss_leaky, aes(weekend, OC43, col="leaky")) +
    geom_line(data=ss_boosting, aes(weekend, OC43, col="boosting")) +
    scale_x_date("Year") +
    scale_y_continuous("Incidence proxy") +
    scale_color_manual(values=c("orange", "black")) +
    ggtitle("A. OC43") +
    theme(
        panel.grid = element_blank(),
        legend.title = element_blank()
    )

g2 <- ggplot(nrevssCDC_ILI) +
    geom_point(aes(WEEKEND, HKU1_ILI), shape=1) +
    geom_line(data=ss_leaky, aes(weekend, HKU1, col="leaky")) +
    geom_line(data=ss_boosting, aes(weekend, HKU1, col="boosting")) +
    scale_x_date("Year") +
    scale_y_continuous("Incidence proxy") +
    scale_color_manual(values=c("orange", "black")) +
    ggtitle("B. HKU1") +
    theme(
        panel.grid = element_blank(),
        legend.title = element_blank()
    )

gcomb <- ggpubr::ggarrange(g1, g2, nrow=1, common.legend = TRUE)

ggsave("figure_beta_hcov.pdf", gcomb, width=8, height=4)
