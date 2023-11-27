library(dplyr)
library(deSolve)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)
source("../R/simulate.R")
source("../R/nllfun.R")

load("../fit/fit_same_R0_theta.rda")

ss_boosting <- simulate_model() %>%
    filter(time/52 > 90)

cc <- coef(fit_same_R0_theta)

ss_leaky <- simulate_model(R01=exp(cc[["log.R0"]]),
                           R02=exp(cc[["log.R0"]]),
                           theta1=plogis(cc[["logit.theta"]]),
                           theta2=plogis(cc[["logit.theta"]]),
                           delta=exp(cc[["log.delta"]]),
                           rho12=plogis(cc[["logit.rho12"]]),
                           rho21=plogis(cc[["logit.rho21"]]),
                           model=model_leaky) %>%
    filter(time/52 > 90)

scale <- calculate_scale(ss_boosting$incidence1, ss_boosting$incidence2,
                         ss_leaky$incidence1, ss_leaky$incidence2)

g1 <- ggplot(ss_boosting) +
    geom_point(aes(time/52, incidence1)) +
    geom_line(data=ss_leaky, aes(time/52, exp(scale)*incidence1)) +
    scale_x_continuous("Year", expand=c(0, 0), limits=c(90, 100), breaks=c(90, 92, 94, 96, 98)) +
    scale_y_continuous("Incidence", expand=c(0, 0), limits=c(0, 0.024)) +
    ggtitle("A. Strain 1") +
    theme(
        panel.grid = element_blank()
    )

g2 <- ggplot(ss_boosting) +
    geom_point(aes(time/52, incidence2)) +
    geom_line(data=ss_leaky, aes(time/52, exp(scale)*incidence2)) +
    scale_x_continuous("Year", expand=c(0, 0), limits=c(90, 100), breaks=c(90, 92, 94, 96, 98)) +
    scale_y_continuous("Incidence", expand=c(0, 0), limits=c(0, 0.024)) +
    ggtitle("B. Strain 2") +
    theme(
        panel.grid = element_blank()
    )

gcomb <- ggarrange(g1, g2, nrow=1, draw=FALSE)

ggsave("figure_fit_same_R0_theta.pdf", gcomb, width=8, height=4)
