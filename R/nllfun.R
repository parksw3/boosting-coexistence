calculate_scale <- function(data_incidence1, data_incidence2,
                            pred_incidence1, pred_incidence2) {
    y <- log(c(data_incidence1, data_incidence2))
    x <- log(c(pred_incidence1, pred_incidence2))

    lfit <- lm(y~1, offset=x)

    coef(lfit)[[1]]
}

nllfun_same_R0_theta_phi <- function(log.R0=log(1.5),
                                 logit.theta=qlogis(0.1),
                                 logit.phi=qlogis(0.5),
                                 log.delta=log(1/52),
                                 logit.rho12=qlogis(0.4),
                                 logit.rho21=qlogis(0.6),
                                 incidence1, incidence2,
                                 tmin=90,
                                 tvec=seq(0, 100*52, by=1),
                                 model=model_leaky) {
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
                          model=model,
                          tvec=tvec) %>%
        filter(time/52 > tmin)

    scale <- calculate_scale(incidence1, incidence2,
                             tmp$incidence1, tmp$incidence2)

    y <- log(c(incidence1, incidence2))
    x <- log(c(tmp$incidence1, tmp$incidence2))

    nll <- sum((y-(scale + x))^2)

    print(nll)

    nll
}
