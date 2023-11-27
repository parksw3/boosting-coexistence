model_boosting <- function(t, y, par) {
    with(as.list(c(y, par)), {
        lambda1 <- R01 * (gamma+mu) * (1 + theta1 * cos(2 * pi * (t-phi1-26)/52)) * (I1 + J1)
        lambda2 <- R02 * (gamma+mu) * (1 + theta2 * cos(2 * pi * (t-phi2-26)/52)) * (I2 + J2)

        dS <- mu - (lambda1 + lambda2 + mu) * S + delta * R1 + delta * R2
        dI1 <- lambda1 * S - (gamma + mu) * I1
        dI2 <- lambda2 * S - (gamma + mu) * I2
        dR1 <- gamma * I1 - (delta + mu) * R1 - lambda2 * R1 + delta * P
        dR2 <- gamma * I2 - (delta + mu) * R2 - lambda1 * R2 + delta * P
        dJ1 <- (1-rho12) * lambda1 * R2 - (gamma + mu) * J1
        dJ2 <- (1-rho21) * lambda2 * R1 - (gamma + mu) * J2
        dP <- gamma * J1 + gamma * J2 - 2 * delta * P - mu * P + rho12 * lambda1 * R2 + rho21 * lambda2 * R1

        N <- S + I1 + I2 + R1 + R2 + J1 + J2 + P

        prevalence1 <- I1 + J1
        prevalence2 <- I2 + J2
        incidence1 <- lambda1 * S + (1-rho12) * lambda1 * R2
        incidence2 <- lambda2 * S + (1-rho21) * lambda2 * R1

        list(c(dS, dI1, dI2, dR1, dR2, dJ1, dJ2, dP),
             N=N, prevalence1=prevalence1, prevalence2=prevalence2,
             incidence1=incidence1, incidence2=incidence2)
    })
}

model_leaky <- function(t, y, par) {
    with(as.list(c(y, par)), {
        lambda1 <- R01 * (gamma+mu) * (1 + theta1 * cos(2 * pi * (t-phi1-26)/52)) * (I1 + J1)
        lambda2 <- R02 * (gamma+mu) * (1 + theta2 * cos(2 * pi * (t-phi2-26)/52)) * (I2 + J2)

        dS <- mu - (lambda1 + lambda2 + mu) * S + delta * R1 + delta * R2
        dI1 <- lambda1 * S - (gamma + mu) * I1
        dI2 <- lambda2 * S - (gamma + mu) * I2
        dR1 <- gamma * I1 - (delta + mu) * R1 - (1-rho21) * lambda2 * R1 + delta * P
        dR2 <- gamma * I2 - (delta + mu) * R2 - (1-rho12) * lambda1 * R2 + delta * P
        dJ1 <- (1-rho12) * lambda1 * R2 - (gamma + mu) * J1
        dJ2 <- (1-rho21) * lambda2 * R1 - (gamma + mu) * J2
        dP <- gamma * J1 + gamma * J2 - 2 * delta * P - mu * P

        N <- S + I1 + I2 + R1 + R2 + J1 + J2 + P

        prevalence1 <- I1 + J1
        prevalence2 <- I2 + J2
        incidence1 <- lambda1 * S + (1-rho12) * lambda1 * R2
        incidence2 <- lambda2 * S + (1-rho21) * lambda2 * R1

        list(c(dS, dI1, dI2, dR1, dR2, dJ1, dJ2, dP),
             N=N, prevalence1=prevalence1, prevalence2=prevalence2,
             incidence1=incidence1, incidence2=incidence2)
    })
}

simulate_model <- function(R01=1.5,
                           R02=1.5,
                           theta1=0.1,
                           theta2=0.1,
                           phi1=26,
                           phi2=26,
                           gamma=1,
                           delta=1/52,
                           rho12=0.4,
                           rho21=0.6,
                           mu=1/80/52,
                           model=model_boosting,
                           tvec=seq(0, 100*52, by=1)) {
    par <- c(R01=R01, R02=R02, theta1=theta1, theta2=theta2, phi1=phi1, phi2=phi2,
             gamma=gamma, delta=delta, rho12=rho12, rho21=rho21,
             mu=mu)

    yini <- c(S=1-2e-6, I1=1e-6, I2=1e-6, R1=0, R2=0, J1=0, J2=0, P=0)

    out <- as.data.frame(ode(yini, tvec, model, par))

    # plot(out$time/52, out$incidence1, type="l", xlim=c(80, 100))
    # lines(out$time/52, out$incidence2, type="l", col=2)
    out
}
