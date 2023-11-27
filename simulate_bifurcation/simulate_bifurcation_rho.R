library(dplyr)
library(deSolve)
source("../R/simulate.R")

meanrho <- 0.5
R0vec <- seq(1.5, 4, length.out=101)
rhodiffvec <- c(0.1, 0.2, 0.3)

tvec <- seq(0, 200*52, by=1)

paramdata <- expand.grid(R0vec, rhodiffvec)

reslist <- vector('list', nrow(paramdata))

for (i in 1:nrow(paramdata)) {
    print(i)
    pp <- paramdata[i,]

    R0 <- pp[[1]]
    rhodiff <- pp[[2]]

    ss_boosting <- simulate_model(R01=R0,
                                  R02=R0,
                                  rho12=meanrho-rhodiff,
                                  rho21=meanrho+rhodiff,
                                  model=model_boosting,
                                  tvec=tvec) %>%
        mutate(model="boosting")

    ss_leaky <- simulate_model(R01=R0,
                               R02=R0,
                               rho12=meanrho-rhodiff,
                               rho21=meanrho+rhodiff,
                               model=model_leaky,
                               tvec=tvec) %>%
        mutate(model="leaky")

    ss_boosting_filter <- ss_boosting %>%
        filter(time/52 >= 190)

    ss_leaky_filter <- ss_leaky %>%
        filter(time/52 >= 190)

    reslist[[i]] <- bind_rows(
        ss_boosting_filter,
        ss_leaky_filter
    ) %>%
        mutate(
            R0=R0,
            rhodiff=rhodiff
        )
}

simulate_bifurcation_rho <- reslist %>%
    bind_rows

save("simulate_bifurcation_rho", file="simulate_bifurcation_rho.rda")
