library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
load("../simulate_bifurcation/simulate_bifurcation_rho.rda")

meanrho <- 0.5
simulate_bifurcation_rho2 <- simulate_bifurcation_rho %>%
    mutate(
        scenario=paste0("rho[12]==", meanrho-rhodiff, "~~rho[21]==", meanrho+rhodiff)
    )

simulate_bifurcation_rho3 <- simulate_bifurcation_rho2 %>%
    filter(time %% 52 ==0)

g1 <- ggplot(simulate_bifurcation_rho3) +
    geom_point(aes(R0, incidence1, col="strain 1"), size=0.5) +
    geom_point(aes(R0, incidence2, col="strain 2"), size=0.5) +
    scale_x_continuous("Basic reproduction number") +
    scale_y_log10("Incidenxe") +
    scale_color_manual(values=c("black", "orange")) +
    facet_grid(model~scenario, labeller = label_parsed) +
    theme(
        panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        strip.background = element_blank()
    )

ggsave("figure_bifurcation.pdf", g1, width=8, height=4)

simulate_bifurcation_rho_min <- simulate_bifurcation_rho2 %>%
    group_by(model, scenario, R0) %>%
    summarize(
        min1=min(prevalence1),
        min2=min(prevalence2)
    )

ggplot(simulate_bifurcation_rho_min) +
    geom_point(aes(R0, min1, col="strain 1"), size=0.5) +
    geom_point(aes(R0, min2, col="strain 2"), size=0.5) +
    scale_x_continuous("Basic reproduction number") +
    scale_y_log10("Epidemic trough") +
    scale_color_manual(values=c("black", "orange")) +
    facet_grid(model~scenario, labeller = label_parsed) +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        strip.background = element_blank()
    )
