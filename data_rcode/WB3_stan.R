library(rstan)
library(tidyverse)


iowa_rec18_dum <- read.csv("~/Recidivism/data_files/iowa18.csv")
# x is for the survival time variable
# z is for the cure rate variable
iowa18_stan <- list(N = nrow(iowa_rec18_dum), t = iowa_rec18_dum$time_to_recurrence,
                    is_censored = iowa_rec18_dum$is_censored,
                    cured = iowa_rec18_dum$cured, x = cbind(1, iowa_rec18_dum[, -c(3,4,5, 32:38, 41, 42)]),
                    z = iowa_rec18_dum[, -c(3,4,5, 32:42)], K = 28, M = 31)

wb_sample <- stan_model("~/Recidivism/simulation/wb3.stan")
wb_sample_fit <- sampling(wb_sample, data = iowa18_stan, seed = 12345, chains = 4, iter = 2000, init = 0,
                          pars = c("betaU", "betaC", "sigma", "lambda_rp", "lp__", "log_lik"))

saveRDS(wb_sample_fit, "~/Recidivism/dataresult/wb3.rds")


