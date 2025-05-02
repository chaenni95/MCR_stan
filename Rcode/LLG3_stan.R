library(rstan)
library(tidyverse)


iowa_rec18_dum <- read.csv("/home/chy20005/Recidivism/Recidivism_rstan/iowa18.csv")
# x is for the survival time variable
# z is for the cure rate variable
iowa18_stan <- list(N = nrow(iowa_rec18_dum), t = iowa_rec18_dum$time_to_recurrence, 
                    is_censored = iowa_rec18_dum$is_censored, 
                    cured = iowa_rec18_dum$cured, x = cbind(1, iowa_rec18_dum[, -c(3,4,5, 32:38, 41, 42)]), 
                    z = iowa_rec18_dum[, -c(3,4,5, 32:42)], K = 28, M = 31)


# llog_sample <- stan_model("/home/chy20005/Recidivism/Recidivism_rstan/aft_llogistic_rpl.stan")
# #llog_sample <- stan_model("/Users/chaenni/Documents/UCONN/rstan/llog_sample.stan")
# llog_sample_fit <- sampling(llog_sample, data = iowa18_stan, seed = 12345, chains = 4, iter = 2000, 
#                             init = 0, pars = c("betaU", "betaC", "sigma","delta", "lp__", "log_lik", "y_rep","pp"))
# 
# saveRDS(llog_sample_fit, "/home/chy20005/res/res_stan_fix_0528/llog_sample3.rds")

llg_sample <- stan_model("/home/chy20005/Recidivism/simulation/LLG3.stan")
llg_sample_fit <- sampling(ln_sample, data = iowa18_stan, seed = 12345, chains = 4, iter = 2000, init = 0, 
                           pars = c("betaU", "betaC", "sigma", "delta", "lp__", "log_lik"))

saveRDS(llg_sample_fit, "/home/chy20005/Recidivism/simulation/res/llg3.rds")
