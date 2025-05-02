library(rstan)
library(tidyverse)



iowa_rec18_dum <- read.csv("/home/chy20005/Recidivism/Recidivism_rstan/iowa18.csv")
# x is for the survival time variable
# z is for the cure rate variable
iowa18_stan <- list(N = nrow(iowa_rec18_dum), t = iowa_rec18_dum$time_to_recurrence, 
                    is_censored = iowa_rec18_dum$is_censored, 
                    cured = iowa_rec18_dum$cured, x = cbind(1, iowa_rec18_dum[, -c(3,4,5, 32:38, 41, 42)]), 
                    z = iowa_rec18_dum[, -c(3,4,5, 32:42)], K = 28, M = 31)


#ln_sample <- stan_model("/Users/chaenni/Documents/UCONN/rstan/ln_sample.stan")
# ln_sample <- stan_model("/home/chy20005/Recidivism/Recidivism_rstan/aft_lognormal_pl.stan")
# ln_sample_fit <- sampling(ln_sample, data = iowa18_stan, seed = 12345, chains = 4, iter = 2000, init = 0, 
#                           pars = c("betaU", "betaC", "delta", "sigma", "lp__", "log_lik", "y_rep", "pp"))
# 
# saveRDS(ln_sample_fit, "/home/chy20005/res/res_stan_fix_0528/ln_sample2.rds")
# 


ln_sample <- stan_model("/home/chy20005/Recidivism/simulation/LN2.stan")
ln_sample_fit <- sampling(ln_sample, data = iowa18_stan, seed = 12345, chains = 4, iter = 2000, init = 0, 
                          pars = c("betaU", "betaC", "sigma", "delta", "lp__", "log_lik"))

saveRDS(ln_sample_fit, "/home/chy20005/Recidivism/simulation/res/ln2.rds")
