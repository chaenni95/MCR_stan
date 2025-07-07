library(rstan)
library(tidyverse)
library(ggplot2)
library(xtable)
library(tictoc)
library(bayesplot)
library(loo)

##### simulation ############################################################################################

n <- 1000
alpha <- 1.5
beta <- c(1,-1,-2)
eta <- c(0.5,0.6)

p <- length(beta)
q <- length(eta)

x <- matrix(rnorm(n*(p-1)),n,p-1)
x <- cbind(1,x)

w <- scale(matrix(runif(n*(q),-1,1),n,q)) # no intercept for logit!

censoring_time<-100


# Number of replications
num_reps <- 100 

models <- list(
  wb1 = stan_model("~~/Rcode/Recidivism/simulation/wb1.stan"),
  wb2 = stan_model("~~/Rcode/Recidivism/simulation/wb2.stan"),
  # wb3 = stan_model("~~/Rcode/Recidivism/simulation/wb3.stan"),
  # wb4 = stan_model("~~/Rcode/Recidivism/simulation/wb4.stan"),
  ln1 = stan_model("~~/Rcode/Recidivism/simulation/LN1.stan"),
  ln2 = stan_model("~~/Rcode/Recidivism/simulation/LN2.stan"),
  # ln3 = stan_model("~~/Rcode/Recidivism/simulation/LN3.stan"),
  # ln4 = stan_model("~~/Rcode/Recidivism/simulation/LN4.stan"),
  llg1 = stan_model("~~/Rcode/Recidivism/simulation/LLG1.stan"),
  llg2 = stan_model("~~/Rcode/Recidivism/simulation/LLG2.stan")
  # llg3 = stan_model("~~/Rcode/Recidivism/simulation/LLG3.stan")
  # llg4 = stan_model("~~/Rcode/Recidivism/simulation/LLG4.stan")
)


model_names <- names(models)
num_models <- length(models)

# Storage
selected_model_loo <- character(num_reps)
selected_model_waic <- character(num_reps)
selected_model_dic <- character(num_reps)

for (i in 1:num_reps) {
  cat("Running simulation", i, "...\n")
  
  # Generate data from the true model
  data <- geraMCM.Loglogistic(n, x, w, censoring_time = 100,
                            beta = beta, eta = eta, alpha = alpha,
                            link = "slogit", tau = 1.5, delta = 1)
  
  stan_data <- list(
    N = length(data$observed_time), 
    t = data$observed_time, 
    is_censored = data$status, 
    cured = data$cured, 
    x = data$x, 
    z = data$w, 
    M = dim(data$x)[2], 
    K = dim(data$w)[2]
  )
  
  scores <- data.frame(model = model_names, loo = NA, waic = NA, dic = NA)
  
  for (j in seq_along(models)) {
    fit <- sampling(models[[j]], data = stan_data, chains = 1, iter = 1000, refresh = 0)
    log_lik <- extract_log_lik(fit)
    
    loo_out <- loo(log_lik)
    waic_out <- waic(log_lik)
    summ <- summary(fit)$summary
    
    lp_mean <- summ["lp__", "mean"]
    lp_med <- summ["lp__", "50%"]
    deviance <- -2 * lp_mean
    pD <- 2 * (lp_mean - lp_med)
    dic <- deviance + pD
    
    scores$loo[j] <- loo_out$estimates["looic", "Estimate"]
    scores$waic[j] <- waic_out$estimates["waic", "Estimate"]
    scores$dic[j] <- dic
  }
  
  selected_model_loo[i] <- scores$model[which.min(scores$loo)]
  selected_model_waic[i] <- scores$model[which.min(scores$waic)]
  selected_model_dic[i] <- scores$model[which.min(scores$dic)]
}

# Summary of model selection frequency
a1 <- table(selected_model_loo)
a2 <- table(selected_model_waic)
a3 <- table(selected_model_dic)


saveRDS(list(a1, a2, a3), "~/Rcode/Recidivism/simulation3/res1000/llg.assess.rds")

