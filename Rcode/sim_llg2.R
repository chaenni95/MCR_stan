library(rstan)
library(tidyverse)
library(ggplot2)
library(xtable)
library(tictoc)
library(bayesplot)
library(loo)

# Generalized loglogistic distribution

gen_loglogistic <- function(n, alpha, lambda) {
  # Generate n random numbers from the Generalized loglogistic distribution
  U <- runif(n)
  T <- lambda * (U / (1 - U))^(1 / alpha)
  return(T)
}

loglogistic_pdf <- function(x, alpha, lambda) {
  pdf <- (alpha / lambda) * (x / lambda)^(alpha - 1) / (1 + (x / lambda)^alpha)^2
  return(pdf)
}

loglogistic_sf <- function(x, alpha, lambda) {
  sf <- 1 / (1 + (x / lambda)^alpha)
  return(sf)
}


cen_dat <- function(perc, T, n, cen_type="left"){
  if(cen_type=="left"){
    aa = sort(T, decreasing=FALSE)
    cutoff<-aa[ceiling(perc*n)]
    
    cc=matrix(1,n,1)
    
    for(i in 1:n){
      cc[i]<-cc[i]*(T[i]>cutoff[1])
      if(cc[i]==0){T[i]=cutoff[1]}
    }
  }else{
    aa = sort(T, decreasing=TRUE)
    cutoff<-aa[ceiling(perc*n)]
    
    cc=matrix(1,n,1)
    
    for(i in 1:n){
      cc[i]<-cc[i]*(T[i]<cutoff[1])
      if(cc[i]==0){T[i]=cutoff[1]}
    }
  }
  return(list(cutoff=cutoff,T=T, cc=cc))
}



logLik.Loglogistic <- function(cc, x, alpha, lambda){
  logLik <- sum(cc*log(loglogistic_pdf(x, alpha, lambda)))+sum((1-cc)*log(loglogistic_sf(x, alpha, lambda)))
  return(logLik)
}


logLik.MCM.Loglogistic <- function(cc, uncured, x, alpha, lambda, theta){
  f1<-f2<-rep(0, length(cc))
  f1 <- ifelse(log(theta*loglogistic_sf(x, alpha, lambda))==-Inf, 0, log(theta*loglogistic_sf(x, alpha, lambda)))
  f2 <- ifelse(log(theta*cc*uncured*loglogistic_pdf(x, alpha, lambda))==-Inf, 0, log(theta*cc*uncured*loglogistic_pdf(x, alpha, lambda)))
  
  logLik <- sum(f1)+
    sum((1-cc)*uncured*(f2))+
    sum((1-cc)*(1-uncured)*log(1-theta))
  return(logLik)
}


################################################################################
## creating links
################################################################################

plogit <- function(tau = 1) {
  
  # Link function: maps from probabilities (mu) to the linear predictor (eta)
  linkfun <- function(mu) {
    (log(mu^(1/tau) / (1 - mu^(1/tau)))) # Standard logit link
  }
  
  # Inverse link function: maps from the linear predictor (eta) to probabilities (mu)
  linkinv <- function(eta) {
    (1 / (1 + exp(-eta)))^tau  # Logistic function raised to power tau
  }
  
  # Derivative of the inverse link function with respect to eta
  mu.eta <- function(eta) {
    logistic <- 1 / (1 + exp(-eta))  # Standard logistic function
    tau * logistic^(tau - 1) * exp(-eta) / (1 + exp(-eta))^2  # Derivative with respect to eta
  }
  
  # Valid eta check (always valid for this transformation)
  valideta <- function(eta) TRUE
  
  # Create and return the custom link list that glm() or vglm() will use
  link <- list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, valideta = valideta, name = paste("plogit(tau=", tau, ")", sep=""))
  class(link) <- "link-glm"
  return(link)
}


####example
#data <- data.frame(x = rnorm(100), y = rbinom(100, 1, prob = 0.5))
#model <- glm(y ~ x, family = binomial(link =plogit(tau = 0.5)), data = data)
#summary(model)

### generating Mixture Cure Model with loglogistic base

geraMCM.Loglogistic <- function(n, x, w, censoring_time=10 , beta=c(1,2), eta=c(1,2), alpha=1.5, link="plogit", tau=1.5, delta = 1){
  #### tau=skewness parameter of power links
  if (link == "logit") {p_cure <- 1 / (1 + exp(-w %*% eta))}
  if (link == "probit") {p_cure <- pnorm(w %*% eta)}
  if (link == "plogit") {p_cure <- (1 / (1 + exp(-w %*% eta)))^tau}
  if (link == "rplogit") {p_cure <- (1 - 1/(1 + exp(w %*% eta))^tau)}
  if (link == "fglogit") {p_cure <- (1 - (1/(1+exp(w %*% eta)))^tau)^(delta)}
  
  # Generate cure status
  cured <- rbinom(n, 1, p_cure)  # Cured = 1, Uncured = 0
  lambda <- exp(x %*% beta)  # Scale parameter
  
  event_time <- rep(NA, n)
  non_cured <- which(cured == 0)  # Identify non-cured individuals
  u <- runif(length(non_cured))  # Uniform(0,1) random variables
  
  # Log-Logistic inverse transform sampling:
  event_time[non_cured] <- ( (u / (1 - u))^(1/alpha) ) * lambda[non_cured]
  censor_time <- runif(n, 0, censoring_time)
  
  observed_time <- pmin(event_time, censor_time, na.rm = TRUE)
  status <- ifelse(cured == 1, 0, ifelse(event_time <= censor_time, 1, 0))  # 0 = censored, 1 = event
  return(list(x = x, w= w, cured = cured, lambda = lambda, event_time = event_time, censor_time = censor_time,
              observed_time = observed_time, status = status))  
}



##### simulation ############################################################################################


n <- 500
alpha <- 1.5
beta <- c(1,-1,-2)
eta <- c(0.5,0.6)

p <- length(beta)
q <- length(eta)

x <- matrix(rnorm(n*(p-1)),n,p-1)
x <- cbind(1,x)

w <- scale(matrix(runif(n*(q),-1,1),n,q)) # no intercept for logit!

censoring_time<-10


# Number of replications
num_reps <- 1000  

# Storage for results
results <- vector("list", num_reps)

# Compile the Stan model once to avoid re-compiling
stanmodel <- stan_model("/home/chy20005/Recidivism/simulation/LLG2.stan")

for (i in 1:num_reps) {
  
  # Step 1: Generate a new dataset (modify this function to fit your simulation setup)
  data <- geraMCM.Loglogistic(n, x, w, censoring_time=censoring_time , beta=beta, eta=eta, alpha=alpha,link="plogit", tau=1.5, delta = 1)  # Custom function to simulate data
  
  # Step 2: Prepare Stan data
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
  
  # Step 3: Fit the model
  fit <- sampling(stanmodel, data = stan_data, seed = 12345, chains = 1, iter = 1000, 
                  init = 0, pars = c("betaU", "betaC", "sigma", "delta", "lp__", "log_lik"))
  
  # Step 4: Extract posterior estimates (Modify `extract_results()` as needed)
  results[[i]] <- fit  
}


saveRDS(results, "/home/chy20005/Recidivism/simulation/res2/simllg2.rds")






# Step 1: Extract posterior summaries from each fitted model
extract_results <- function(fit) {
  
  # Extract posterior summaries
  summary_fit <- summary(fit)$summary
  
  # Compute 95% credible interval (CI)
  estimates <- summary_fit[, "mean"]   # Posterior mean
  median_estimates <- summary_fit[, "50%"]  # Posterior median
  lower_ci <- summary_fit[, "2.5%"]   # 2.5% quantile
  upper_ci <- summary_fit[, "97.5%"]  # 97.5% quantile
  
  # Compute LOO-CV
  log_lik <- extract_log_lik(fit, merge_chains = FALSE)
  loo_results <- loo(log_lik)
  
  # Compute WAIC
  waic_results <- waic(log_lik)
  
  # Compute DIC (Deviance Information Criterion)
  deviance <- -2 * summary_fit["lp__", "mean"]
  pD <- 2 * (summary_fit["lp__", "mean"] - summary_fit["lp__", "50%"])  # Effective parameters
  dic <- deviance + pD  # DIC Formula
  
  # Return a list with all metrics
  return(list(
    estimates = estimates,
    median_estimates = median_estimates,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    loo = loo_results,
    waic = waic_results,
    dic = dic
  ))
}

# Step 2: Apply the function to all simulation results
simulation_summaries <- lapply(results, extract_results)

# Step 3: Aggregate the results across all simulations
summarize_simulations <- function(simulation_summaries) {
  
  # Convert list to a data frame
  estimates_list <- do.call(rbind, lapply(simulation_summaries, function(res) res$estimates))
  lower_ci_list <- do.call(rbind, lapply(simulation_summaries, function(res) res$lower_ci))
  upper_ci_list <- do.call(rbind, lapply(simulation_summaries, function(res) res$upper_ci))
  
  # Compute mean & standard deviation across simulations
  summary_df <- data.frame(
    Parameter = colnames(estimates_list),
    Mean_Estimate = colMeans(estimates_list, na.rm = TRUE),
    SD_Estimate = apply(estimates_list, 2, sd, na.rm = TRUE),
    Lower_95_CI = colMeans(lower_ci_list, na.rm = TRUE),
    Upper_95_CI = colMeans(upper_ci_list, na.rm = TRUE)
  )
  
  # Compute mean LOO, WAIC, and DIC across simulations
  mean_loo <- mean(sapply(simulation_summaries, function(res) res$loo)["looic", ] %>% unlist(), na.rm = T)
  mean_waic <- mean(sapply(simulation_summaries, function(res) res$waic)["waic", ] %>% unlist(),  na.rm = T)
  mean_dic <- mean(sapply(simulation_summaries, function(res) res$dic), na.rm = T)
  
  # Print summary statistics
  print(summary_df[!grepl("^log_lik", summary_df$Parameter), ])
  
  # Print model comparison metrics
  cat("\nModel Evaluation Metrics (Averaged Over Simulations):\n")
  cat("Mean LOO: ", mean_loo, "\n")
  cat("Mean WAIC: ", mean_waic, "\n")
  cat("Mean DIC: ", mean_dic, "\n")
}

# Run the summary function
summarize_simulations(simulation_summaries)


saveRDS(summarize_simulations(simulation_summaries), "/home/chy20005/Recidivism/simulation/res2/simllg2_summary.rds")




