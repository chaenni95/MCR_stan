library(rstan)
library(tidyverse)
library(ggplot2)
library(xtable)
library(tictoc)
library(bayesplot)
library(loo)

# Step 1: Extract posterior summaries from each fitted model
extract_results <- function(fit) {
  
  # Extract posterior summaries
  fit <- fit$fit
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




extract_cs_residual <- function(result, dist = "weibull", link = "logit") {
  fit_post <- rstan::extract(result$fit)
  data <- result$data
  
  n <- length(data$time)
  L <- length(fit_post$sigma)
  S_mix_mat <- matrix(NA, n, L)
  
  for (l in 1:L) {
    betaU_l <- fit_post$betaU[l, ]
    betaC_l <- fit_post$betaC[l, ]
    sigma_l <- fit_post$sigma[l]
    
    delta_l <- if (!is.null(fit_post$delta)) fit_post$delta[l] else 1
    lambda_link_l <- if (!is.null(fit_post$lambda)) fit_post$lambda[l] else 1
    
    muU_l <- as.numeric(data$X %*% betaU_l)
    zbeta <- data$Z %*% betaC_l
    
    p_cure <- switch(link,
                     "logit" = 1 / (1 + exp(-zbeta)),
                     "slogit" = (1 / (1 + exp(-zbeta)))^delta_l,
                     "rplogit" = 1 - (1 / (1 + exp(zbeta)))^delta_l,
                     "fglogit" = (1 - (1 / (1 + exp(zbeta)))^lambda_link_l)^delta_l,
                     stop("Unsupported link")
    )
    
    Su_l <- switch(dist,
                   "weibull" = {
                     alpha_l <- 1 / sigma_l
                     lambda_l <- exp(muU_l)
                     exp(-(data$time / lambda_l)^alpha_l)
                   },
                   "lognormal" = {
                     lambda_l <- muU_l
                     1 - plnorm(data$time, meanlog = lambda_l, sdlog = sigma_l)
                   },
                   "loglogistic" = {
                     shape <- 1 / sigma_l
                     scale <- exp(muU_l)
                     1 / (1 + (data$time / scale)^shape)
                   },
                   stop("Unsupported distribution")
    )
    
    S_mix_mat[, l] <- p_cure + (1 - p_cure) * Su_l
  }
  
  r_cs <- -log(rowMeans(S_mix_mat))
  list(residuals = r_cs, status = data$status)
}






