library(rstan)
library(ggplot2)
library(bayesplot)
library(xtable)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)

## calculate coverage of the true parameter values


compute_coverage <- function(sim_results, true_vals) {
  
  # Build one long data frame from all simulation summaries
  summary_long_df <- map_dfr(seq_along(sim_results), function(i) {
    res <- sim_results[[i]]
    
    # Filter unwanted parameters
    keep_params <- !grepl("log_lik|lp_", names(res$estimates))
    param_names <- names(res$estimates)[keep_params]
    
    tibble(
      sim_id = i,
      parameter = param_names,
      mean = res$estimates[keep_params],
      lower = res$lower_ci[keep_params],
      upper = res$upper_ci[keep_params]
    )
  })
  
  # Filter only parameters with known true values
  coverage_df <- summary_long_df %>%
    filter(parameter %in% names(true_vals)) %>%
    mutate(
      true_value = true_vals[parameter],
      covered = (lower <= true_value) & (upper >= true_value)
    )
  
  # Summarize coverage rate
  coverage_summary <- coverage_df %>%
    group_by(parameter) %>%
    summarise(
      coverage = mean(covered),
      .groups = "drop"
    )
  
  return(coverage_summary)
}

## simulation summary results 

simulation_summaries.wb1.1000 <- lapply(results.sim.wb1.1000, extract_results) # weibull logit n = 1000 result extract
simulation_summaries.wb2.1000 <- lapply(results.sim.wb2.1000, extract_results) # weibull slogit n = 1000 result extract
simulation_summaries.wb3.1000 <- lapply(results.sim.wb3.1000, extract_results) # weibull rplogit n = 1000 result extract
simulation_summaries.wb4.1000 <- lapply(results.sim.wb4.1000, extract_results) # weibull fglogit n = 1000 result extract



true_wb.l <- c(1, -1, -2, 0.5, 0.6, 1 / 1.5)
true_ln.l <- c(1, -1, -2, 0.5, 0.6, 1.5)
true_llg.l <- c(1, -1, -2, 0.5, 0.6, 1 / 1.5)

true_wb.s <- c(1, -1, -2, 0.5, 0.6, 1.5, 1 / 1.5)
true_ln.s <- c(1, -1, -2, 0.5, 0.6, 1.5, 1.5)
true_llg.s <- c(1, -1, -2, 0.5, 0.6, 1.5, 1 / 1.5)

true_wb.rp <- c(1, -1, -2, 0.5, 0.6, 1.5, 1 / 1.5)
true_ln.rp <- c(1, -1, -2, 0.5, 0.6, 1.5, 1.5)
true_llg.rp <- c(1, -1, -2, 0.5, 0.6, 1.5, 1 / 1.5)

true_wb.fg <- c(1, -1, -2, 0.5, 0.6, 0.5, 1.5, 1 / 1.5)
true_ln.fg <- c(1, -1, -2, 0.5, 0.6, 0.5, 1.5, 1.5)
true_llg.fg <- c(1, -1, -2, 0.5, 0.6, 0.5, 1.5, 1 / 1.5)


# example use ###################################################################################################################
true_wbfg <- c(1, -1, -2, 0.5, 0.6, 0.5, 1.5, 1 / 1.5)
names(true_wbfg) <- c("betaU[1]", "betaU[2]", "betaU[3]", "betaC[1]", "betaC[2]", "alpha_fg", "lambda_fg", "sigma")

coverage_wbfg <- lapply(list(simulation_summaries.wb4.500), compute_coverage, true_vals = true_wbfg)

print(coverage_wbfg)




