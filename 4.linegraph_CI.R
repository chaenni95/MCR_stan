library(rstan)
library(tidyverse)
library(ggplot2)
library(xtable)
library(tictoc)
library(bayesplot)
library(loo)

## data frame of true values 

true_wb.srp <- c(1, -1, -2, 0.5, 0.6, 1.5, 1 / 1.5)
true_ln.srp <- c(1, -1, -2, 0.5, 0.6, 1.5, 1.5)
true_llg.srp <- c(1, -1, -2, 0.5, 0.6, 1.5, 1 / 1.5)


true_wb.l <- c(1, -1, -2, 0.5, 0.6, 1 / 1.5)
true_ln.l <- c(1, -1, -2, 0.5, 0.6, 1.5)
true_llg.l <- c(1, -1, -2, 0.5, 0.6, 1 / 1.5)


true_wb.fg <- c(1, -1, -2, 0.5, 0.6, 0.5, 1.5, 1 / 1.5)
true_ln.fg <- c(1, -1, -2, 0.5, 0.6, 0.5, 1.5, 1.5)
true_llg.fg <- c(1, -1, -2, 0.5, 0.6, 0.5, 1.5, 1 / 1.5)



parameter.l <- c("betaU[1]", "betaU[2]", "betaU[3]", "betaC[1]", "betaC[2]", "sigma")
true.wb.l <- data.frame(TrueValue = true_wb.l, Parameter = parameter.l)
true.ln.l <- data.frame(TrueValue = true_ln.l, Parameter = parameter.l)
true.llg.l <- data.frame(TrueValue = true_llg.l, Parameter = parameter.l)

parameter.s <- c("betaU[1]", "betaU[2]", "betaU[3]", "betaC[1]", "betaC[2]", "alpha_s", "sigma")
true.wb.s <- data.frame(TrueValue = true_wb.s, Parameter = parameter.s)
true.ln.s <- data.frame(TrueValue = true_ln.s, Parameter = parameter.s)
true.llg.s <- data.frame(TrueValue = true_llg.s, Parameter = parameter.s)


parameter.rp <- c("betaU[1]", "betaU[2]", "betaU[3]", "betaC[1]", "betaC[2]", "lambda_rp", "sigma")
true.wb.rp <- data.frame(TrueValue = true_wb.rp, Parameter = parameter.rp)
true.ln.rp <- data.frame(TrueValue = true_ln.rp, Parameter = parameter.rp)
true.llg.rp <- data.frame(TrueValue = true_llg.rp, Parameter = parameter.rp)


parameter.fg <- c("betaU[1]", "betaU[2]", "betaU[3]", "betaC[1]", "betaC[2]", "alpha_fg", "lambda_fg", "sigma")
true.wb.fg <- data.frame(TrueValue = true_wb.fg, Parameter = parameter.fg)
true.ln.fg <- data.frame(TrueValue = true_ln.fg, Parameter = parameter.fg)
true.llg.fg <- data.frame(TrueValue = true_llg.fg, Parameter = parameter.fg)




lines_plot_simulation_estimates_by_sample_size <- function(sim_list, sim_labels, true.values, exclude_pattern = c("^log_lik", "lp__")) {
  if (length(sim_list) != length(sim_labels)) {
    stop("`sim_list` and `sim_labels` must be the same length.")
  }
  
  # Step 1: Gather all simulation summary data
  summary_data <- do.call(rbind, lapply(seq_along(sim_list), function(i) {
    sim <- sim_list[[i]]
    label <- sim_labels[i]
    
    do.call(rbind, lapply(sim, function(res) {
      data.frame(
        Parameter = names(res$estimates),
        Mean = res$estimates,
        Lower = res$lower_ci,
        Upper = res$upper_ci,
        SampleSize = label,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  # Step 2: Exclude parameters by pattern
  if (!is.null(exclude_pattern)) {
    pattern_combined <- paste(exclude_pattern, collapse = "|")
    summary_data <- subset(summary_data, !grepl(pattern_combined, Parameter))
  }
  
  # Step 3: Aggregate by parameter and sample size
  library(dplyr)
  summary_df <- summary_data %>%
    group_by(Parameter, SampleSize) %>%
    summarise(
      Mean = mean(Mean, na.rm = TRUE),
      Lower = mean(Lower, na.rm = TRUE),
      Upper = mean(Upper, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  # Step 5: Plot
  library(ggplot2)
  p <- ggplot(summary_df, aes(x = Parameter, y = Mean, color = SampleSize)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
                  width = 0.2,
                  position = position_dodge(width = 0.5)) +
    geom_point(data = true.values,
               aes(x = Parameter, y = TrueValue),
               shape = 18, size = 2, color = "red", inherit.aes = FALSE) +
    geom_segment(data = true.values,
                 aes(x = as.numeric(factor(Parameter)) - 0.3,
                     xend = as.numeric(factor(Parameter)) + 0.3,
                     y = TrueValue, yend = TrueValue),
                 color = "red", linewidth = 0.7, inherit.aes = FALSE) +
    scale_color_grey(start = 0.3, end = 0.8) +
    labs(title = "",
         x = "Parameter", y = "Estimate", color = "Sample Size") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  return(p)
}




# extracted result from simulation ###############################################################################################
simulation_summaries.wb1.1000 <- lapply(results.sim.wb1.1000, extract_results) # weibull logit n = 1000 result extract
simulation_summaries.wb2.1000 <- lapply(results.sim.wb2.1000, extract_results) # weibull slogit n = 1000 result extract
simulation_summaries.wb3.1000 <- lapply(results.sim.wb3.1000, extract_results) # weibull rplogit n = 1000 result extract
simulation_summaries.wb4.1000 <- lapply(results.sim.wb4.1000, extract_results) # weibull fglogit n = 1000 result extract


simulation_summaries.ln1.1000 <- lapply(results.sim.ln1.1000, extract_results) # lognormal logit n = 1000 result extract
simulation_summaries.ln2.1000 <- lapply(results.sim.ln2.1000, extract_results) # lognormal slogit n = 1000 result extract
simulation_summaries.ln3.1000 <- lapply(results.sim.ln3.1000, extract_results) # lognormal rplogit n = 1000 result extract
simulation_summaries.ln4.1000 <- lapply(results.sim.ln4.1000, extract_results) # lognormal fglogit n = 1000 result extract


simulation_summaries.llg1.1000 <- lapply(results.sim.llg1.1000, extract_results) # loglogistic logit n = 1000 result extract
simulation_summaries.llg2.1000 <- lapply(results.sim.llg2.1000, extract_results) # loglogistic slogit n = 1000 result extract
simulation_summaries.llg3.1000 <- lapply(results.sim.llg3.1000, extract_results) # loglogistic rplogit n = 1000 result extract
simulation_summaries.llg4.1000 <- lapply(results.sim.llg4.1000, extract_results) # loglogistic fglogit n = 1000 result extract



simulation_summaries.wb1.500 <- lapply(results.sim.wb1.500, extract_results) # weibull logit n = 500 result extract
simulation_summaries.wb2.500 <- lapply(results.sim.wb2.500, extract_results) # weibull slogit n = 500 result extract
simulation_summaries.wb3.500 <- lapply(results.sim.wb3.500, extract_results) # weibull rplogit n = 500 result extract
simulation_summaries.wb4.500 <- lapply(results.sim.wb4.500, extract_results) # weibull fglogit n = 500 result extract


simulation_summaries.ln1.500 <- lapply(results.sim.ln1.500, extract_results) # lognormal logit n = 500 result extract
simulation_summaries.ln2.500 <- lapply(results.sim.ln2.500, extract_results) # lognormal slogit n = 500 result extract
simulation_summaries.ln3.500 <- lapply(results.sim.ln3.500, extract_results) # lognormal rplogit n = 500 result extract
simulation_summaries.ln4.500 <- lapply(results.sim.ln4.500, extract_results) # lognormal fglogit n = 500 result extract


simulation_summaries.llg1.500 <- lapply(results.sim.llg1.500, extract_results) # loglogistic logit n = 500 result extract
simulation_summaries.llg2.500 <- lapply(results.sim.llg2.500, extract_results) # loglogistic slogit n = 500 result extract
simulation_summaries.llg3.500 <- lapply(results.sim.llg3.500, extract_results) # loglogistic rplogit n = 500 result extract
simulation_summaries.llg4.500 <- lapply(results.sim.llg4.500, extract_results) # loglogistic fglogit n = 500 result extract



### Weibull
wb1.plot <- lines_plot_simulation_estimates_by_sample_size(
  sim_list = list(simulation_summaries.wb1.500, simulation_summaries.wb1.1000),
  sim_labels = c("n = 500", "n = 1000"), true.values = true.wb.l) + 
  labs(title = "logit link", x = "Parameter", y = "Estimate", color = "Sample Size")


wb2.plot <- lines_plot_simulation_estimates_by_sample_size(
  sim_list = list(simulation_summaries.wb2.500, simulation_summaries.wb2.1000),
  sim_labels = c("n = 500", "n = 1000"), true.values = true.wb.srp)+ 
  labs(title = "slogit link", x = "Parameter", y = "Estimate", color = "Sample Size")


wb3.plot <- lines_plot_simulation_estimates_by_sample_size(
  sim_list = list(simulation_summaries.wb3.500, simulation_summaries.wb3.1000),
  sim_labels = c("n = 500", "n = 1000"), true.values = true.wb.srp)+ 
  labs(title = "rplogit link", x = "Parameter", y = "Estimate", color = "Sample Size")



wb4.plot <- lines_plot_simulation_estimates_by_sample_size(
  sim_list = list(simulation_summaries.wb4.500, simulation_summaries.wb4.500),
  sim_labels = c("n = 500", "n = 1000"), true.values = true.wb.fg)+ 
  labs(title = "fglogit link", x = "Parameter", y = "Estimate", color = "Sample Size")




### Lognormal
ln1.plot <- lines_plot_simulation_estimates_by_sample_size(
  sim_list = list(simulation_summaries.ln1.500, simulation_summaries.ln1.1000),
  sim_labels = c("n = 500", "n = 1000"), true.values = true.ln.l)+ 
  labs(title = "logit link", x = "Parameter", y = "Estimate", color = "Sample Size")


ln2.plot <- lines_plot_simulation_estimates_by_sample_size(
  sim_list = list(simulation_summaries.ln2.500, simulation_summaries.ln2.1000),
  sim_labels = c("n = 500", "n = 1000"), true.values = true.ln.srp)+ 
  labs(title = "slogit link", x = "Parameter", y = "Estimate", color = "Sample Size")


ln3.plot <- lines_plot_simulation_estimates_by_sample_size(
  sim_list = list(simulation_summaries.ln3.500, simulation_summaries.ln3.1000),
  sim_labels = c("n = 500", "n = 1000"), true.values = true.ln.srp)+ 
  labs(title = "rplogit link", x = "Parameter", y = "Estimate", color = "Sample Size")


ln4.plot <- lines_plot_simulation_estimates_by_sample_size(
  sim_list = list(simulation_summaries.ln4.500, simulation_summaries.ln4.1000),
  sim_labels = c("n = 500", "n = 1000"), true.values = true.ln.fg)+ 
  labs(title = "fglogit link", x = "Parameter", y = "Estimate", color = "Sample Size")




### Loglogistic
llg1.plot <- lines_plot_simulation_estimates_by_sample_size(
  sim_list = list(simulation_summaries.llg1.500, simulation_summaries.llg1.1000),
  sim_labels = c("n = 500", "n = 1000"), true.values = true.llg.l)+ 
  labs(title = "logit link", x = "Parameter", y = "Estimate", color = "Sample Size")


llg2.plot <- lines_plot_simulation_estimates_by_sample_size(
  sim_list = list(simulation_summaries.llg2.500, simulation_summaries.llg2.1000),
  sim_labels = c("n = 500", "n = 1000"), true.values = true.llg.srp)+ 
  labs(title = "slogit link", x = "Parameter", y = "Estimate", color = "Sample Size")


llg3.plot <- lines_plot_simulation_estimates_by_sample_size(
  sim_list = list(simulation_summaries.llg3.500, simulation_summaries.llg3.1000),
  sim_labels = c("n = 500", "n = 1000"), true.values = true.llg.srp)+ 
  labs(title = "rplogit link", x = "Parameter", y = "Estimate", color = "Sample Size")


llg4.plot <- lines_plot_simulation_estimates_by_sample_size(
  sim_list = list(simulation_summaries.llg4.500, simulation_summaries.llg4.1000),
  sim_labels = c("n = 500", "n = 1000"), true.values = true.llg.fg)+ 
  labs(title = "fglogit link", x = "Parameter", y = "Estimate", color = "Sample Size")



