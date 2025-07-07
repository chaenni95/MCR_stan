library(rstan)
library(ggplot2)
library(bayesplot)
library(xtable)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)


plot_coxsnell_envelope <- function(dist, link, dir_path, pdf_path = NULL) {
  # Generate tag and filename
  tag <- paste(dist, link, sep = "_")
  filename <- paste0("coxsnell_", tag, ".RData")
  filepath <- file.path(dir_path, filename)
  
  # Load data
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }
  load(filepath)
  
  # Extract tagged variables
  grid_times <- get(paste0("grid_times_", tag))
  lower <- get(paste0("lower_", tag))
  upper <- get(paste0("upper_", tag))
  km_df <- get(paste0("km_df_", tag))
  exp_curve <- get(paste0("exp_curve_", tag))
  
  # Save to PDF if pdf_path is given
  if (!is.null(pdf_path)) {
    pdf(file = pdf_path, width = 7, height = 5)
    on.exit(dev.off())  # ensure device closes even if error occurs
  }
  
  # Plot
  plot(km_df$time, km_df$surv, type = "s",
       xlab = "Cox-Snell Residual", ylab = "Survival",
       main = paste("Pooled KM Cox-Snell Residuals with 95% Envelope\n(", tag, ")"),
       xlim = range(grid_times), ylim = c(0, 1),
       col = "black", lwd = 2)
  
  curve(exp(-x), add = TRUE, col = "red", lty = 2, lwd = 2)
  lines(grid_times, lower, col = "gray", lty = 3)
  lines(grid_times, upper, col = "gray", lty = 3)
  
  legend("topright",
         legend = c("KM of Residuals", "Exponential(1)", "95% Envelope"),
         col = c("black", "red", "gray"),
         lty = c(1, 2, 3), lwd = c(2, 2, 1), bty = "n")
}



### example with simulation of n = 500 ############################################################################

plot_coxsnell_envelope(
  dist = "weibull",
  link = "logit",
  dir_path = "~/Recidivism/simulation3/res1000",
  pdf_path = "~/Recidivism/simulation3/res1000/coxsnell_weibull_logit.500.pdf"
)


plot_coxsnell_envelope(
  dist = "weibull",
  link = "slogit",
  dir_path = "~/Recidivism/simulation3/res1000",
  pdf_path = "~/Recidivism/simulation3/res1000/coxsnell_weibull_slogit.500.pdf"
)

plot_coxsnell_envelope(
  dist = "weibull",
  link = "rplogit",
  dir_path = "~/Recidivism/simulation3/res1000",
  pdf_path = "~/Recidivism/simulation3/res1000/coxsnell_weibull_rplogit.500.pdf"
)

plot_coxsnell_envelope(
  dist = "weibull",
  link = "fglogit",
  dir_path = "~/Recidivism/simulation3/res1000",
  pdf_path = "~/Recidivism/simulation3/res1000/coxsnell_weibull_fglogit.500.pdf"
)





plot_coxsnell_envelope(
  dist = "lognormal",
  link = "logit",
  dir_path = "~/Recidivism/simulation3/res1000",
  pdf_path = "~/Recidivism/simulation3/res1000/coxsnell_lognormal_logit.500.pdf"
)


plot_coxsnell_envelope(
  dist = "lognormal",
  link = "slogit",
  dir_path = "~/Recidivism/simulation3/res1000",
  pdf_path = "~/Recidivism/simulation3/res1000/coxsnell_lognormal_slogit.500.pdf"
)

plot_coxsnell_envelope(
  dist = "lognormal",
  link = "rplogit",
  dir_path = "~/Recidivism/simulation3/res1000",
  pdf_path = "~/Recidivism/simulation3/res1000/coxsnell_lognormal_rplogit.500.pdf"
)

plot_coxsnell_envelope(
  dist = "lognormal",
  link = "fglogit",
  dir_path = "~/Recidivism/simulation3/res1000",
  pdf_path = "~/Recidivism/simulation3/res1000/coxsnell_lognormal_fglogit.500.pdf"
)






plot_coxsnell_envelope(
  dist = "loglogistic",
  link = "logit",
  dir_path = "~/Recidivism/simulation3/res1000",
  pdf_path = "~/Recidivism/simulation3/res1000/coxsnell_loglogistic_logit.500.pdf"
)


plot_coxsnell_envelope(
  dist = "loglogistic",
  link = "slogit",
  dir_path = "~/Recidivism/simulation3/res1000",
  pdf_path = "~/Recidivism/simulation3/res1000/coxsnell_loglogistic_slogit.500.pdf"
)

plot_coxsnell_envelope(
  dist = "loglogistic",
  link = "rplogit",
  dir_path = "~/Recidivism/simulation3/res1000",
  pdf_path = "~/Recidivism/simulation3/res1000/coxsnell_loglogistic_rplogit.500.pdf"
)

plot_coxsnell_envelope(
  dist = "loglogistic",
  link = "fglogit",
  dir_path = "~/Recidivism/simulation3/res1000",
  pdf_path = "~/Recidivism/simulation3/res1000/coxsnell_loglogistic_fglogit.500.pdf"
)
