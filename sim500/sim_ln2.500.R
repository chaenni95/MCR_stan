library(rstan)
library(tidyverse)
library(ggplot2)
library(xtable)
library(tictoc)
library(bayesplot)
library(loo)


##### simulation ############################################################################################


n <- 500
alpha <- 1.5
beta <- c(1,-1,-2)
eta <- c(0.5,0.6)

p <- length(beta)
q <- length(eta)

x <- matrix(rnorm(n*(p-1)),n,p-1)
x <- cbind(1,x)

w <- scale(matrix(runif(n*(q),-1,1),n,q)) # no intercept for slogit!

censoring_time<-100


# Number of replications
num_reps <- 1000  

# Storage for results
results <- vector("list", num_reps)

# Compile the Stan model once to avoid re-compiling
stanmodel <- stan_model("~/Recidivism/simulation/LN2.stan")

# stanmodel <- stan_model("/Users/chaenni/Documents/UCONN/Research/Paper1/simulation/LN2.stan")

for (i in 1:num_reps) {
  
  # Step 1: Generate a new dataset (modify this function to fit your simulation setup)
  data <- geraMCM.Lognormal(n, x, w, censoring_time=censoring_time , beta=beta, eta=eta, alpha=alpha,link="slogit", tau=1.5, delta = 1)  # Custom function to simulate data
  
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
                  init = 0, pars = c("betaU", "betaC", "sigma", "alpha_s" "lp__", "log_lik"))
  
  # Step 4: Extract posterior estimates (Modify `extract_results()` as needed)
  results[[i]] <- list(
    fit = fit,
    data = list(
      X = data$x,
      Z = data$w,
      time = data$observed_time,
      status = data$status,
      event_time = data$event_time  # For quantile check
    )
  )
}



saveRDS(results, "~/Recidivism/res/simln2.rds")



# Step 2: Apply the function to all simulation results
simulation_summaries <- lapply(results, extract_results)
saveRDS(simulation_summaries, "~/Recidivism/res500/simln2_summary.rds")

# Step 3: Aggregate the results across all simulations
summarize_simulations(simulation_summaries)



# Step 4: Extract result for cox-snell residual plot 
cs_all <- lapply(results, extract_cs_residual, dist = "lognormal", link = "slogit")
all_r <- unlist(lapply(cs_all, `[[`, "residuals"))
all_status <- unlist(lapply(cs_all, `[[`, "status"))

km_all <- survfit(Surv(all_r, all_status) ~ 1)
km_all_df <- data.frame(time = km_all$time, surv = km_all$surv)

grid_times <- seq(0, max(all_r), length.out = 100)

# Compute theoretical exponential curve on same grid
exp_curve <- data.frame(time = grid_times, surv = exp(-grid_times))

# Matrix to store interpolated KM survival at each grid time for each simulation
km_surv_matrix <- matrix(NA, nrow = length(grid_times), ncol = length(cs_all))

for (i in seq_along(cs_all)) {
  km_i <- survfit(Surv(cs_all[[i]]$residuals, cs_all[[i]]$status) ~ 1)
  
  # Interpolate KM survival onto the fixed grid
  interp <- approx(km_i$time, km_i$surv, xout = grid_times, method = "constant", rule = 2)$y
  km_surv_matrix[, i] <- interp
}

lower <- apply(km_surv_matrix, 1, quantile, probs = 0.025, na.rm = TRUE)
upper <- apply(km_surv_matrix, 1, quantile, probs = 0.975, na.rm = TRUE)

# Enforce monotonicity for aesthetic survival envelope
lower <- cummin(lower)
upper <- cummin(upper)


dist <- "lognormal"
link <- "slogit"
tag <- paste(dist, link, sep = "_")  # e.g., "lognormal_slogit"

assign(paste0("grid_times_", tag), grid_times)
assign(paste0("lower_", tag), lower)
assign(paste0("upper_", tag), upper)
assign(paste0("km_df_", tag), km_all_df)
assign(paste0("exp_curve_", tag), exp_curve)

outdir <- "~/Recidivism/res500"
filename <- paste0("coxsnell_", tag, ".RData")
filepath <- file.path(outdir, filename)

# Save only these variables
save(list = c(paste0("grid_times_", tag),
              paste0("lower_", tag),
              paste0("upper_", tag),
              paste0("km_df_", tag),
              paste0("exp_curve_", tag)),
     file = filepath)

