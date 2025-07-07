library(survival)
library(bayesplot)
library(posterior)
library(survminer)

iowa_rec18_dum <- read.csv("~/Recidivism/data_files/iowa18.csv")

wbl.fit <- readRDS("~/Recidivism/dataresult/wb1.rds")
wbsl.fit <- readRDS("~/Recidivism/dataresult/wb2.rds")
wbrpl.fit <- readRDS("~/Recidivism/dataresult/wb3.rds")
wbfg.fit <- readRDS("~/Recidivism/dataresult/wb4.rds")


lnl.fit <- readRDS("~/Recidivism/dataresult/ln1.rds")
lnsl.fit <- readRDS("~/Recidivism/dataresult/ln2.rds")
lnrpl.fit <- readRDS("~/Recidivism/dataresult/ln3.rds")
lnfg.fit <- readRDS("~/Recidivism/dataresult/ln4.rds")


llgl.fit <- readRDS("~/Recidivism/dataresult/llg1.rds")
llgsl.fit <- readRDS("~/Recidivism/dataresult/llg2.rds")
llgrpl.fit <- readRDS("~/Recidivism/dataresult/llg3.rds")
llgfg.fit <- readRDS("~/Recidivism/dataresult/llg4.rds")



iowa18_stan <- list(N = nrow(iowa_rec18_dum), t = iowa_rec18_dum$time_to_recurrence, 
                    is_censored = iowa_rec18_dum$is_censored, 
                    cured = iowa_rec18_dum$cured, x = cbind(1, iowa_rec18_dum[, -c(3,4,5, 32:38, 41, 42)]), 
                    z = iowa_rec18_dum[, -c(3,4,5, 32:42)], K = 28, M = 31)

X <- as.matrix(cbind(1, iowa_rec18_dum[, -c(3,4,5, 32:38, 41, 42)]))
Z <- as.matrix(iowa_rec18_dum[, -c(3,4,5, 32:42)])

time <- iowa_rec18_dum$time_to_recurrence
status <- iowa_rec18_dum$is_censored




post.fit.wbl <- rstan::extract(wbl.fit)
post.fit.wbsl <- rstan::extract(wbsl.fit)
post.fit.wbrpl <- rstan::extract(wbrpl.fit)
post.fit.wbfg <- rstan::extract(wbfg.fit)


post.fit.lnl <- rstan::extract(lnl.fit)
post.fit.lnsl <- rstan::extract(lnsl.fit)
post.fit.lnrpl <- rstan::extract(lnrpl.fit)
post.fit.lnfg <- rstan::extract(lnfg.fit)


post.fit.llgl <- rstan::extract(llgl.fit)
post.fit.llgsl <- rstan::extract(llgsl.fit)
post.fit.llgrpl <- rstan::extract(llgrpl.fit)
post.fit.llgfg <- rstan::extract(llgfg.fit)

extract_cs_residual <- function(result, dist = "weibull", link = "rplogit") {
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

cs_all <- lapply(results, extract_cs_residual, dist = "lognormal", link = "rplogit")
all_r <- unlist(lapply(cs_all, `[[`, "residuals"))
all_status <- unlist(lapply(cs_all, `[[`, "status"))

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


plot(km_all, xlab = "Cox-Snell Residual", ylab = "Survival",
     main = "Pooled KM Cox-Snell Residuals with 95% Envelope",
     xlim = range(grid_times), ylim = c(0, 1))

curve(exp(-x), add = TRUE, col = "red", lty = 2)
lines(grid_times, lower, col = "gray", lty = 3)
lines(grid_times, upper, col = "gray", lty = 3)

legend("topright",
       legend = c("KM of Residuals", "Exponential(1)", "95% Envelope"),
       col = c("black", "red", "gray"), lty = c(1, 2, 3), bty = "n")



plot_cs_envelope <- function(results,
                             data,
                             time,
                             status,
                             X,
                             Z,
                             dist = "weibull",
                             link = "logit",
                             tag = NULL,
                             save_path = NULL) {
  
  if (is.null(tag)) tag <- paste(dist, link, sep = "_")
  
  # --- Helper: Extract Cox-Snell residuals from one Stan fit result ---
  extract_cs_residual <- function(result) {
    fit_post <- rstan::extract(result)
    
    n <- length(time)
    L <- length(fit_post$sigma)
    S_mix_mat <- matrix(NA, n, L)
    
    for (l in 1:L) {
      betaU_l <- fit_post$betaU[l, ]
      betaC_l <- fit_post$betaC[l, ]
      sigma_l <- fit_post$sigma[l]
      delta_l <- if (!is.null(fit_post$delta)) fit_post$delta[l] else 1
      lambda_link_l <- if (!is.null(fit_post$lambda)) fit_post$lambda[l] else 1
      
      muU_l <- as.numeric(X %*% betaU_l)
      zbeta <- Z %*% betaC_l
      
      p_cure <- switch(link,
                       "logit" = 1 / (1 + exp(-zbeta)),
                       "slogit" = (1 / (1 + exp(-zbeta)))^delta_l,
                       "rplogit" = 1 - (1 / (1 + exp(zbeta)))^delta_l,
                       "fglogit" = (1 - (1 / (1 + exp(zbeta)))^lambda_link_l)^delta_l,
                       stop("Unsupported link"))
      
      Su_l <- switch(dist,
                     "weibull" = {
                       alpha_l <- 1 / sigma_l
                       lambda_l <- exp(muU_l)
                       exp(-(time / lambda_l)^alpha_l)
                     },
                     "lognormal" = {
                       lambda_l <- muU_l
                       1 - plnorm(time, meanlog = lambda_l, sdlog = sigma_l)
                     },
                     "loglogistic" = {
                       shape <- 1 / sigma_l
                       scale <- exp(muU_l)
                       1 / (1 + (time / scale)^shape)
                     },
                     stop("Unsupported distribution"))
      
      S_mix_mat[, l] <- p_cure + (1 - p_cure) * Su_l
    }
    
    r_cs <- -log(rowMeans(S_mix_mat))
    list(residuals = r_cs, status = status)
  }
  
  # --- Step 1: Compute Cox-Snell residuals across simulations ---
  cs_all <- lapply(results, extract_cs_residual)
  all_r <- unlist(lapply(cs_all, `[[`, "residuals"))
  all_status <- unlist(lapply(cs_all, `[[`, "status`"))
  
  # --- Step 2: Pooled KM estimate ---
  km_all <- survfit(Surv(all_r, all_status) ~ 1)
  km_all_df <- data.frame(time = km_all$time, surv = km_all$surv)
  
  # --- Step 3: Posterior predictive envelope ---
  grid_times <- seq(0, max(all_r), length.out = 100)
  exp_curve <- data.frame(time = grid_times, surv = exp(-grid_times))
  
  km_surv_matrix <- matrix(NA, nrow = length(grid_times), ncol = length(cs_all))
  for (i in seq_along(cs_all)) {
    km_i <- survfit(Surv(cs_all[[i]]$residuals, cs_all[[i]]$status) ~ 1)
    interp <- approx(km_i$time, km_i$surv, xout = grid_times, method = "constant", rule = 2)$y
    km_surv_matrix[, i] <- interp
  }
  
  lower <- apply(km_surv_matrix, 1, quantile, probs = 0.025, na.rm = TRUE)
  upper <- apply(km_surv_matrix, 1, quantile, probs = 0.975, na.rm = TRUE)
  lower <- cummin(lower)
  upper <- cummin(upper)
  
  # --- Step 4: Plot ---
  if (!is.null(save_path)) {
    pdf(file.path(save_path, paste0("coxsnell_", tag, ".pdf")), width = 7, height = 5)
  }
  
  plot(km_all, xlab = "Cox-Snell Residual", ylab = "Survival",
       main = paste0("Pooled KM Cox-Snell Residuals with 95% Envelope (", tag, ")"),
       xlim = range(grid_times), ylim = c(0, 1))
  curve(exp(-x), add = TRUE, col = "red", lty = 2)
  lines(grid_times, lower, col = "gray", lty = 3)
  lines(grid_times, upper, col = "gray", lty = 3)
  
  legend("topright",
         legend = c("KM of Residuals", "Exponential(1)", "95% Envelope"),
         col = c("black", "red", "gray"), lty = c(1, 2, 3), bty = "n")
  
  if (!is.null(save_path)) dev.off()
  
  # --- Return results ---
  invisible(list(
    tag = tag,
    residuals = all_r,
    status = all_status,
    grid_times = grid_times,
    lower = lower,
    upper = upper,
    km_df = km_all_df,
    exp_curve = exp_curve
  ))
}


plot_cs_envelope(
  results = wbl.fit,
  data = NULL,  # kept for compatibility
  time = time,
  status = status,
  X = X,
  Z = Z,
  dist = "weibull",
  link = "logit",
  save_path = NULL
)



plot_cs_envelope_single <- function(fit,
                                    time,
                                    status,
                                    X,
                                    Z,
                                    dist = "weibull",
                                    link = "logit",
                                    tag = NULL,
                                    save_path = NULL,
                                    envelope = TRUE) {
  if (is.null(tag)) tag <- paste(dist, link, sep = "_")
  
  fit_post <- rstan::extract(fit)
  n <- length(time)
  L <- length(fit_post$sigma)
  S_mix_mat <- matrix(NA, n, L)
  
  for (l in 1:L) {
    betaU_l <- fit_post$betaU[l, ]
    betaC_l <- fit_post$betaC[l, ]
    sigma_l <- fit_post$sigma[l]
    
    delta_l <- if (!is.null(fit_post$delta)) fit_post$delta[l] else 1
    lambda_link_l <- if (!is.null(fit_post$lambda)) fit_post$lambda[l] else 1
    
    muU_l <- as.numeric(X %*% betaU_l)
    zbeta <- Z %*% betaC_l
    
    p_cure <- switch(link,
                     "logit" = 1 / (1 + exp(-zbeta)),
                     "slogit" = (1 / (1 + exp(-zbeta)))^delta_l,
                     "rplogit" = 1 - (1 / (1 + exp(zbeta)))^delta_l,
                     "fglogit" = (1 - (1 / (1 + exp(zbeta)))^lambda_link_l)^delta_l,
                     stop("Unsupported link"))
    
    Su_l <- switch(dist,
                   "weibull" = {
                     alpha_l <- 1 / sigma_l
                     lambda_l <- exp(muU_l)
                     exp(-(time / lambda_l)^alpha_l)
                   },
                   "lognormal" = {
                     lambda_l <- muU_l
                     1 - plnorm(time, meanlog = lambda_l, sdlog = sigma_l)
                   },
                   "loglogistic" = {
                     shape <- 1 / sigma_l
                     scale <- exp(muU_l)
                     1 / (1 + (time / scale)^shape)
                   },
                   stop("Unsupported distribution"))
    
    S_mix_mat[, l] <- p_cure + (1 - p_cure) * Su_l
  }
  
  r_cs_pp <- -log(rowMeans(S_mix_mat))
  km_fit <- survfit(Surv(r_cs_pp, status) ~ 1)
  km_df <- data.frame(time = km_fit$time, surv = km_fit$surv)
  
  grid_times <- seq(0, max(r_cs_pp), length.out = 100)
  exp_curve <- data.frame(time = grid_times, surv = exp(-grid_times))
  
  library(ggplot2)
  
  # Core KM data
  km_df <- data.frame(
    time = km_fit$time,
    surv = km_fit$surv,
    lower = km_fit$lower,
    upper = km_fit$upper
  )
  grid_times <- seq(0, max(r_cs_pp), length.out = 100)
  exp_curve <- data.frame(time = grid_times, surv = exp(-grid_times))
  
  # Envelope (if requested)
  if (envelope) {
    km_surv_matrix <- matrix(NA, nrow = length(grid_times), ncol = L)
    for (l in 1:L) {
      r_cs_l <- -log(S_mix_mat[, l])
      km_l <- survfit(Surv(r_cs_l, status) ~ 1)
      interp <- approx(km_l$time, km_l$surv, xout = grid_times, method = "linear", rule = 2)$y
      km_surv_matrix[, l] <- interp
    }
    
    lower <- apply(km_surv_matrix, 1, quantile, probs = 0.025, na.rm = TRUE)
    upper <- apply(km_surv_matrix, 1, quantile, probs = 0.975, na.rm = TRUE)
    envelope_df <- data.frame(time = grid_times, lower = lower, upper = upper)
  } else {
    envelope_df <- NULL
  }
  
  # Construct base ggplot
  p <- ggplot() +
    # Posterior envelope
    {if (!is.null(envelope_df)) geom_ribbon(data = envelope_df,
                                            aes(x = time, ymin = lower, ymax = upper, fill = "95% Envelope"),
                                            alpha = 0.5)} +
    # KM confidence interval (dotted lines)
    geom_step(data = km_df, aes(x = time, y = lower), color = "black", linetype = "dotted") +
    geom_step(data = km_df, aes(x = time, y = upper), color = "black", linetype = "dotted") +
    # KM curve
    geom_step(data = km_df, aes(x = time, y = surv, color = "KM of Residuals"), linewidth = 1) +
    # Exp(1) curve
    geom_line(data = exp_curve, aes(x = time, y = surv, color = "Exponential(1)"),
              linetype = "dashed", linewidth = 1) +
    labs(x = "Cox–Snell Residual", y = "Survival") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(name = NULL,
                       values = c("KM of Residuals" = "black",
                                  "Exponential(1)" = "red")) +
    scale_fill_manual(name = NULL,
                      values = c("95% Envelope" = "gray80")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top")
  
  
  # Save if needed
  if (!is.null(save_path)) {
    ggsave(filename = file.path(save_path, paste0("coxsnell_", tag, ".pdf")),
           plot = p, width = 7, height = 5)
  }
}



plot.wb1 <- plot_cs_envelope_single(
  fit = wbl.fit,
  time = iowa_rec18_dum$time_to_recurrence,
  status = iowa_rec18_dum$is_censored,
  X = X,
  Z = Z,
  dist = "weibull",
  link = "logit",
  envelope = TRUE, 
  save_path = "~/Research/Paper1/dataresult/", 
  tag = "weibull_logit"
)


plot.wb2 <- plot_cs_envelope_single(
  fit = wbsl.fit,
  time = iowa_rec18_dum$time_to_recurrence,
  status = iowa_rec18_dum$is_censored,
  X = X,
  Z = Z,
  dist = "weibull",
  link = "slogit",
  envelope = TRUE, 
  save_path = "~/Research/Paper1/dataresult/", 
  tag = "weibull_slogit"
)

plot.wb3 <- plot_cs_envelope_single(
  fit = wbrpl.fit,
  time = iowa_rec18_dum$time_to_recurrence,
  status = iowa_rec18_dum$is_censored,
  X = X,
  Z = Z,
  dist = "weibull",
  link = "rplogit",
  envelope = TRUE, 
  save_path = "~/Research/Paper1/dataresult/", 
  tag = "weibull_rplogit"
)


plot.wb4 <- plot_cs_envelope_single(
  fit = wbfg.fit,
  time = iowa_rec18_dum$time_to_recurrence,
  status = iowa_rec18_dum$is_censored,
  X = X,
  Z = Z,
  dist = "weibull",
  link = "fglogit",
  envelope = TRUE, 
  save_path = "~/Research/Paper1/dataresult/", 
  tag = "weibull_fglogit"
)


plot.ln1 <- plot_cs_envelope_single(
  fit = lnl.fit,
  time = iowa_rec18_dum$time_to_recurrence,
  status = iowa_rec18_dum$is_censored,
  X = X,
  Z = Z,
  dist = "lognormal",
  link = "logit",
  envelope = TRUE, 
  save_path = "~/Research/Paper1/dataresult/", 
  tag = "lognormal_logit"
)


plot.ln2 <- plot_cs_envelope_single(
  fit = lnsl.fit,
  time = iowa_rec18_dum$time_to_recurrence,
  status = iowa_rec18_dum$is_censored,
  X = X,
  Z = Z,
  dist = "lognormal",
  link = "slogit",
  envelope = TRUE, 
  save_path = "~/Research/Paper1/dataresult/", 
  tag = "lognormal_slogit"
)

plot.ln3 <- plot_cs_envelope_single(
  fit = lnrpl.fit,
  time = iowa_rec18_dum$time_to_recurrence,
  status = iowa_rec18_dum$is_censored,
  X = X,
  Z = Z,
  dist = "lognormal",
  link = "rplogit",
  envelope = TRUE, 
  save_path = "~/Research/Paper1/dataresult/", 
  tag = "lognormal_rplogit"
)


plot.ln4 <- plot_cs_envelope_single(
  fit = lnfg.fit,
  time = iowa_rec18_dum$time_to_recurrence,
  status = iowa_rec18_dum$is_censored,
  X = X,
  Z = Z,
  dist = "lognormal",
  link = "fglogit",
  envelope = TRUE, 
  save_path = "~/Research/Paper1/dataresult/", 
  tag = "lognormal_fglogit"
)



plot.llg1 <- plot_cs_envelope_single(
  fit = llgl.fit,
  time = iowa_rec18_dum$time_to_recurrence,
  status = iowa_rec18_dum$is_censored,
  X = X,
  Z = Z,
  dist = "loglogistic",
  link = "logit",
  envelope = TRUE, 
  save_path = "~/Research/Paper1/dataresult/", 
  tag = "loglogistic_logit"
)


plot.llg2 <- plot_cs_envelope_single(
  fit = llgsl.fit,
  time = iowa_rec18_dum$time_to_recurrence,
  status = iowa_rec18_dum$is_censored,
  X = X,
  Z = Z,
  dist = "loglogistic",
  link = "slogit",
  envelope = TRUE, 
  save_path = "~/Research/Paper1/dataresult/", 
  tag = "loglogistic_slogit"
)

plot.llg3 <- plot_cs_envelope_single(
  fit = llgrpl.fit,
  time = iowa_rec18_dum$time_to_recurrence,
  status = iowa_rec18_dum$is_censored,
  X = X,
  Z = Z,
  dist = "loglogistic",
  link = "rplogit",
  envelope = TRUE, 
  save_path = "~/Research/Paper1/dataresult/", 
  tag = "final"
)


plot.llg4 <- plot_cs_envelope_single(
  fit = llgfg.fit,
  time = iowa_rec18_dum$time_to_recurrence,
  status = iowa_rec18_dum$is_censored,
  X = X,
  Z = Z,
  dist = "loglogistic",
  link = "fglogit",
  envelope = TRUE, 
  save_path = "~/Research/Paper1/dataresult/", 
  tag = "loglogistic_fglogit"
)





compute_pp_value <- function(S_mix_mat, status, grid_times = seq(0, 2, length.out = 100)) {
  L <- ncol(S_mix_mat)
  T_rep <- numeric(L)
  
  # Compute T_obs
  r_obs <- -log(rowMeans(S_mix_mat))
  km_obs <- survfit(Surv(r_obs, status) ~ 1)
  s_obs <- approx(km_obs$time, km_obs$surv, xout = grid_times, method = "constant", rule = 2)$y
  T_obs <- sum((s_obs - exp(-grid_times))^2)
  
  # Replicates
  for (l in 1:L) {
    r_l <- -log(S_mix_mat[, l])
    km_l <- survfit(Surv(r_l, status) ~ 1)
    s_l <- approx(km_l$time, km_l$surv, xout = grid_times, method = "constant", rule = 2)$y
    T_rep[l] <- sum((s_l - exp(-grid_times))^2)
  }
  
  ppp <- mean(T_rep >= T_obs)
  return(list(ppp = ppp, T_obs = T_obs, T_rep = T_rep))
}


result.wb1 <- compute_pp_value(plot.wb1$S_mix_mat, iowa_rec18_dum$is_censored)
result.wb1$ppp  # Posterior predictive p-value
0.48675

result.wb2 <- compute_pp_value(plot.wb2$S_mix_mat, iowa_rec18_dum$is_censored)
result.wb2$ppp  # Posterior predictive p-value
0.39575

result.wb3 <- compute_pp_value(plot.wb3$S_mix_mat, iowa_rec18_dum$is_censored)
result.wb3$ppp  # Posterior predictive p-value
0.736

result.wb4 <- compute_pp_value(plot.wb4$S_mix_mat, iowa_rec18_dum$is_censored)
result.wb4$ppp  # Posterior predictive p-value
0.438


result.ln1 <- compute_pp_value(plot.ln1$S_mix_mat, iowa_rec18_dum$is_censored)
result.ln1$ppp  # Posterior predictive p-value
0.16675

result.ln2 <- compute_pp_value(plot.ln2$S_mix_mat, iowa_rec18_dum$is_censored)
result.ln2$ppp  # Posterior predictive p-value
0.18125

result.ln3 <- compute_pp_value(plot.ln3$S_mix_mat, iowa_rec18_dum$is_censored)
result.ln3$ppp  # Posterior predictive p-value
0.209

result.ln4 <- compute_pp_value(plot.ln4$S_mix_mat, iowa_rec18_dum$is_censored)
result.ln4$ppp  # Posterior predictive p-value
0.2



result.llg1 <- compute_pp_value(plot.llg1$S_mix_mat, iowa_rec18_dum$is_censored)
result.llg1$ppp  # Posterior predictive p-value
0.3835

result.llg2 <- compute_pp_value(plot.llg2$S_mix_mat, iowa_rec18_dum$is_censored)
result.llg2$ppp  # Posterior predictive p-value
0.2775

result.llg3 <- compute_pp_value(plot.llg3$S_mix_mat, iowa_rec18_dum$is_censored)
result.llg3$ppp  # Posterior predictive p-value
0.25725

result.llg4 <- compute_pp_value(plot.llg4$S_mix_mat, iowa_rec18_dum$is_censored)
result.llg4$ppp  # Posterior predictive p-value
0.2825