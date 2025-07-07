library(rstan)
library(tidyverse)
library(ggplot2)
library(xtable)
library(tictoc)
library(bayesplot)
library(loo)

stanmodel.wb1 <- stan_model("~/Recidivism/simulation/wb1.stan")
stanmodel.wb2 <- stan_model("~/Recidivism/simulation/wb2.stan")
stanmodel.wb3 <- stan_model("~/Recidivism/simulation/wb3.stan")
stanmodel.wb4 <- stan_model("~/Recidivism/simulation/wb4.stan")

stanmodel.ln1 <- stan_model("~/Recidivism/simulation/LN1.stan")
stanmodel.ln2 <- stan_model("~/Recidivism/simulation/LN2.stan")
stanmodel.ln3 <- stan_model("~/Recidivism/simulation/LN3.stan")
stanmodel.ln4 <- stan_model("~/Recidivism/simulation/LN4.stan")

stanmodel.llg1 <- stan_model("~/Recidivism/simulation/LLG1.stan")
stanmodel.llg2 <- stan_model("~/Recidivism/simulation/LLG2.stan")
stanmodel.llg3 <- stan_model("~/Recidivism/simulation/LLG3.stan")
stanmodel.llg4 <- stan_model("~/Recidivism/simulation/LLG4.stan")




# make stan data out of data 
generate_stan_data <- function(data) {
  list(
    N = length(data$observed_time), 
    t = data$observed_time, 
    is_censored = data$status,
    x = data$x, 
    z = data$w, 
    M = ncol(data$x), 
    K = ncol(data$w)
  )
}


### single simulation cox-snell residual KM plot 
cox_snell_pp_residuals <- function(data, post, dist = "weibull", link = "logit",
                                   envelope = FALSE, filename = NULL) {
  data_name <- deparse(substitute(data))
  n <- length(data$observed_time)
  L <- length(post$sigma)
  S_mix_mat <- matrix(NA, n, L)
  
  for (l in 1:L) {
    betaU_l <- post$betaU[l, ]
    betaC_l <- post$betaC[l, ]
    sigma_l <- post$sigma[l]
    
    # Optional link parameters
    delta_l <- if (!is.null(post$delta)) post$delta[l] else 1
    lambda_link_l <- if (!is.null(post$lambda)) post$lambda[l] else 1
    
    muU_l <- as.numeric(data$x %*% betaU_l)
    zbeta <- data$w %*% betaC_l
    
    # Cure probabilities
    p_cure <- switch(link,
                     "logit" = 1 / (1 + exp(-zbeta)),
                     "slogit" = (1 / (1 + exp(-zbeta)))^delta_l,
                     "rplogit" = 1 - (1 / (1 + exp(zbeta)))^delta_l,
                     "fglogit" = (1 - (1 / (1 + exp(zbeta)))^lambda_link_l)^delta_l,  
                     stop("Unsupported link function")
    )
    
    # Uncured survival
    Su_l <- switch(dist,
                   "weibull" = {
                     alpha_l <- 1 / sigma_l
                     lambda_l <- exp(muU_l)
                     exp(-(data$observed_time / lambda_l)^alpha_l)
                   },
                   "lognormal" = {
                     lambda_l <- muU_l
                     1 - plnorm(data$observed_time, meanlog = lambda_l, sdlog = sigma_l)
                   },
                   "loglogistic" = {
                     shape <- 1 / sigma_l
                     scale <- exp(muU_l)
                     1 / (1 + (data$observed_time / scale)^shape)
                   },
                   stop("Unsupported distribution")
    )
    
    S_mix_mat[, l] <- p_cure + (1 - p_cure) * Su_l
  }
  
  S_mix_mean <- rowMeans(S_mix_mat)
  r_cs_pp <- -log(S_mix_mean)
  
  library(survival)
  library(survminer)
  
  surv_obj <- Surv(r_cs_pp, data$status)
  km_fit <- survfit(surv_obj ~ 1)
  
  plot_data <- data.frame(
    time = km_fit$time,
    surv = km_fit$surv,
    lower = km_fit$lower,
    upper = km_fit$upper
  )
  
  p <- ggplot(plot_data, aes(x = time, y = surv)) +
    # KM estimate
    geom_step(aes(linetype = "KM"), color = "black", size = 0.6) +
    # KM 95% CI
    geom_step(aes(y = lower, linetype = "KM 95% CI"), color = "black", size = 0.5) +
    geom_step(aes(y = upper, linetype = "KM 95% CI"), color = "black", size = 0.5) +
    # Exponential(1) line
    stat_function(fun = function(x) exp(-x), aes(linetype = "Exponential(1)"), color = "red", size = 0.6) +
    scale_linetype_manual(
      name = NULL,
      values = c(
        "KM" = "solid",
        "KM 95% CI" = "dotted",
        "Exponential(1)" = "dashed"
      )
    ) +
    labs(
      x = "Cox–Snell Residual",
      y = "Survival"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom", legend.box = "horizontal")
  
  if (envelope) {
    times <- seq(0, max(r_cs_pp), length.out = 100)
    surv_curves <- matrix(NA, nrow = length(times), ncol = L)
    
    for (l in 1:L) {
      r_cs_l <- -log(S_mix_mat[, l])
      km_l <- survfit(Surv(r_cs_l, data$status) ~ 1)
      # interpolate at fixed times
      surv_interp <- approx(km_l$time, km_l$surv, xout = times, method = "linear", rule = 2)$y
      surv_curves[, l] <- surv_interp
    }
    
    lower <- apply(surv_curves, 1, quantile, probs = 0.025, na.rm = TRUE)
    upper <- apply(surv_curves, 1, quantile, probs = 0.975, na.rm = TRUE)
    envelope_df <- data.frame(time = times, lower = lower, upper = upper)
    
    # Add shaded envelope correctly
    
    p <- p + geom_ribbon(
      data = envelope_df,
      mapping = aes(x = time, ymin = lower, ymax = upper, fill = "Posterior Envelope"),
      inherit.aes = FALSE,
      alpha = 0.5
    ) +
      scale_fill_manual(values = c("Posterior Envelope" = "gray80"), guide = guide_legend(override.aes = list(alpha = 0.5))) +
      guides(linetype = guide_legend(order = 1), fill = guide_legend(order = 2))
  }
  
  if (!is.null(filename)) ggsave(filename, plot = p, width = 8, height = 6)
  
  return(list(residuals = r_cs_pp, km = km_fit, plot = p))
}


#### simulation example with Weibull distribution ###############################################################################################
n <- 1000
alpha <- 1.5
beta <- c(1,-1,-2)
eta <- c(0.5,0.6)
tau <- 1.5  
delta <- 0.5

p <- length(beta)
q <- length(eta)

x <- matrix(rnorm(n*(p-1)),n,p-1)
x <- cbind(1,x)

w <- scale(matrix(runif(n*(q),-1,1),n,q)) # no intercept for logit!

censoring_time <- 100



data.wb1.1000 <- geraMCM.WB(n, x, w, censoring_time=censoring_time , 
                            beta=beta, eta=eta, alpha=alpha,link="logit", tau=1, delta = 1)
data.wb2.1000 <- geraMCM.WB(n, x, w, censoring_time=censoring_time , 
                            beta=beta, eta=eta, alpha=alpha,link="slogit", tau=1.5, delta = 1)
data.wb3.1000 <- geraMCM.WB(n, x, w, censoring_time=censoring_time , 
                            beta=beta, eta=eta, alpha=alpha,link="rplogit", tau=1.5, delta = 1)
data.wb4.1000 <- geraMCM.WB(n, x, w, censoring_time=censoring_time , 
                            beta=beta, eta=eta, alpha=alpha,link="fglogit", tau=1.5, delta = 0.5)



### true logit ################################################################################
# fitted with logit
gwb1 <- generate_stan_data(data.wb1.1000)
fit.wb1.1000 <- sampling(stanmodel.wb1, data = gwb1, seed = 12345, chains = 1, iter = 2000, 
                         init = 0)

gwb1.post <- rstan::extract(fit.wb1.1000)


# fitted with slogit
fit.wb1.1000.s <- sampling(stanmodel.wb2, data = gwb1, seed = 12345, chains = 1, iter = 2000, 
                           init = 0)
gwb1.post.s <- rstan::extract(fit.wb1.1000.s)


# fitted with rplogit
fit.wb1.1000.rp <- sampling(stanmodel.wb3, data = gwb1, seed = 12345, chains = 1, iter = 2000, 
                            init = 0)
gwb1.post.rp <- rstan::extract(fit.wb1.1000.rp)


# fitted with fglogit
fit.wb1.1000.fg <- sampling(stanmodel.wb4, data = gwb1, seed = 12345, chains = 1, iter = 2000, 
                            init = 0)
gwb1.post.fg <- rstan::extract(fit.wb1.1000.fg)


### true slogit ##################################################################################

# logit
fit.wb2.1000.l <- sampling(stanmodel.wb1, data = gwb2, seed = 12345, chains = 1, iter = 2000, 
                           init = 0)
gwb2.post.l <- rstan::extract(fit.wb2.1000.l)

# slogit (true)
gwb2 <- generate_stan_data(data.wb2.1000)
fit.wb2.1000 <- sampling(stanmodel.wb2, data = gwb2, seed = 12345, chains = 1, iter = 2000, 
                         init = 0)

gwb2.post <- rstan::extract(fit.wb2.1000)


# rplogit
fit.wb2.1000.rp <- sampling(stanmodel.wb3, data = gwb2, seed = 12345, chains = 1, iter = 2000, 
                            init = 0)
gwb2.post.rp <- rstan::extract(fit.wb2.1000.rp)


# fglogit
fit.wb2.1000.fg <- sampling(stanmodel.wb4, data = gwb2, seed = 12345, chains = 1, iter = 2000, 
                            init = 0)
gwb2.post.fg <- rstan::extract(fit.wb2.1000.fg)


### true rplogit ##################################################################################
gwb3 <- generate_stan_data(data.wb3.1000)
# logit
fit.wb3.1000.l <- sampling(stanmodel.wb1, data = gwb3, seed = 12345, chains = 1, iter = 2000, 
                           init = 0)
gwb3.post.l <- rstan::extract(fit.wb3.1000.l)


# slogit
fit.wb3.1000.s <- sampling(stanmodel.wb2, data = gwb3, seed = 12345, chains = 1, iter = 2000, 
                           init = 0)
gwb3.post.s <- rstan::extract(fit.wb3.1000.s)

#rplogit (true)

fit.wb3.1000 <- sampling(stanmodel.wb3, data = gwb3, seed = 12345, chains = 1, iter = 2000, 
                         init = 0)

gwb3.post <- rstan::extract(fit.wb3.1000)

# fglogit
fit.wb3.1000.fg <- sampling(stanmodel.wb4, data = gwb3, seed = 12345, chains = 1, iter = 2000, 
                            init = 0)
gwb3.post.fg <- rstan::extract(fit.wb3.1000.fg)




### true fglogit ##################################################################################

# logit
fit.wb4.1000.l <- sampling(stanmodel.wb1, data = gwb4, seed = 12345, chains = 1, iter = 2000, 
                           init = 0)
gwb4.post.l <- rstan::extract(fit.wb4.1000.l)


# slogit
fit.wb4.1000.s <- sampling(stanmodel.wb2, data = gwb4, seed = 12345, chains = 1, iter = 2000, 
                           init = 0)
gwb4.post.s <- rstan::extract(fit.wb4.1000.s)

#rplogit
fit.wb4.1000.rp <- sampling(stanmodel.wb3, data = gwb4, seed = 12345, chains = 1, iter = 2000, 
                            init = 0)

gwb4.post.rp <- rstan::extract(fit.wb4.1000.rp)

# fglogit (true)
gwb4 <- generate_stan_data(data.wb4.1000)
fit.wb4.1000 <- sampling(stanmodel.wb4, data = gwb4, seed = 12345, chains = 1, iter = 2000, 
                         init = 0)

gwb4.post <- rstan::extract(fit.wb4.1000)



### lognormal fglogit fit to weibull fglogit
fit.wb4.1000.lnfg <- sampling(stanmodel.ln4, data = gwb4, seed = 12345, chains = 1, iter = 2000, 
                              init = 0)

gwb4.post.lnfg <- rstan::extract(fit.wb4.1000.lnfg)


### loglogistic fglogit fit to weibull fglogit
fit.wb4.1000.llgfg <- sampling(stanmodel.llg4, data = gwb4, seed = 12345, chains = 1, iter = 2000, 
                               init = 0)

gwb4.post.llgfg <- rstan::extract(fit.wb4.1000.llgfg)




####### cox-snell residual km plot example ##############################################################

### true model is Weibull logit 
cs.wbl.l <- cox_snell_pp_residuals(data.wb1.1000, gwb1.post, dist = "weibull", link = "logit", envelope = TRUE) 
cs.wbl.s <- cox_snell_pp_residuals(data.wb1.1000, gwb1.post.s, dist = "weibull", link = "slogit", envelope = TRUE) 
cs.wbl.rp <- cox_snell_pp_residuals(data.wb1.1000, gwb1.post.rp, dist = "weibull", link = "rplogit", envelope = TRUE)
cs.wbl.fg <- cox_snell_pp_residuals(data.wb1.1000, gwb1.post.fg, dist = "weibull", link = "fglogit", envelope = TRUE)
# Arrange plots
ggarrange(
  cs.wbl.l$plot + theme(legend.text = element_text(size = 14)),
  cs.wbl.s$plot + theme(legend.text = element_text(size = 14)),
  cs.wbl.rp$plot + theme(legend.text = element_text(size = 14)),
  cs.wbl.fg$plot + theme(legend.text = element_text(size = 14)),
  ncol = 2, nrow = 2,
  common.legend = TRUE,
  labels = c("a logit", "b slogit", "c rplogit", "d fglogit"),
  label.y = 1.03  # increase to move labels higher (default is 1)
)

### true model is Weibull slogit
cs.wbs.l <- cox_snell_pp_residuals(data.wb2.1000, gwb2.post.l, dist = "weibull", link = "logit", envelope = TRUE)
cs.wbs.s <- cox_snell_pp_residuals(data.wb2.1000, gwb2.post, dist = "weibull", link = "slogit", envelope = TRUE)
cs.wbs.rp <- cox_snell_pp_residuals(data.wb2.1000, gwb2.post.rp, dist = "weibull", link = "rplogit", envelope = TRUE)
cs.wbs.fg <- cox_snell_pp_residuals(data.wb2.1000, gwb2.post.fg, dist = "weibull", link = "fglogit", envelope = TRUE)
ggarrange(
  cs.wbs.l$plot + theme(legend.text = element_text(size = 14)),
  cs.wbs.s$plot + theme(legend.text = element_text(size = 14)),
  cs.wbs.rp$plot + theme(legend.text = element_text(size = 14)),
  cs.wbs.fg$plot + theme(legend.text = element_text(size = 14)),
  ncol = 2, nrow = 2,
  common.legend = TRUE,
  labels = c("logit", "slogit", "rplogit", "fglogit"),
  label.y = 1.03  # increase to move labels higher (default is 1)
)


### true model is Weibull rplogit
cs.wbrp.l <- cox_snell_pp_residuals(data.wb3.1000, gwb3.post.l, dist = "weibull", link = "logit", envelope = TRUE)
cs.wbrp.s <- cox_snell_pp_residuals(data.wb3.1000, gwb3.post.s, dist = "weibull", link = "slogit", envelope = TRUE)
cs.wbrp.rp <- cox_snell_pp_residuals(data.wb3.1000, gwb3.post, dist = "weibull", link = "rplogit", envelope = TRUE)
cs.wbrp.fg <- cox_snell_pp_residuals(data.wb3.1000, gwb3.post.fg, dist = "weibull", link = "fglogit", envelope = TRUE)
ggarrange(
  cs.wbrp.l$plot + theme(legend.text = element_text(size = 14)),
  cs.wbrp.s$plot + theme(legend.text = element_text(size = 14)),
  cs.wbrp.rp$plot + theme(legend.text = element_text(size = 14)),
  cs.wbrp.fg$plot + theme(legend.text = element_text(size = 14)),
  ncol = 2, nrow = 2,
  common.legend = TRUE,
  labels = c("logit", "slogit", "rplogit", "fglogit"),
  label.y = 1.03  # increase to move labels higher (default is 1)
)


### true model is Weibull fglogit
cs.wbfg.l <- cox_snell_pp_residuals(data.wb4.1000, gwb4.post.l, dist = "weibull", link = "logit", envelope = TRUE)
cs.wbfg.s <- cox_snell_pp_residuals(data.wb4.1000, gwb4.post.s, dist = "weibull", link = "slogit", envelope = TRUE)
cs.wbfg.rp <- cox_snell_pp_residuals(data.wb4.1000, gwb4.post.rp, dist = "weibull", link = "rplogit", envelope = TRUE)
cs.wbfg.fg <- cox_snell_pp_residuals(data.wb4.1000, gwb4.post, dist = "weibull", link = "fglogit", envelope = TRUE)
ggarrange(
  cs.wbfg.l$plot + theme(legend.text = element_text(size = 14)),
  cs.wbfg.s$plot + theme(legend.text = element_text(size = 14)),
  cs.wbfg.rp$plot + theme(legend.text = element_text(size = 14)),
  cs.wbfg.fg$plot + theme(legend.text = element_text(size = 14)),
  ncol = 2, nrow = 2,
  common.legend = TRUE,
  labels = c("logit", "slogit", "rplogit", "fglogit"),
  label.y = 1.03  # increase to move labels higher (default is 1)
)


### true model is Weibull fglogit but fitted with different distributions 
cs.wbfg.lnfg <- cox_snell_pp_residuals(data.wb4.1000, gwb4.post.lnfg, dist = "lognormal", link = "fglogit", envelope = TRUE)
cs.wbfg.llgfg <- cox_snell_pp_residuals(data.wb4.1000, gwb4.post.llgfg, dist = "loglogistic", link = "fglogit", envelope = TRUE)
cs.wbfg.fg <- cox_snell_pp_residuals(data.wb4.1000, gwb4.post, dist = "weibull", link = "fglogit", envelope = TRUE)
ggarrange(
  cs.wbfg.lnfg$plot+ theme(legend.text = element_text(size = 14)),
  cs.wbfg.llgfg$plot+ theme(legend.text = element_text(size = 14)),
  cs.wbfg.fg$plot+ theme(legend.text = element_text(size = 14)), # true
  ncol = 2, nrow = 2,
  common.legend = TRUE,
  labels = c("lognormal", "loglogistic", "weibull"),
  label.y = 1.03  # increase to move labels higher (default is 1)
)



### logit model fit to other link functions 
cs.wbl.l <- cox_snell_pp_residuals(data.wb1.1000, gwb1.post, dist = "weibull", link = "logit", envelope = TRUE)
cs.wbl.s <- cox_snell_pp_residuals(data.wb2.1000, gwb2.post.l, dist = "weibull", link = "logit", envelope = TRUE)
cs.wbl.rp <- cox_snell_pp_residuals(data.wb3.1000, gwb3.post.l, dist = "weibull", link = "logit", envelope = TRUE)
cs.wbl.fg <- cox_snell_pp_residuals(data.wb4.1000, gwb4.post.l, dist = "weibull", link = "logit", envelope = TRUE)
ggarrange(
  cs.wbl.l$plot + theme(legend.text = element_text(size = 14)),
  cs.wbl.s$plot + theme(legend.text = element_text(size = 14)),
  cs.wbl.rp$plot + theme(legend.text = element_text(size = 14)),
  cs.wbl.fg$plot + theme(legend.text = element_text(size = 14)),
  ncol = 2, nrow = 2,
  common.legend = TRUE,
  labels = c("a logit", "b slogit", "c rplogit", "d fglogit"),
  label.y = 1.03  # increase to move labels higher (default is 1)
)


