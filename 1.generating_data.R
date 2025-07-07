library(rstan)
library(tidyverse)
library(ggplot2)
library(xtable)
library(tictoc)
library(bayesplot)
library(loo)



# Generalized Weibull distribution

gen_weibull <- function(n, alpha, lambda) {
  # Generate n random numbers from the Generalized Weibull distribution
  U <- runif(n)
  T <- (-(log(1-U))/lambda)^(1/alpha)
  return(T)
}

weibull_pdf <- function(x, alpha, lambda) {
  pdf <- alpha * lambda * x^(alpha - 1) * exp(- lambda* x^alpha)
  return(pdf)
}

weibull_sf <- function(x, alpha, lambda) {
  sf <- exp(-lambda*x^alpha)
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



logLik.Weibull <- function(cc, x, alpha, lambda){
  logLik <- sum(cc*log(weibull_pdf(x, alpha, lambda)))+sum((1-cc)*log(weibull_sf(x, alpha, lambda)))
  return(logLik)
}


logLik.MCM.Weibull <- function(cc, uncured, x, alpha, lambda, theta){
  f1<-f2<-rep(0, length(cc))
  f1 <- ifelse(log(theta*weibull_sf(x, alpha, lambda))==-Inf, 0, log(theta*weibull_sf(x, alpha, lambda)))
  f2 <- ifelse(log(theta*cc*uncured*weibull_pdf(x, alpha, lambda))==-Inf, 0, log(theta*cc*uncured*weibull_pdf(x, alpha, lambda)))
  
  logLik <- sum(f1)+
    sum((1-cc)*uncured*(f2))+
    sum((1-cc)*(1-uncured)*log(1-theta))
  return(logLik)
}




# Generate random numbers from the Lognormal distribution
gen_lognormal <- function(n, meanlog, sdlog) {
  # Generate n random numbers from the Lognormal distribution
  T <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)
  return(T)
}

# Lognormal Probability Density Function (PDF)
lognormal_pdf <- function(x, meanlog, sdlog) {
  pdf <- (1 / (x * sdlog * sqrt(2 * pi))) * exp(-((log(x) - meanlog)^2) / (2 * sdlog^2))
  return(pdf)
}

# Lognormal Survival Function (SF)
lognormal_sf <- function(x, meanlog, sdlog) {
  sf <- 1 - pnorm(log(x), mean = meanlog, sd = sdlog)
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



logLik.Lognormal <- function(cc, x, meanlog, sdlog){
  logLik <- sum(cc*log(lognormal_pdf(x, meanlog, sdlog)))+sum((1-cc)*log(lognormal_sf(x, meanlog, sdlog)))
  return(logLik)
}


logLik.MCM.Lognormal <- function(cc, uncured, x, meanlog, sdlog, theta){
  f1<-f2<-rep(0, length(cc))
  f1 <- ifelse(log(theta*lognormal_sf(x, meanlog, sdlog))==-Inf, 0, log(theta*lognormal_sf(x, meanlog, sdlog)))
  f2 <- ifelse(log(theta*cc*uncured*lognormal_pdf(x, meanlog, sdlog))==-Inf, 0, log(theta*cc*uncured*lognormal_pdf(x, meanlog, sdlog)))
  
  logLik <- sum(f1)+
    sum((1-cc)*uncured*(f2))+
    sum((1-cc)*(1-uncured)*log(1-theta))
  return(logLik)
}




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

slogit <- function(tau = 1) {
  
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
  link <- list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, valideta = valideta, name = paste("slogit(tau=", tau, ")", sep=""))
  class(link) <- "link-glm"
  return(link)
}


### generating Weibull model ##############################################################################################################

# tau stands for alpha_s and lambda_rp, alpha_fg respectively for each link functions slogit, rplogit, and fglogit in paper
# delta stands for lambda_fg for the fglogit link in paper 

geraMCM.WB <- function(n, x, w, censoring_time=10 , beta=c(1,2), eta=c(1,2), alpha=1.5, link="logit", tau=1, delta = 1){
  #### tau=skewness parameter of power links
  if(link=="logit"){p_cure <- 1 / (1 + exp(-w%*%eta))}
  if(link=="probit"){p_cure <- pnorm(w%*%eta)} 
  if(link=="slogit"){p_cure<-(1/(1+exp(-w %*% eta)))^tau}
  if(link == "rplogit"){p_cure <- (1 - 1/(1 + exp(w%*%eta))^tau)}
  if(link == "fglogit"){p_cure <- (1 - (1/(1+exp(w%*%eta)))^tau)^(delta)}
  cured <- rbinom(n, 1, p_cure)     
  lambda <- exp(x%*%beta) 
  event_time <- rep(NA, n)
  non_cured <- which(cured == 0) 
  u <- runif(length(non_cured))  # uniform random variables for inverse transform sampling
  event_time[non_cured] <- (- log(1 - u) )^(1 / alpha) * lambda[non_cured]
  censor_time <- runif(n, 0, censoring_time)
  
  observed_time <- pmin(event_time, censor_time, na.rm = TRUE)
  status <- ifelse(cured == 1, 0, ifelse(event_time <= censor_time, 1, 0))  # 0 = censored, 1 = event
  return(list(x = x, w= w, cured = cured, lambda = lambda, event_time = event_time, censor_time = censor_time,
              observed_time = observed_time, status = status))  
}





### generating Lognormal model ##############################################################################################################
geraMCM.Lognormal <- function(n, x, w, censoring_time = 10, beta = c(1, 2), eta = c(1, 2), 
                              alpha = 1, link = "logit", tau = 1, delta = 1) {
  #### tau=skewness parameter of power links
  if (link == "logit") {p_cure <- 1 / (1 + exp(-w %*% eta))}
  if (link == "probit") {p_cure <- pnorm(w %*% eta)}
  if (link == "slogit") {p_cure <- (1 / (1 + exp(-w %*% eta)))^tau}
  if(link == "rplogit"){p_cure <- (1 - 1/(1 + exp(w%*%eta))^tau)}
  if(link == "fglogit"){p_cure <- (1 - (1/(1+exp(w%*%eta)))^tau)^(delta)}
  cured <- rbinom(n, 1, p_cure)
  lambda <- exp(x %*% beta)
  event_time <- rep(NA, n)
  non_cured <- which(cured == 0)
  # Generate event times for non-cured individuals using Lognormal distribution
  u <- runif(length(non_cured))  # uniform random variables for inverse transform sampling
  event_time[non_cured] <- exp(qnorm(u)*alpha + log(lambda[non_cured]))
  
  
  # Generate random censoring times
  censor_time <- runif(n, 0, censoring_time)
  # Compute observed time and censoring status
  observed_time <- pmin(event_time, censor_time, na.rm = TRUE)
  status <- ifelse(cured == 1, 0, ifelse(event_time <= censor_time, 1, 0))  # 0 = censored, 1 = event
  
  return(list(x = x, w = w, cured = cured, lambda = lambda, event_time = event_time, 
              censor_time = censor_time, observed_time = observed_time, status = status))
}




### generating Mixture Cure Model with loglogistic base

geraMCM.Loglogistic <- function(n, x, w, censoring_time=10 , beta=c(1,2), eta=c(1,2), alpha=1.5, link="logit", tau=1, delta = 1){
  #### tau=skewness parameter of power links
  if (link == "logit") {p_cure <- 1 / (1 + exp(-w %*% eta))}
  if (link == "probit") {p_cure <- pnorm(w %*% eta)}
  if (link == "slogit") {p_cure <- (1 / (1 + exp(-w %*% eta)))^tau}
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


