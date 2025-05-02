functions{
  
  // real loglogistic_lpdf (real y, real alpha, real beta){
  //   real lalpha = log(alpha);
  //   real numerator = log(beta) - lalpha + (beta - 1) * (log(y) - lalpha);
  //   return numerator - 2 * (log1p( (y/alpha)^beta ));
  // }
  
  // real loglogistic_lcdf (real y, real alpha, real beta) {
  //   return -log1p((y / alpha) ^-beta);
  // }
  
  // real loglogistic_cdf (real y, real alpha, real beta) {
  //   return 1/(1 + (y / alpha) ^-beta);
  // }
  
  // real loglogistic_lccdf (real y, real alpha, real beta){
  //   return -log1p( (y / alpha)^beta);
  // }
  // 
  // real loglogistic_rng (real alpha, real beta, real lb, real ub) {
  //   real p_ub = loglogistic_cdf(ub, alpha, beta);
  //   real p_lb = loglogistic_cdf(lb, alpha, beta);
  //   real r = uniform_rng(p_lb, p_ub);
  //   return alpha * (r / (1 - r) )^(1/beta);
  // }
  
 real lik_bern_lpmf(int cured, row_vector z, vector betaC, real delta){
    real lik_b;
    real curedd = cured/1.0; 
    
    lik_b = curedd * log(pow(inv_logit(z * betaC), delta)) + 
        (1-curedd) * log(1-pow(inv_logit(z * betaC), delta));
    
    return lik_b;
  }
  

  real cen_loglogistic_lpdf( real t, int cured, int is_censored, real scale, real shape, row_vector z, vector betaC,  real delta) {
   real log_lik_cen;
    real log_lik_uncen;
    real prop; 
    real prob = 0;
      
    real lik_b = lik_bern_lpmf(cured| z, betaC, delta); 
      prop = pow(inv_logit(z * betaC), delta); 
    
    log_lik_cen = (1-is_censored) * log(prop + (1-prop) * (1-loglogistic_cdf(t|scale, shape)));
    log_lik_uncen = is_censored * log1m(prop) + is_censored * loglogistic_lpdf(t|scale, shape);

    prob = log_lik_cen + log_lik_uncen ; 
    return(prob) ;
  }
  
}
//data block
data{
  int N;  // number of observations
  vector <lower = 0> [N] t; // observed times
  array [N] int<lower = 0, upper = 1>  is_censored; //censoring (0 = censored., 1 = observed)
  array [N] int<lower = 0, upper = 1>  cured;
  int M; // number of covariates
  matrix[N, M] x; //matrix of covariates (N rows, M columns)
  
  int K;
  matrix[N, K] z; 
}

//parameters block
parameters{
  vector [M] betaU; //coeff. in the linear prediction for survival function
  vector [K] betaC; // coeff. in the linear prediction for logistic regression
  real <lower = 0> sigma; // scale parameter sigma = 1/shape
  real <lower = 0> delta;
}

//transformed parameters block
transformed parameters{
  vector[N] muU;
  
  muU = exp(x * betaU);

}

// model block
model{
  sigma ~ cauchy(0, 2.5);
  betaU ~ normal(0, 1);
  betaC ~ normal(0, 1);
  delta ~ lognormal(0, 1); 
   
  for (i in 1:N) {
    real log_lik_t = cen_loglogistic_lpdf(t[i]| cured[i], is_censored[i], muU[i], 1/sigma, z[i], betaC, delta);
    target += log_lik_t;
  }
  
}

generated quantities{
  vector[N] log_lik; 
  real curedd; 
  
  for (i in 1:N) {
    log_lik[i] = cen_loglogistic_lpdf(t[i]| cured[i], is_censored[i], muU[i], 1/sigma,  z[i], betaC, delta);
    }
}
  
