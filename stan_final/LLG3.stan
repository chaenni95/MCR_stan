functions{

  real cen_loglogistic_lpdf( real t, real is_censored, real scale, real shape, row_vector z, vector betaC,  real lambda_rp) {
   real log_lik_cen;
    real log_lik_uncen;
    real prop; 
    real prob = 0;
      
      prop = 1-pow(1-inv_logit(z * betaC), lambda_rp); 
    
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
  real <lower = 0> lambda_rp;
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
  lambda_rp ~ lognormal(0, 1); 
   
  for (i in 1:N) {
    real log_lik_t = cen_loglogistic_lpdf(t[i]|is_censored[i], muU[i], 1/sigma, z[i], betaC, lambda_rp);
    target += log_lik_t;
  }
  
}

generated quantities{
  vector[N] log_lik; 
  for (i in 1:N) {
    log_lik[i] = cen_loglogistic_lpdf(t[i]|is_censored[i], muU[i], 1/sigma, z[i], betaC, lambda_rp);
    
  }
}
  
