functions{
   real cen_lognormal_lpdf(real t, int is_censored, real muU, real sigma, row_vector z, vector betaC) {
      real log_lik_cen;
      real log_lik_uncen;
      real prop; 
      real prob = 0;
      prop = inv_logit(z * betaC); 
    
      log_lik_cen = (1-is_censored) * log(prop + (1-prop) * exp(lognormal_lccdf(t|muU, sigma)));
      log_lik_uncen = is_censored * log1m(prop) + is_censored * lognormal_lpdf(t|muU, sigma);

      prob = log_lik_cen + log_lik_uncen ; 
      return(prob) ;
  }
} 

//data block
data{
  int N;  // number of observations
  vector <lower = 0> [N] t; // observed times
  array[N] int <lower = 0, upper = 1> is_censored; //censoring (0 = censored., 1 = observed)
  int M; // number of covariates
  matrix[N, M] x; //matrix of covariates (N rows, M columns)
  
  int K;
  matrix[N, K] z; 
}

//parameters block
parameters{
  vector [M] betaU; //coeff. in the linear prediction for survival function
  vector [K] betaC; // coeff. in the linear prediction for logistic regression
  real<lower = 0> sigma; // scale parameter sigma
}
//transformed parameters block
transformed parameters{
  vector[N] muU;
  
  muU = x * betaU ; // linear predictor
  
}
// model block
model{
  sigma ~ cauchy(0, 2.5);
  betaU ~ normal(0, 1);
  betaC ~ normal(0, 1);
  
  for (i in 1:N) {
    real log_lik_t = cen_lognormal_lpdf(t[i]| is_censored[i], muU[i], sigma, z[i], betaC);
    target += log_lik_t;
  }
}
// generated quantities block
generated quantities{
  vector[N] log_lik; 
  for(i in 1:N){
    log_lik[i] = cen_lognormal_lpdf(t[i]|is_censored[i], muU[i], sigma, z[i], betaC);
    
  }
  
}
