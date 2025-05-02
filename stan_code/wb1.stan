functions{

  real cen_weibull_lpdf(real t, int cured, int is_censored, real shape, real scale, row_vector z, vector betaC) {
    
    real prop = inv_logit(dot_product(z, betaC)); // Cure probability

    real log_lik_cen = (1 - is_censored) * log(prop + (1 - prop) * exp(weibull_lccdf(t | shape, scale))); 
    real log_lik_uncen = is_censored * log1m(prop) + is_censored * weibull_lpdf(t | shape, scale);

    return log_lik_cen + log_lik_uncen;
  }

}


//data block
data{
  int N;  // number of observations
  vector <lower = 0> [N] t; // observed times
  array[N] int <lower = 0, upper = 1> is_censored; //censoring (0 = censored., 1 = observed)
  array[N] int <lower = 0, upper = 1>  cured; // cured = 1, uncured = 0 
  int M; // number of covariates
  matrix[N, M] x; //matrix of covariates for survival(N rows, M columns)
  
  int K;
  matrix[N, K] z; //matrix of covariates for incidence 
}
//parameters block
parameters{
  vector [M] betaU; //coeff. in the linear prediction for survival function
  vector [K] betaC; // coeff. in the linear prediction for logistic regression
  real <lower = 0> sigma; 
  // vector <lower = 0, upper = 1> [N] cured; //cured(1 = cured, 0 = no cured)
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
   
  for (i in 1:N) {
    real log_lik_t = cen_weibull_lpdf(t[i]| cured[i], is_censored[i], 1/sigma, muU[i], z[i], betaC);
    target += cen_weibull_lpdf(t[i]| cured[i], is_censored[i], 1/sigma, muU[i], z[i], betaC);
  }
  
}

generated quantities{
  vector[N] log_lik;

  for(i in 1:N){
    log_lik[i] = cen_weibull_lpdf(t[i]| cured[i], is_censored[i], 1/sigma, muU[i], z[i], betaC);
  }
}
