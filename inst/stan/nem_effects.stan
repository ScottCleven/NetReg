data{
  int<lower=0> N;          // Number of observations
  int<lower=0> K;          // Number of predictors
  vector[N] y;             // Response variable
  matrix[N, K] X;          // Predictor matrix
  matrix[N, N] A;          // Network adjacency/weights matrix

  // Prior Parameters
  vector[K] b_mean;        // Beta prior means
  matrix[K, K] b_sigma;    // Beta prior covariance
  real s2_shape;           // Inverse-Gamma shape for s2
  real s2_scale;           // Inverse-Gamma scale for s2
  real rho_mean;           // Mean for rho
  real rho_sd;             // SD for rho
}

parameters{
  vector[K] beta;
  real<lower=0> s2;
  real<lower=-1, upper=1> rho; // Truncated between -1 and 1
}

transformed parameters{
  real sigma = sqrt(s2);
}

model{
  matrix[N, N] I = diag_matrix(rep_vector(1.0, N));
  matrix[N, N] W = I - rho * A;
  vector[N] mu = W * (X * beta); // The mean structure

  // Priors
  beta ~ multi_normal(b_mean, b_sigma);
  s2 ~ inv_gamma(s2_shape, s2_scale);
  rho ~ normal(rho_mean, rho_sd);

  // Likelihood
  // We use the log-determinant because y is transformed by (I - rho*A)
  target += log_determinant(W);
  y ~ normal(mu, sigma);
}
