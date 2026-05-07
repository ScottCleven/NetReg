data {
  int<lower=0> N;
  int<lower=0> K;
  vector[N] y;
  matrix[N, K] X;
  matrix[N, N] A;
  // Prior Parameters
  vector[K] b_mean;
  matrix[K, K] b_sigma;
  real s2_shape;
  real s2_scale;
  real rho_mean;
  real rho_sd;
}

parameters {
  vector[K] beta;
  real<lower=0> s2;
  real<lower=-1, upper=1> rho;
}

transformed parameters {
  real sigma = sqrt(s2);
}

model {
  matrix[N, N] I = diag_matrix(rep_vector(1.0, N));
  matrix[N, N] W = I - rho * A;

  vector[N] epsilon = W * (y - X * beta);

  // Priors
  beta ~ multi_normal(b_mean, b_sigma);
  s2 ~ inv_gamma(s2_shape, s2_scale);
  rho ~ normal(rho_mean, rho_sd);

  // Likelihood
  target += log_determinant(W);
  epsilon ~ normal(0, sigma);
}
