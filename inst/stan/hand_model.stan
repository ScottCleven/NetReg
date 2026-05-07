data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> D;
  vector[N] y;
  matrix[N, K] X;
  matrix[N, N] A;
  matrix[N, D] Lambda;
  matrix[N, N] Omega;
  matrix[D, D] Psi;

  // Priors
  vector[K] b_mean;
  matrix[K, K] b_sigma;
  vector[D] g_mean;
  matrix[D, D] g_sigma;
  real s2_shape;
  real s2_scale;
  real rho_mean;
  real rho_sd;
}

parameters {
  vector[K] beta;
  vector[D] gamma;
  real<lower=0> s2;
  real<lower=-1, upper=1> rho;
}

model {
  matrix[N, N] W = diag_matrix(rep_vector(1.0, N)) - rho * A;
  matrix[N, N] M = inverse(W);

  real gPg = dot_product(gamma, Psi * gamma);
  matrix[N, N] Sigma = gPg * Omega + M * (s2 * diag_matrix(rep_vector(1.0, N))) * M';

  vector[N] mu = X * beta + Lambda * gamma;


  y ~ multi_normal(mu, Sigma);
}
