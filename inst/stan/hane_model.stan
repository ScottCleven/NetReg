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
  matrix[N, N] I = diag_matrix(rep_vector(1.0, N));
  matrix[N, N] W = I - rho * A;

  // Structured Covariance: (gamma' * Psi * gamma) * Omega + s2 * I
  real gPg = dot_product(gamma, Psi * gamma);
  matrix[N, N] Sigma = gPg * Omega + s2 * I;

  // Mean: X*beta + Lambda*gamma
  // Note: We model W*y ~ Normal(mean, Sigma) which implies Jacobian adjustment
  vector[N] mu = X * beta + Lambda * gamma;

  // Priors
  beta ~ multi_normal(b_mean, b_sigma);
  gamma ~ multi_normal(g_mean, g_sigma);
  s2 ~ inv_gamma(s2_shape, s2_scale);
  rho ~ normal(rho_mean, rho_sd);

  // Likelihood with Jacobian adjustment
  target += log_determinant(W);
  W * y ~ multi_normal(mu, Sigma);
}
