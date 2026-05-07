data {
  int<lower=1> N;
  vector[N] Y;
  int<lower=1> K;
  matrix[N, K] X;
  int<lower=1> N_1; // Number of communities
  int<lower=1> M_1; // Number of moderated variables
  array[N] int<lower=1> J_1;
  matrix[N, M_1] Z_1;

  // Prior parameters
  vector[K] b_mean;
  matrix[K, K] b_sigma;
  real psi_mean;
  real psi_sd;
  real lkj_shape;
}
parameters {
  vector[K] b;
  real<lower=0> sigma;
  vector[M_1] psi;
  matrix[M_1, N_1] z_1;
  cholesky_factor_corr[M_1] L_1;
}
transformed parameters {
  vector[M_1] lambda = exp(psi);
  matrix[N_1, M_1] r_1;
  r_1 = transpose(diag_pre_multiply(lambda * sigma, L_1) * z_1);
}
model {
  vector[N] mu = X * b;
  for (n in 1:N) {
    mu[n] += dot_product(r_1[J_1[n]], Z_1[n]);
  }
  // Priors using user-defined parameters
  b ~ multi_normal(b_mean, b_sigma);
  sigma ~ student_t(3, 0, 2.5); // Global scale
  psi ~ normal(psi_mean, psi_sd);
  L_1 ~ lkj_corr_cholesky(lkj_shape);
  to_vector(z_1) ~ std_normal();

  Y ~ normal(mu, sigma);
}
