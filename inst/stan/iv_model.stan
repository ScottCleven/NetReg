data {
  int<lower=0> N;
  int<lower=0> K;
  vector[N] y;
  matrix[N, K] X;
  matrix[N, N] A;

  // Priors
  vector[K] b_mean;
  matrix[K, K] b_sigma;
  real s2_shape;
  real s2_scale;
  real rho_mean;
  real rho_sd;
  real lambda_mean;
  real lambda_sd;
}

parameters {
  vector[K] beta;
  real<lower=0> s2;
  real<lower=-1, upper=1> rho;
  real<lower=-1, upper=1> lambda;
}

transformed parameters {
  real sigma = sqrt(s2);
}

model {
  matrix[N, N] I = diag_matrix(rep_vector(1.0, N));
  matrix[N, N] W_rho = I - rho * A;
  matrix[N, N] W_lambda = I - lambda * A;

  // The 'Innovation' vector (epsilon)
  // epsilon = (I - lambda*A) * ( (I - rho*A)*y - X*beta )
  vector[N] epsilon = W_lambda * (W_rho * y - X * beta);

  // Priors
  beta ~ multi_normal(b_mean, b_sigma);
  s2 ~ inv_gamma(s2_shape, s2_scale);
  rho ~ normal(rho_mean, rho_sd);
  lambda ~ normal(lambda_mean, lambda_sd);

  // Likelihood with dual Jacobian adjustments
  target += log_determinant(W_rho);
  target += log_determinant(W_lambda);

  // Since epsilon is the isolated N(0, sigma) noise:
  epsilon ~ normal(0, sigma);
}










