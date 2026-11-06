#' Homophily Adjusted Network Effects Model (HANE)
#' @param formula A symbolic formula object mapping response to predictors.
#' @param data A data frame containing the variables in the model. Defaults to NULL.
#' @param network A matrix or graph object representing the network structure.
#' @param method Character string specifying estimation ("normal" or "stan").
#' @param beta_prior List containing prior parameters for beta.
#' @param gamma_prior List containing prior parameters for gamma.
#' @param s2_prior List containing prior parameters for s2.
#' @param rho_prior List containing prior parameters for rho.
#' @param U_post A list with elements (Lambda, Omega, Psi) or a 3D array of size (samples, n, D) where size is the number of samples, n is the number of nodes in the network, and D is the number of latent dimensions. If NULL, defaults to a matrix-normal approximation of a spectral embedding.
#' @param rownorm Logical indicating if network matrices should be row-normalized.
#' @param ... Additional arguments.

#' @import mvtnorm
#' @importFrom truncnorm dtruncnorm
#' @importFrom invgamma dinvgamma
#' @importFrom sna lnam
#' @importFrom cmdstanr cmdstan_model write_stan_file
#' @importFrom stats rgamma rnorm
#' @export
hane <- function(formula, data=NULL, network, method,
                 U_post=NULL,
                 D = 2,
                 beta_prior,
                 gamma_prior,
                 s2_prior = list(shape=0.0001, scale=0.0001),
                 rho_prior = list(mean=0.36, sd=0.7),
                 rownorm=TRUE, compile = TRUE, ...){

  # Scavenger for Variables
  vars_in_formula <- all.vars(formula)
  scavenged_list <- list()
  missing_vars <- c()

  pdat <- convert_network(network, "matrix")
  net_df <- pdat$df
  A <- pdat$A

  for(v in vars_in_formula){
    if(!is.null(data) && (v %in% colnames(data))){
      scavenged_list[[v]] <- data[[v]]
    }else if(exists(v, envir = environment(formula))){
      scavenged_list[[v]] <- get(v, envir = environment(formula))
    }else if (!is.null(net_df) && (v %in% colnames(net_df))){
      scavenged_list[[v]] <- net_df[[v]]
    }else{
      missing_vars <- c(missing_vars, v)
    }
  }

  if(length(missing_vars) > 0){
    stop(paste("The following variable(s) required by the formula could not be
               found in 'data', the global environment, or the network node
               attributes:", paste(paste0("'", missing_vars, "'"),
                                   collapse = ", ")))
  }

  data <- as.data.frame(scavenged_list)

  X <- model.matrix(formula, data=data)
  y <- model.frame(formula, data)[,1]
  n <- NROW(y)
  p <- NCOL(X)
  if(rownorm) A <- row_norm(A)

  # U_post Processing
  if(is.null(U_post)){
    Usample <- default_latent_embedding(A, D = D, samples = 100)
    Uapprox <- Uprior_appox(Usample)
    Lambda  <- Uapprox$Lambda
    Omega   <- Uapprox$Omega
    Psi     <- Uapprox$Psi
    D       <- D
  }else if(is.list(U_post)){
    if(!all(c("Lambda", "Omega", "Psi") %in% names(U_post))){
      stop("Custom list inputs for U_post must contain 'Lambda',
           'Omega', and 'Psi'.")
    }
    Lambda <- U_post$Lambda
    Omega  <- U_post$Omega
    Psi    <- U_post$Psi
    D      <- NCOL(Lambda)
  }else if(is.array(U_post) && length(dim(U_post)) == 3){
    Uapprox <- Uprior_appox(U_post)
    Lambda  <- Uapprox$Lambda
    Omega   <- Uapprox$Omega
    Psi     <- Uapprox$Psi
    D       <- dim(U_post)[3]
  }else{
    stop("U_post format unrecognized. Must be NULL, a list of matrices,
         or a 3D sample array.")
  }

  # Build Scale-Invariant Prior for Gamma if argument is omitted
  if(missing(gamma_prior)){
    message("\nThe gamma_prior parameter is not specified. It will be set automatically based on data scales.")
    s_Gamma <- apply(Lambda, 2, sd)
    s_Gamma[s_Gamma == 0] <- 1
    Sigma_gamma_mat <- diag((s_Gamma / 2.5)^2, nrow = D)

    gamma_prior <- list(
      mean = rep(0, D),
      sigma = Sigma_gamma_mat
    )
  }
  # Default Prior for Beta
  if(missing(beta_prior)){
    message("\nThe beta_prior parameter is not specified. It will be set automatically based on data scales.")
    s_j <- apply(X[, -1, drop = FALSE], 2, sd)
    s_j[s_j == 0] <- 1
    s_y <- sd(y)
    if(ncol(X) > 1){
      intercept_var <- max(2.5^2, (2.5 * s_y)^2)
      slope_vars    <- 2.5^2 * s_j^2
      Sigma_beta_mat <- diag(c(intercept_var, slope_vars), nrow = ncol(X))
    }else{
      Sigma_beta_mat <- matrix(2.5^2, nrow = 1, ncol = 1)
    }
    beta_prior <- list(
      mean  = rep(0, ncol(X)),
      sigma = Sigma_beta_mat
    )
  }

  # Validate importance sampling inputs
  is_importance <- tolower(method) %in% c("is", "importance sampling")
  if(is_importance && (typeof(beta_prior) != "list" ||
                        typeof(s2_prior) != "list" ||
                        typeof(rho_prior) != "list" ||
                        typeof(gamma_prior) != "list")){
    stop("importance sampling methods require list input of prior parameters,
         prior functions cannot be changed.")
  }

  # Consolidate Data
  hanedat <- list(y=y, X=X, A=A, p=p, n=n, D=D,
                  Lambda=Lambda, Omega=Omega, Psi=Psi,
                  beta_prior=beta_prior, gamma_prior=gamma_prior,
                  s2_prior=s2_prior, rho_prior=rho_prior,
                  compile = compile, ...)

  # Dispatch
  if(tolower(method) == "normal"){
    haneout <- hane_normal(hanedat)
  }else if(tolower(method)=="stan"){
    haneout <- hane_stan(hanedat)
  }else if(tolower(method) %in% c("is", "importance sampling")){
    haneout <- hane_IS(hanedat)
  }else{
    stop("Unknown estimation method.")
  }

  # Early exit check if user only wants the uncompiled Stan script string
  if(tolower(method) == "stan" && !compile){
    return(haneout)
  }

  # 6. Format Output
  final_out <- list(
    estimates = haneout,
    formula   = formula,
    model     = "hane",
    method    = tolower(method),
    data      = list(y=y, X=X, network=A),
    priors    = list(beta=beta_prior, gamma=gamma_prior, s2=s2_prior,
                     rho=rho_prior,
                     U=list(Lambda=Lambda, Omega=Omega, Psi=Psi)),
    call      = match.call()
  )
  class(final_out) <- c("hane", "netreg")
  return(final_out)
}




# Likelihood of the HANE model.
hane_likelihood <- function(theta, A, y, X, Lambda, Omega, Psi, log=TRUE){
  p <- NCOL(X)
  D <- NCOL(Lambda)
  n <- NROW(y)

  # Extract parameters from theta
  beta <- theta[1:p]
  gamma <- theta[(p+1):(p+D)]
  s2 <- theta[p+D+1]
  rho <- theta[p+D+2]

  I <- diag(n)

  # The 'L' matrix (I - rho*A)^{-1}
  L <- qr.solve(I - rho*A)

  # Calculate the combined variance structure from latent space approximation
  gPg <- as.numeric(t(gamma) %*% Psi %*% gamma)

  # Total internal variance before the network transformation
  Sigma_internal <- gPg * Omega + s2 * I

  # The mean structure: L %*% (X %*% beta + Lambda %*% gamma)
  mu_hane <- L %*% (X %*% beta + Lambda %*% gamma)

  # The total covariance: L %*% Sigma_internal %*% L'
  V_hane <- tcrossprod(L %*% Sigma_internal, L)

  mvtnorm::dmvnorm(y, mean = as.vector(mu_hane), sigma = V_hane, log = log)
}




# Posterior of the Homophily Adjusted Network Effects model
hane_posterior <- function(theta, A, y, X, Lambda, Omega, Psi, log = TRUE,
                           bprior, gprior, s2prior, rhoprior){
  p <- NCOL(X)
  n <- NROW(y)

  # Unpack parameters
  beta  <- theta[1:p]
  gamma <- theta[p + 1]
  s2    <- theta[p + 2]
  rho   <- theta[p + 3]

  # Enforce parameter stationarity and variance bounds
  if(s2 <= 0 || rho <= -1 || rho >= 1){
    if(log){
      return(-Inf)
    }else{
      return(0)
    }
  }

  # Calculate log-likelihood
  log_lik <- hane_likelihood(theta, A, y, X, Lambda, Omega, Psi, log = TRUE)
  if(!is.finite(log_lik)){
    if(log){
      return(-Inf)
    }else{
      return(0)
    }
  }

  # Calculate log-posterior
  log_post <- log_lik + bprior(beta) + gprior(gamma) + s2prior(s2) +
    rhoprior(rho)

  if(log){
    return(log_post)
  }else{
    return(exp(log_post))
  }
}



# Function for approximating the prior for U
Uprior_appox <- function(Usample, maxiter = 100, eps = 1e-8){
  K <- dim(Usample)[1]
  n <- dim(Usample)[2]
  D <- dim(Usample)[3]

  Lambda <- apply(Usample, MARGIN=c(2,3), FUN=mean)

  Psi    <- diag(1, D)
  Omega  <- diag(1, n)

  ll    <- numeric(maxiter)
  ll[1] <- llMatNorm(Usample, Lambda, Omega, Psi)

  for(it in 2:maxiter){
    Psi_inv <- chol2inv(chol(Psi))
    S   <- matrix(0, nrow=n, ncol=n)
    for(j in 1:K){
      S <- S + tcrossprod((Usample[j,,] - Lambda)%*%Psi_inv,
                          Usample[j,,] - Lambda)
    }
    Omega  <- S/S[1,1]
    eta    <- S[1,1]/(K*D)
    Omega[-1,-1] <- eta*Omega[-1,-1] + (1-eta)*tcrossprod(Omega[-1, 1])
    Omega_inv <- chol2inv(chol(Omega))

    Psi <- matrix(0, nrow=D, ncol=D)
    for(j in 1:K){
      Psi <- Psi + crossprod((Usample[j,,] - Lambda),
                             Omega_inv%*%Usample[j,,] - Lambda)/(K*n)
    }

    ll[it] <- llMatNorm(Usample, Lambda, Omega, Psi)
    if( abs((ll[it] - ll[it-1])/ll[it-1]) < eps ) break
  }


  return(list(Lambda = Lambda,
              Omega  = Omega,
              Psi    = Psi,
              conv   = (it <= maxiter)))
}

llMatNorm <- function(X, M, U, V){
  K <- dim(X)[1]
  n <- dim(X)[2]
  D <- dim(X)[3]

  U_inv <- chol2inv(chol(U))
  V_inv <- chol2inv(chol(V))

  0.5*n*K*determinant(V, logarithm = T)$modulus +
    0.5*D*K*determinant(U, logarithm = T)$modulus -
    0.5*sum( apply(X, 1,
                   function(x){sum(diag(tcrossprod(V_inv, x - M)%*%
                                          U_inv%*%(x - M)))}) )
}

default_latent_embedding_final <- function(A, D = 3, samples = 1000) {
  # 1. Initialization and Format Handling ------------------------------------
  if (!inherits(A, "sparseMatrix")) {
    A <- Matrix::Matrix(A, sparse = TRUE)
  }

  n <- nrow(A)
  if (n <= D) stop("n must be strictly greater than D.")

  if (exists("row_norm")) {
    A <- row_norm(row_norm(A), reverse = TRUE)
  }

  is_directed <- !Matrix::isSymmetric(A)

  # 2. Spectral Space Extraction via RSpectra ---------------------------------
  if (is_directed) {
    svd_decomp <- RSpectra::svds(A, k = D)
    r_singular_values <- pmax(svd_decomp$d, 1e-5)
    U_base <- svd_decomp$u %*% diag(sqrt(r_singular_values), nrow = D)
  } else {
    eig_decomp <- RSpectra::eigs(A, k = D, which = "LM")
    r_eigenvalues <- pmax(eig_decomp$values, 1e-5)
    U_base <- eig_decomp$vectors %*% diag(sqrt(r_eigenvalues), nrow = D)
  }

  # 3. Precompute Tensors (Dense Math Operations) -----------------------------
  Delta <- (t(U_base) %*% U_base) / n
  Delta_inv <- solve(Delta)

  S1 <- array(0, dim = c(D, D, D))
  for (j in 1:D) {
    S1[,,j] <- (t(U_base) %*% (U_base[, j] * U_base)) / n
  }

  S2 <- array(0, dim = c(D, D, D, D))
  for (j in 1:D) {
    for (k in 1:D) {
      S2[,,j,k] <- (t(U_base) %*% ((U_base[, j] * U_base[, k]) * U_base)) / n
    }
  }

  # 4. Independent Mixture Mapping & Unique Variance Cache -------------------
  set.seed(123)

  # Independent sample paths: rows = samples, columns = nodes
  h_matrix <- matrix(sample(1:n, size = samples * n, replace = TRUE),
                     nrow = samples, ncol = n)

  # Flatten h_matrix column-by-column to perfectly match a 2D noise design
  h_vector <- as.vector(h_matrix)
  unique_active_h <- unique(h_vector)

  # Calculate unique covariance matrices and cache their Cholesky factors
  chol_list <- list()
  for (h in unique_active_h) {
    x_hat_h <- U_base[h, ]

    inner_M1 <- matrix(0, nrow = D, ncol = D)
    for (j in 1:D) {
      inner_M1 <- inner_M1 + x_hat_h[j] * S1[,,j]
    }

    inner_M2 <- matrix(0, nrow = D, ncol = D)
    for (j in 1:D) {
      for (k in 1:D) {
        inner_M2 <- inner_M2 + (x_hat_h[j] * x_hat_h[k]) * S2[,,j,k]
      }
    }

    Sigma_h <- Delta_inv %*% (inner_M1 - inner_M2) %*% Delta_inv
    V <- Sigma_h / n

    min_ev <- min(eigen(V, only.values = TRUE)$values)
    if (min_ev < 1e-10) {
      diag(V) <- diag(V) + (abs(min_ev) + 1e-10)
    }
    V <- 0.5 * (V + t(V))

    chol_list[[as.character(h)]] <- chol(V)
  }

  # 5. Flat 2D Noise Transformation & Reshaping -------------------------------
  N_total <- samples * n
  raw_noise <- matrix(rnorm(N_total * D), nrow = N_total, ncol = D)
  transformed_noise <- matrix(0, nrow = N_total, ncol = D)

  # Map and transform noise using a fast block assignment on the 2D layout
  for (h in unique_active_h) {
    matched_indices <- which(h_vector == h)
    if (length(matched_indices) > 0) {
      transformed_noise[matched_indices, ] <- raw_noise[matched_indices, , drop = FALSE] %*% chol_list[[as.character(h)]]
    }
  }

  # Reshape flat 2D noise matrix into a structured 3D array (samples, n, D)
  # R fills arrays column-first, which perfectly matches our column-flattened h_vector
  transformed_noise_array <- array(transformed_noise, dim = c(samples, n, D))

  # Broadcast the transformed noise array over the initial baseline coordinate means
  Usample <- array(0, dim = c(samples, n, D))
  for (s in 1:samples) {
    Usample[s, , ] <- U_base + transformed_noise_array[s, , ]
  }

  dimnames(Usample) <- list(
    sample = paste0("S", 1:samples),
    node   = paste0("N", 1:n),
    dim    = paste0("D", 1:D)
  )

  return(Usample)
}







# Code for when method="normal"
hane_normal <- function(hanedat){
  p <- hanedat$p

  # Beta Prior
  if(typeof(hanedat$beta_prior) == "list"){
    if(length(hanedat$beta_prior$mean) == 1){
      bmean <- rep(hanedat$beta_prior$mean, p)
    }else{
      bmean <- hanedat$beta_prior$mean
    }
    if(length(hanedat$beta_prior$sigma) == 1){
      bsigma <- diag(p) * hanedat$beta_prior$sigma
    }else{
      bsigma <- hanedat$beta_prior$sigma
    }
    bprior <- function(beta){
      mvtnorm::dmvnorm(beta, bmean, bsigma, log = TRUE)
    }
  }else{
    bprior <- hanedat$beta_prior
  }

  # Gamma Prior
  if(typeof(hanedat$gamma_prior) == "list"){
    gprior <- function(gamma) dnorm(gamma, hanedat$gamma_prior$mean,
                                    sqrt(hanedat$gamma_prior$sigma), log = TRUE)
  }else{
    gprior <- hanedat$gamma_prior
  }

  # s2 Prior
  if(typeof(hanedat$s2_prior) == "list"){
    s2prior <- function(s2){
      invgamma::dinvgamma(s2, hanedat$s2_prior$shape,
                          hanedat$s2_prior$scale, log = TRUE)
    }
  }else{
    s2prior <- hanedat$s2_prior
  }

  # Rho Prior
  if(typeof(hanedat$rho_prior) == "list"){
    rhoprior <- function(rho){
      log(truncnorm::dtruncnorm(rho, a = -1, b = 1,
                                mean = hanedat$rho_prior$mean,
                                sd = hanedat$rho_prior$sd))
    }
  }else{
    rhoprior <- hanedat$rho_prior
  }


  # Initial coordinates
  init_betas <- as.numeric(lm.fit(hanedat$X, hanedat$y)$coefficients)
  init_par   <- c(init_betas, 0, 1, 0.36) # Beta, Gamma=0, s2=1, rho=0.36

  # Checking if additional Arguments were input
  if(!is.null(hanedat$gr)){
    gr=hanedat$gr
  }else{
    gr=NULL
  }

  if(!is.null(hanedat$lower)){
    lower <- hanedat$lower
  }else{
    lower <- c(rep(-Inf, p), -Inf, 1e-6, -0.999) # Goes c(Beta, Gamma, s2, rho)
  }

  if(!is.null(hanedat$upper)){
    upper <- hanedat$upper
  }else{
    upper <- c(rep(Inf, p), Inf, Inf, 0.999)
  }

  if(!is.null(hanedat$control)){
    control=hanedat$control
  }else{
    control=list(fnscale = -1)
  }

  # Execute optimization using L-BFGS-B
  myopt <- optim(par = init_par, fn = hane_posterior, gr = gr,
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = control,
    hessian = TRUE,

    A       = hanedat$A,
    y       = hanedat$y,
    X       = hanedat$X,
    Lambda  = hanedat$Lambda,
    Omega   = hanedat$Omega,
    Psi     = hanedat$Psi,
    log     = TRUE,
    bprior  = bprior,
    gprior  = gprior,
    s2prior = s2prior,
    rhoprior = rhoprior
  )

  beta_ests <- myopt$par[1:p]
  names(beta_ests) <- colnames(hanedat$X)

  list(
    beta      = beta_ests,
    gamma     = myopt$par[p + 1],
    s2        = myopt$par[p + 2],
    rho       = myopt$par[p + 3],
    hessian   = myopt$hessian,
    loglik    = hane_likelihood(myopt$par, hanedat$A, hanedat$y,
                                hanedat$X, hanedat$Lambda, hanedat$Omega,
                                hanedat$Psi, log = TRUE),
    logpost   = myopt$value,
    counts    = myopt$counts,
    message   = myopt$message,
    converged = myopt$convergence == 0
  )
}



# Code for when method="stan"
hane_stan <- function(hanedat){
  if(!requireNamespace("cmdstanr", quietly = TRUE)){
    stop("Package 'cmdstanr' is required for this method. Please install it.",
         call. = FALSE)
  }

  if(is.null(hanedat$iterations)){
    iterations <- 2000
  }else{
    iterations <- hanedat$iterations
  }
  if(is.null(hanedat$chains)){
    chains <- 4
  }else{
    chains <- hanedat$chains
  }
  K <- hanedat$p

  # Format beta prior line
  if(is.list(hanedat$beta_prior)){
    if(length(hanedat$beta_prior$mean) == 1){
      b_mu <- rep(hanedat$beta_prior$mean, K)
    }else{
      b_mu <- as.vector(hanedat$beta_prior$mean)
    }
    if(length(hanedat$beta_prior$sigma) == 1){
      b_sig <- diag(K) * hanedat$beta_prior$sigma
    }else{
      b_sig <- as.matrix(hanedat$beta_prior$sigma)
    }

    b_mean_str <- paste0("[", paste(b_mu, collapse = ", "), "]'")
    b_sigma_rows <- apply(b_sig, 1,
                          function(row) paste0("[",
                                               paste(row, collapse = ", "),
                                               "]"))
    b_sigma_str <- paste0("[", paste(b_sigma_rows, collapse = ", "), "]")

    beta_prior_line <- paste0("beta ~ multi_normal(", b_mean_str, ", ",
                              b_sigma_str, ");")
  }else if(is.character(hanedat$beta_prior)){
    beta_prior_line <- hanedat$beta_prior
  }else{
    stop("beta_prior must be a list of parameters
         or a valid Stan character string.")
  }

  # Format gamma prior line
  if(is.list(hanedat$gamma_prior)){
    if(length(hanedat$gamma_prior$mean) == 1){
      g_mu <- rep(hanedat$gamma_prior$mean, hanedat$D)
    }else{
      g_mu <- as.vector(hanedat$gamma_prior$mean)
    }
    if(length(hanedat$gamma_prior$sigma) == 1){
      g_sig <- diag(hanedat$D) * hanedat$gamma_prior$sigma
    }else{
      g_sig <- as.matrix(hanedat$gamma_prior$sigma)
    }

    g_mean_str <- paste0("[", paste(g_mu, collapse = ", "), "]'")
    g_sigma_rows <- apply(g_sig, 1,
                          function(row) paste0("[",
                                               paste(row, collapse = ", "),
                                               "]"))
    g_sigma_str <- paste0("[", paste(g_sigma_rows, collapse = ", "), "]")

    gamma_prior_line <- paste0("gamma ~ multi_normal(", g_mean_str, ", ",
                               g_sigma_str, ");")
  }else if(is.character(hanedat$gamma_prior)){
    gamma_prior_line <- hanedat$gamma_prior
  }else{
    stop("gamma_prior must be a list of parameters or a valid Stan character string.")
  }

  # Format s2 prior line
  if(is.list(hanedat$s2_prior)){
    s2_prior_line <- paste0("s2 ~ inv_gamma(", hanedat$s2_prior$shape, ", ",
                            hanedat$s2_prior$scale, ");")
  }else if(is.character(hanedat$s2_prior)){
    s2_prior_line <- hanedat$s2_prior
  }else{
    stop("s2_prior must be a list of parameters
         or a valid Stan character string.")
  }

  # Format rho prior line
  if(is.list(hanedat$rho_prior)){
    rho_prior_line <- paste0("rho ~ normal(", hanedat$rho_prior$mean, ", ",
                             hanedat$rho_prior$sd, ");")
  }else if(is.character(hanedat$rho_prior)){
    rho_prior_line <- hanedat$rho_prior
  }else{
    stop("rho_prior must be a list of parameters
         or a valid Stan character string.")
  }

  # Uncollapsed Stan Code
  stan_template <- paste0("
data {
  int<lower=0> N;                // Number of nodes
  int<lower=0> K;                // Number of predictors
  int<lower=0> D;                // Latent position dimensions
  vector[N] y;                   // Response variable
  matrix[N, K] X;                // Predictor matrix
  matrix[N, N] A;                // Adjacency matrix

  // Matrix-Normal Prior Parameters for Homophily Correction
  matrix[N, D] Lambda;           // Mean matrix positions
  matrix[N, N] Omega;            // Row-covariance matrix
  matrix[D, D] Psi;              // Column-covariance matrix
}

transformed data {
  // Pre-compute Cholesky factors for the Kronecker structure once
  matrix[N, N] L_Omega = cholesky_decompose(Omega);
  matrix[D, D] L_Psi = cholesky_decompose(Psi);
}

parameters {
  vector[K] beta;
  vector[D] gamma;               // Mappings for the D latent features
  real<lower=0> s2;
  real<lower=-1, upper=1> rho;
  matrix[N, D] Z;                // Standardized normal latent positions
}

transformed parameters {
  real sigma = sqrt(s2);

  // Efficient Matrix-Normal mapping: U = Lambda + L_Omega * Z * L_Psi'
  matrix[N, D] U = Lambda + L_Omega * Z * L_Psi';
}

model {
  matrix[N, N] I = diag_matrix(rep_vector(1.0, N));
  matrix[N, N] W = I - rho * A;

  // Base regression transformation block matching Equation 1:
  // (I - rho*A)y = X*beta + U*gamma + epsilon
  vector[N] y_transformed = W * y;
  vector[N] mu = X * beta + U * gamma;

  // Priors dynamically inserted via R
  ", beta_prior_line, "
  ", gamma_prior_line, "
  ", s2_prior_line, "
  ", rho_prior_line, "

  // Sampling the standardized matrix elements implies U ~ MN(Lambda, Omega, Psi)
  to_vector(Z) ~ std_normal();

  // Likelihood with spatial Jacobian correction
  target += log_determinant(W);
  y_transformed ~ normal(mu, sigma);
}

")

  if(!is.null(hanedat$compile) && !hanedat$compile){
    return(stan_template)
  }

  stan_file_tmp <- cmdstanr::write_stan_file(stan_template)
  compiled_model <- cmdstanr::cmdstan_model(stan_file_tmp)

  stan_data <- list(
    N      = hanedat$n,
    K      = K,
    D      = hanedat$D,
    y      = hanedat$y,
    X      = hanedat$X,
    A      = hanedat$A,
    Lambda = hanedat$Lambda,
    Omega  = hanedat$Omega,
    Psi    = hanedat$Psi
  )

  # Extract stan arguments
  all_args <- names(formals(cmdstanr:::CmdStanModel$public_methods$sample))
  valid_sample_args <- setdiff(all_args, "data")
  user_sampling_params <- hanedat[names(hanedat) %in% valid_sample_args]

  # Sample
  final_args <- c(
    list(
      data = stan_data,
      iter_sampling = iterations,
      iter_warmup = iterations,
      chains = chains,
      refresh = 0
    ),
    user_sampling_params            # Overrides or adds choices from user
  )
  final_args <- final_args[!duplicated(names(final_args), fromLast = TRUE)]
  cmdstan_fit <- do.call(compiled_model$sample, final_args)

  draws_matrix <- cmdstan_fit$draws(format = "draws_matrix")
  attr(draws_matrix, "stan_code") <- stan_template

  return(draws_matrix)
}





lnam2sls_effect = function(y,X,A){
  ZZ = cbind(X,A%*%y)
  if(mean(X[,1] == 1) == 1){
    X1 = X[,-1]
  }else{
    X1 = X
  }
  AX = A%*%X1
  HH = model.matrix(~X1+AX+A%*%AX)
  PP = tcrossprod(HH%*%chol2inv(chol(crossprod(HH))),HH)
  ZHat = PP%*%ZZ
  ZtZInv= qr.solve(crossprod(ZHat))
  Coefs = drop(ZtZInv%*%crossprod(ZHat,y))

  s2Hat = drop(crossprod(y-ZZ%*%Coefs))/nrow(X)

  return(list(coefs=Coefs,s2=s2Hat))
}



#' Summary method for homophily adjusted network effects objects
#' @method summary hane
#' @export
summary.hane <- function(object, probs = c(0.025, 0.975), ...){
  p <- ncol(object$data$X)
  prob_names <- paste0(probs * 100, "%")

  if(object$method == "stan"){
    draws <- object$estimates

    # FIX: Use grep to capture all dimensions of gamma (e.g., gamma[1], gamma[2])
    beta_cols  <- grep("^beta\\[", colnames(draws))
    gamma_cols <- grep("^gamma\\[", colnames(draws))
    s2_col     <- which(colnames(draws) == "s2")
    rho_col    <- which(colnames(draws) == "rho")
    lp_col     <- which(colnames(draws) == "lp__")

    target_cols <- c(beta_cols, gamma_cols, s2_col, rho_col)

    means <- colMeans(draws[, target_cols, drop = FALSE])
    quants <- apply(draws[, target_cols, drop = FALSE], 2, quantile,
                    probs = probs, na.rm = TRUE)

    if(length(target_cols) == 1){
      coef_table <- matrix(c(means, quants), nrow = 1)
    }else{
      coef_table <- cbind(means, t(quants))
    }

    colnames(coef_table) <- c("Estimate", prob_names)

    # FIX: Dynamically generate row names matching the number of gammas
    gamma_names <- paste0("gamma[", 1:length(gamma_cols), "]")
    rownames(coef_table) <- c(colnames(object$data$X), gamma_names, "sigma2", "rho")

    logpost_val <- mean(draws[, lp_col])
    theta_star  <- as.numeric(means)
    loglik_val  <- hane_likelihood(theta_star, object$data$network,
                                   object$data$y, object$data$X,
                                   object$priors$U$Lambda,
                                   object$priors$U$Omega, object$priors$U$Psi,
                                   log=TRUE)

  }else{
    # Normal Approximation
    # FIX: Dynamically handle length of gamma in the normal estimation block too
    D_gamma <- length(object$estimates$gamma)
    theta_star <- c(object$estimates$beta, object$estimates$gamma,
                    object$estimates$s2, object$estimates$rho)
    hinv <- tryCatch(solve(-object$estimates$hessian),
                     error = function(e) matrix(NA, length(theta_star),
                                                length(theta_star)))
    std_errs <- sqrt(diag(hinv))

    z_crit <- qnorm(probs)
    intervals <- matrix(NA, nrow = length(theta_star), ncol = length(probs))
    for(i in 1:length(probs)) {
      intervals[,i] <- theta_star + z_crit[i] * std_errs
    }

    coef_table <- cbind(Estimate = theta_star, intervals)
    colnames(coef_table) <- c("Estimate", prob_names)

    gamma_names <- paste0("gamma[", 1:D_gamma, "]")
    rownames(coef_table) <- c(colnames(object$data$X), gamma_names, "sigma2", "rho")

    loglik_val  <- object$estimates$loglik
    logpost_val <- object$estimates$logpost
  }

  res <- list(
    call         = object$call,
    coefficients = coef_table,
    loglik       = loglik_val,
    logpost      = logpost_val,
    marg_lik     = NA,
    model        = "hane",
    method       = object$method,
    converged = if(object$method == "normal") object$estimates$converged else T
  )

  class(res) <- "summary.netreg"
  return(res)
}



#' Print method for hane objects
#' @method print hane
#' @export
print.hane <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  sm <- summary(x)
  coefs <- sm$coefficients[, "Estimate"]

  cat("Coefficients & Parameters:\n")
  print(round(coefs, 5))
  cat("\n")

  invisible(x)
}


# Internal utility to calculate fitted values and residuals for HANE.
# Keeps NetReg::: completely clean of S3 clutter.
calc_hane_diags <- function(object) {
  sm    <- summary(object)
  coefs <- sm$coefficients

  # 1. Extract parameter values cleanly using unified names
  beta  <- coefs[colnames(object$data$X), "Estimate"]
  rho   <- coefs["rho", "Estimate"]

  # Grab gamma dimensions dynamically from row names
  gamma_rows <- grep("^gamma\\[", rownames(coefs), value = TRUE)
  gamma      <- coefs[gamma_rows, "Estimate"]

  X      <- object$data$X
  y      <- object$data$y
  W      <- object$data$network
  Lambda <- object$priors$U$Lambda

  # 2. Compute diagnostics accounting for latent homophily correction
  fitted_vals <- as.vector(X %*% beta + Lambda %*% gamma + rho * (W %*% y))
  resids      <- as.vector(y - fitted_vals)

  return(list(fitted = fitted_vals, residuals = resids))
}


#' Plot method for hane objects
#' @method plot hane
#' @export
plot.hane <- function(x, ...) {
  diags    <- calc_hane_diags(x)
  fit_vals <- diags$fitted
  res_vals <- diags$residuals

  # Minimal diagnostic scatter plot inheriting base user margins
  plot(fit_vals, res_vals,
       xlab = "Fitted values",
       ylab = "Residuals",
       main = "Residuals vs Fitted",
       panel.first = grid(),
       pch = 20, col = "darkgray", ...)

  # Horizontal reference line at 0
  abline(h = 0, lty = 3, col = "gray40", lwd = 1.5)

  # Solid black trend line
  if (length(res_vals) > 3) {
    lines(lowess(fit_vals, res_vals), col = "black", lwd = 1.5)
  }

  invisible(x)
}














