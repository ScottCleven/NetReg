#' Homophily Adjusted Network Disturbances Model (HAND)
#' @param formula A symbolic formula object mapping response to predictors.
#' @param data A data frame containing the variables in the model. Defaults to NULL.
#' @param network A matrix or graph object representing the network structure.
#' @param U_post A list with elements (Lambda, Omega, Psi) or a 3D array of size (samples, n, D). If NULL, defaults to a matrix-normal approximation of a spectral embedding.
#' @param method Character string specifying estimation ("normal" or "stan").
#' @param beta_prior List containing prior parameters for beta.
#' @param gamma_prior List containing prior parameters for gamma.
#' @param s2_prior List containing prior parameters for s2.
#' @param rho_prior List containing prior parameters for rho.
#' @param rownorm Logical indicating if network matrices should be row-normalized.
#' @param ... Additional arguments.
#' @import mvtnorm
#' @importFrom truncnorm dtruncnorm
#' @importFrom invgamma dinvgamma
#' @importFrom sna lnam
#' @importFrom cmdstanr cmdstan_model write_stan_file
#' @importFrom stats rgamma rnorm
#' @export
hand <- function(formula, data=NULL, network, method,
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

  # Look for each variable one by one in order of data>environment>network
  for(v in vars_in_formula){
    if(!is.null(data) && (v %in% colnames(data))){
      scavenged_list[[v]] <- data[[v]]
    }else if(exists(v, envir = environment(formula))){
      scavenged_list[[v]] <- get(v, envir = environment(formula))
    }else if (!is.null(net_df) && (v %in% colnames(net_df))){
      scavenged_list[[v]] <- net_df[[v]]
    }else{
      # Variable is missing everywhere
      missing_vars <- c(missing_vars, v)
    }
  }

  # Missing variables error
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

  # U_post
  if(is.null(U_post)){
    Usample <- default_latent_embedding(A, D = D, samples = 100)
    Uapprox <- Uprior_appox(Usample)
    Lambda  <- Uapprox$Lambda
    Omega   <- Uapprox$Omega
    Psi     <- Uapprox$Psi
    D       <- D
  }else if(is.list(U_post)){
    if(!all(c("Lambda", "Omega", "Psi") %in% names(U_post))){
      stop("Custom list inputs for U_post must contain
           'Lambda', 'Omega', and 'Psi'.")
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

  # importance sampling prior inputs
  is_importance <- tolower(method) %in% c("is", "importance sampling")
  if (is_importance && (typeof(beta_prior) != "list" ||
                        typeof(s2_prior) != "list" ||
                        typeof(rho_prior) != "list" ||
                        typeof(gamma_prior) != "list")) {
    stop("importance sampling methods require list input of prior parameters,
         prior functions cannot be changed.")
  }

  # Consolidate Data
  handdat <- list(y=y, X=X, A=A, p=p, n=n, D=D,
                  Lambda=Lambda, Omega=Omega, Psi=Psi,
                  beta_prior=beta_prior, gamma_prior=gamma_prior,
                  s2_prior=s2_prior, rho_prior=rho_prior,
                  compile = compile, ...)


  if(tolower(method) == "normal"){
    handout <- hand_normal(handdat)
  }else if(tolower(method)=="stan"){
    handout <- hand_stan(handdat)
  }else{
    stop("Unknown estimation method. HAND supports 'normal' or 'stan'.")
  }

  # Early exit check if user only wanted the uncompiled Stan script string
  if(tolower(method) == "stan" && !compile){
    return(handout)
  }

  # Format Output
  final_out <- list(
    estimates = handout,
    formula   = formula,
    model     = "hand",
    method    = tolower(method),
    data      = list(y=y, X=X, network=A),
    priors    = list(beta=beta_prior, gamma=gamma_prior, s2=s2_prior,
                     rho=rho_prior,
                     U=list(Lambda=Lambda, Omega=Omega, Psi=Psi)),
    call      = match.call()
  )
  class(final_out) <- c("hand", "netreg")
  return(final_out)
}



# Log-Likelihood for HAND
hand_likelihood <- function(theta, A, y, X, Lambda, Omega, Psi, log = TRUE){
  p <- NCOL(X)
  n <- NROW(y)
  D <- NCOL(Lambda) # Detect the number of latent space dimensions

  beta  <- theta[1:p]
  gamma <- theta[(p + 1):(p + D)]
  s2    <- theta[p + D + 1]
  rho   <- theta[p + D + 2]

  # Guard clauses for parameter boundaries
  if(s2 <= 0 || rho <= -1 || rho >= 1){
    return(if(log) -Inf else 0)
  }

  # log|I - rho*A| = sum(log(1 - rho * eigenvalues))
  ev <- eigen(A, only.values = TRUE)$values
  log_det_W <- sum(log(1 - rho * Re(ev)))

  if(is.na(log_det_W) || is.nan(log_det_W)){
    return(if(log) -Inf else 0)
  }

  # nu = y - X*beta - Lambda*gamma
  nu <- as.numeric(y - X %*% beta - Lambda %*% gamma)

  # Apply spatial transformation matrix components: (I - rho*A) %*% nu
  epsilon <- as.numeric(nu - rho * as.numeric(A %*% nu))

  # Calculate log-likelihood
  log_lik <- sum(dnorm(epsilon, mean = 0, sd = sqrt(s2), log = TRUE)) + log_det_W

  if(log){
    return(log_lik)
  }else{
    return(exp(log_lik))
  }
}


# Log-Posterior for HAND
hand_posterior <- function(theta, A, y, X, Lambda, Omega, Psi, log = TRUE,
                           bprior, gprior, s2prior, rhoprior){
  p <- NCOL(X)
  s2    <- theta[p + 2]
  rho   <- theta[p + 3]

  if(s2 <= 0 || rho <= -1 || rho >= 1){
    if(log){
      return(-Inf)
    }else{
      return(0)
    }
  }

  log_lik <- hand_likelihood(theta, A, y, X, Lambda, Omega, Psi, log = TRUE)
  if(!is.finite(log_lik)){
    if(log){
      return(-Inf)
    }else{
      return(0)
    }
  }
  log_post <- log_lik + bprior(theta[1:p]) + gprior(theta[p + 1]) +
    s2prior(s2) + rhoprior(rho)
  if(log){
    return(log_post)
  }else{
    return(exp(log_post))
  }
}


# Code for when method="normal"
hand_normal <- function(handdat){
  p <- handdat$p

  if(typeof(handdat$beta_prior) == "list"){
    if(length(handdat$beta_prior$mean) == 1){
      bmean <- rep(handdat$beta_prior$mean, p)
    }else{
      bmean <- handdat$beta_prior$mean
    }
    if(length(handdat$beta_prior$sigma) == 1){
      bsigma <- diag(p) * handdat$beta_prior$sigma
    }else{
      bsigma <- handdat$beta_prior$sigma
    }
    bprior <- function(beta){
      mvtnorm::dmvnorm(beta, bmean, bsigma, log = TRUE)
    }
  }else{
    bprior <- handdat$beta_prior
  }

  if(typeof(handdat$gamma_prior) == "list"){
    gprior <- function(gamma){
      dnorm(gamma, handdat$gamma_prior$mean, sqrt(handdat$gamma_prior$sigma),
            log = TRUE)
    }
  }else{
    gprior <- handdat$gamma_prior
  }

  if(typeof(handdat$s2_prior) == "list"){
    s2prior <- function(s2) invgamma::dinvgamma(s2, handdat$s2_prior$shape,
                                                handdat$s2_prior$scale,
                                                log = TRUE)
  }else{
    s2prior <- handdat$s2_prior
  }

  if(typeof(handdat$rho_prior) == "list"){
    rhoprior <- function(rho){
      log(truncnorm::dtruncnorm(rho, a = -1, b = 1,
                                mean = handdat$rho_prior$mean,
                                sd = handdat$rho_prior$sd))
    }
  }else{
    rhoprior <- handdat$rho_prior
  }


  init_betas <- as.numeric(lm.fit(handdat$X, handdat$y)$coefficients)
  init_par   <- c(init_betas, 0, 1, 0.36)

  # Checking if additional Arguments were input
  if(!is.null(handdat$gr)){
    gr=handdat$gr
  }else{
    gr=NULL
  }

  if(!is.null(handdat$lower)){
    lower <- handdat$lower
  }else{
    lower <- c(rep(-Inf, p), -Inf, 1e-6, -0.999) # Goes c(Beta, Gamma, s2, rho)
  }

  if(!is.null(handdat$upper)){
    upper <- handdat$upper
  }else{
    upper <- c(rep(Inf, p), Inf, Inf, 0.999)
  }

  if(!is.null(handdat$control)){
    control=handdat$control
  }else{
    control=list(fnscale = -1)
  }

  myopt <- optim(par = init_par, fn = hand_posterior, gr = gr,
    method   = "L-BFGS-B",
    lower    = lower,
    upper    = upper,
    control  = list(fnscale = -1),
    hessian  = TRUE,

    A        = handdat$A,
    y        = handdat$y,
    X        = handdat$X,
    Lambda   = handdat$Lambda,
    Omega    = handdat$Omega,
    Psi      = handdat$Psi,
    log      = TRUE,
    bprior   = bprior,
    gprior   = gprior,
    s2prior  = s2prior,
    rhoprior = rhoprior
  )

  beta_ests <- myopt$par[1:p]
  names(beta_ests) <- colnames(handdat$X)

  list(
    beta = beta_ests, gamma = myopt$par[p + 1], s2 = myopt$par[p + 2],
    rho = myopt$par[p + 3], hessian = myopt$hessian,
    loglik = hand_likelihood(myopt$par, handdat$A, handdat$y, handdat$X,
                             handdat$Lambda, handdat$Omega, handdat$Psi,
                             log = TRUE),
    logpost = myopt$value, converged = myopt$convergence == 0
  )
}


# Code for when method="stan"
hand_stan <- function(handdat){
  if(!requireNamespace("cmdstanr", quietly = TRUE)){
    stop("Package 'cmdstanr' is needed for this method. Please install it.",
         call. = FALSE)
  }
  if(is.null(handdat$iterations)){
    iterations <- 2000
  }else{
    iterations <- handdat$iterations
  }
  if(is.null(handdat$chains)){
    chains <- 4
  }else{
    chains <- handdat$chains
  }
  K <- handdat$p

  # beta prior
  if(is.list(handdat$beta_prior)){
    if(length(handdat$beta_prior$mean) == 1){
      b_mu <- rep(handdat$beta_prior$mean, K)
    }else{
      b_mu <- as.vector(handdat$beta_prior$mean)
    }
    if(length(handdat$beta_prior$sigma) == 1){
      b_sig <- diag(K) * handdat$beta_prior$sigma
    }else{
      b_sig <- as.matrix(handdat$beta_prior$sigma)
    }

    b_mean_str <- paste0("[", paste(b_mu, collapse = ", "), "]'")
    b_sigma_rows <- apply(b_sig, 1,
                          function(row) paste0("[",
                                               paste(row, collapse = ", "),
                                               "]"))
    b_sigma_str <- paste0("[", paste(b_sigma_rows, collapse = ", "), "]")

    beta_prior_line <- paste0("beta ~ multi_normal(", b_mean_str, ", ",
                              b_sigma_str, ");")
  }else if(is.character(handdat$beta_prior)){
    beta_prior_line <- handdat$beta_prior
  }else{
    stop("beta_prior must be a list of parameters
         or a valid Stan character string.")
  }

  # gamma prior
  if(is.list(handdat$gamma_prior)){
    if(length(handdat$gamma_prior$mean) == 1){
      g_mu <- rep(handdat$gamma_prior$mean, handdat$D)
    }else{
      g_mu <- as.vector(handdat$gamma_prior$mean)
    }
    if(length(handdat$gamma_prior$sigma) == 1){
      g_sig <- diag(handdat$D) * handdat$gamma_prior$sigma
    }else{
      g_sig <- as.matrix(handdat$gamma_prior$sigma)
    }

    g_mean_str <- paste0("[", paste(g_mu, collapse = ", "), "]'")
    g_sigma_rows <- apply(g_sig, 1,
                          function(row) paste0("[",
                                               paste(row, collapse = ", "),
                                               "]"))
    g_sigma_str <- paste0("[", paste(g_sigma_rows, collapse = ", "), "]")

    gamma_prior_line <- paste0("gamma ~ multi_normal(", g_mean_str, ", ",
                               g_sigma_str, ");")
  }else if(is.character(handdat$gamma_prior)){
    gamma_prior_line <- handdat$gamma_prior
  }else{
    stop("gamma_prior must be a list of parameters or a valid Stan character string.")
  }

  # s2 prior
  if(is.list(handdat$s2_prior)){
    s2_prior_line <- paste0("s2 ~ inv_gamma(",
                            handdat$s2_prior$shape, ", ",
                            handdat$s2_prior$scale, ");")
  }else if(is.character(handdat$s2_prior)){
    s2_prior_line <- handdat$s2_prior
  }else{
    stop("s2_prior must be a list of parameters
         or a valid Stan character string.")
  }

  # rho prior
  if(is.list(handdat$rho_prior)){
    rho_prior_line <- paste0("rho ~ normal(",
                             handdat$rho_prior$mean, ", ",
                             handdat$rho_prior$sd, ");")
  }else if(is.character(handdat$rho_prior)){
    rho_prior_line <- handdat$rho_prior
  }else{
    stop("rho_prior must be a list of parameters
         or a valid Stan character string.")
  }

  # Clean template
  stan_template <- paste0("
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
}

transformed data {
  matrix[N, N] L_Omega = cholesky_decompose(Omega);
  matrix[D, D] L_Psi = cholesky_decompose(Psi);
}

parameters {
  vector[K] beta;
  vector[D] gamma;               // FIX: Changed real to vector[D]
  real<lower=0> s2;
  real<lower=-1, upper=1> rho;
  matrix[N, D] Z;
}

transformed parameters {
  real sigma = sqrt(s2);
  matrix[N, D] U = Lambda + L_Omega * Z * L_Psi';
}

model {
  matrix[N, N] I = diag_matrix(rep_vector(1.0, N));
  matrix[N, N] W = I - rho * A;

  // Compute the homophily-adjusted disturbance errors vector (Equation 2)
  vector[N] nu = y - X * beta - U * gamma;
  vector[N] epsilon = W * nu;

  // Priors dynamically inserted via R
  ", beta_prior_line, "
  ", gamma_prior_line, "
  ", s2_prior_line, "
  ", rho_prior_line, "
  to_vector(Z) ~ std_normal();

  // Corrected Likelihood tracking the network spatial transformation matrix W
  target += log_determinant(W);
  epsilon ~ normal(0, sigma);
}
")


  if(!is.null(handdat$compile) && !handdat$compile){
    return(stan_template)
  }

  stan_file_tmp <- cmdstanr::write_stan_file(stan_template)
  compiled_model <- cmdstanr::cmdstan_model(stan_file_tmp)

  stan_data <- list(
    N = handdat$n, K = K, D = handdat$D, y = handdat$y, X = handdat$X,
    A = handdat$A,
    Lambda = handdat$Lambda, Omega = handdat$Omega, Psi = handdat$Psi
  )

  # Extract stan arguments
  all_args <- names(formals(cmdstanr:::CmdStanModel$public_methods$sample))
  valid_sample_args <- setdiff(all_args, "data")
  user_sampling_params <- handdat[names(handdat) %in% valid_sample_args]

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




#' Summary method for homophily adjusted network disturbances objects
#' @method summary hand
#' @export
summary.hand <- function(object, probs = c(0.025, 0.975), ...){
  p <- ncol(object$data$X)
  prob_names <- paste0(probs * 100, "%")

  if(object$method == "stan"){
    draws <- object$estimates

    beta_cols  <- grep("^beta\\[", colnames(draws))
    gamma_cols <- grep("^gamma\\[", colnames(draws))
    s2_col     <- which(colnames(draws) == "s2")
    rho_col    <- which(colnames(draws) == "rho")
    lp_col     <- which(colnames(draws) == "lp__")

    target_cols <- c(beta_cols, gamma_cols, s2_col, rho_col)

    means <- colMeans(draws[, target_cols, drop = FALSE])
    quants <- apply(draws[, target_cols, drop = FALSE], 2,
                    quantile, probs = probs, na.rm = TRUE)

    if(length(target_cols) == 1){
      coef_table <- matrix(c(means, quants), nrow = 1)
    }else{
      coef_table <- cbind(means, t(quants))
    }

    colnames(coef_table) <- c("Estimate", prob_names)

    # FIX: Dynamically generate bracketed names for all D dimensions of gamma
    gamma_names <- paste0("gamma[", 1:length(gamma_cols), "]")
    rownames(coef_table) <- c(colnames(object$data$X), gamma_names, "sigma2", "rho")

    logpost_val <- mean(draws[, lp_col])
    theta_star  <- as.numeric(means)
    loglik_val  <- hand_likelihood(theta_star, object$data$network,
                                   object$data$y, object$data$X,
                                   object$priors$U$Lambda,
                                   object$priors$U$Omega, object$priors$U$Psi,
                                   log=TRUE)
  }else{
    # Normal Approximation
    D_gamma <- length(object$estimates$gamma)
    theta_star <- c(object$estimates$beta, object$estimates$gamma,
                    object$estimates$s2, object$estimates$rho)
    hinv <- tryCatch(solve(-object$estimates$hessian),
                     error = function(e) matrix(NA, length(theta_star),
                                                length(theta_star)))
    std_errs <- sqrt(diag(hinv))

    z_crit <- qnorm(probs)
    intervals <- matrix(NA, nrow = length(theta_star), ncol = length(probs))
    for(i in 1:length(probs)) intervals[,i] <- theta_star + z_crit[i] * std_errs

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
    model        = "hand",
    method       = object$method,
    converged = if(object$method == "normal") object$estimates$converged else T)

  class(res) <- "summary.netreg"
  return(res)
}


#' Print method for hand objects
#' @method print hand
#' @export
print.hand <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  sm <- summary(x)
  coefs <- sm$coefficients[, "Estimate"]

  cat("Coefficients & Parameters:\n")
  print(round(coefs, 5))
  cat("\n")

  invisible(x)
}



# Internal utility to calculate fitted values and residuals for HAND.
# Keeps NetReg::: completely clean of S3 clutter.
calc_hand_diags <- function(object) {
  sm    <- summary(object)
  coefs <- sm$coefficients

  # 1. Extract parameter values cleanly using unified names
  beta  <- coefs[colnames(object$data$X), "Estimate"]

  # Grab gamma dimensions dynamically from row names
  gamma_rows <- grep("^gamma\\[", rownames(coefs), value = TRUE)
  gamma      <- coefs[gamma_rows, "Estimate"]

  X      <- object$data$X
  y      <- object$data$y
  Lambda <- object$priors$U$Lambda

  # 2. Compute diagnostics accounting for latent homophily correction
  fitted_vals <- as.vector(X %*% beta + Lambda %*% gamma)
  resids      <- as.vector(y - fitted_vals)

  return(list(fitted = fitted_vals, residuals = resids))
}

#' Plot method for hand objects
#' @method plot hand
#' @export
plot.hand <- function(x, ...) {
  diags    <- calc_hand_diags(x)
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

