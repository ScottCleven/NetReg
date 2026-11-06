#' Heterogeneous Susceptibilities Network Effects Model
#' @param formula A symbolic formula object mapping response to predictors.
#' @param W_formula A one-sided formula (e.g., ~ x1 + x3) defining susceptibility covariates. Defaults to NULL.
#' @param data A data frame containing the variables in the model. Defaults to NULL.
#' @param network A matrix or graph object representing the network structure.
#' @param method Character string specifying estimation ("normal" or "stan").
#' @param beta_prior List containing prior parameters for beta.
#' @param s2_prior List containing prior parameters for s2.
#' @param gamma_prior List containing prior parameters for the gamma susceptibility coefficients.
#' @param rownorm Logical indicating if network matrices should be row-normalized.
#' @param ... Additional arguments.
#' @import mvtnorm
#' @importFrom truncnorm dtruncnorm
#' @importFrom invgamma dinvgamma
#' @importFrom sna lnam
#' @importFrom cmdstanr cmdstan_model write_stan_file
#' @importFrom stats rgamma rnorm
#' @export
hs_nem <- function(formula, W_formula = NULL, data = NULL, network, method,
                   beta_prior,
                   s2_prior = list(shape = 0.0001, scale = 0.0001),
                   gamma_prior,
                   rownorm = TRUE, compile = TRUE, ...){

  # Extract all individual variables needed across both formulas
  main_vars <- all.vars(formula)
  if(!is.null(W_formula)){
    w_vars <- all.vars(W_formula)
  }else{
    w_vars <- c()
  }
  vars_in_formula <- unique(c(main_vars, w_vars))

  # Initialize lists for scavenger engine
  scavenged_list <- list()
  missing_vars <- c()

  pdat <- convert_network(network, "matrix")
  net_df <- pdat$df
  network <- pdat$A

  # Look for each variable one by one in order of data>environment>network
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
    stop(paste("The following variable(s) required by the model could not be
               found in 'data', the global environment, or the network node
               attributes:", paste(paste0("'", missing_vars, "'"),
                                   collapse = ", ")))
  }

  data <- as.data.frame(scavenged_list)

  # Construct the core data matrices
  X <- model.matrix(formula, data = data)
  y <- model.frame(formula, data)[, 1]
  n <- NROW(y)
  p <- NCOL(X)

  # Set up W matrix
  if(is.null(W_formula)){
    W_tilde <- X
  }else{
    W_tilde <- model.matrix(W_formula, data = data)
  }
  q <- NCOL(W_tilde)

  if(rownorm){
    network <- row_norm(network)
  }

  # Build Adaptive Scale-Invariant Zellner's g-prior for Gamma (g = n)
  if(missing(gamma_prior)){
    message("\nThe gamma_prior parameter is not specified. It will be set automatically based on data scales.")
    W_cross <- t(W_tilde) %*% W_tilde
    # Calculate the scale-invariant Zellner's covariance matrix: n * (W'W)^-1
    Sigma_gamma_mat <- n * solve(W_cross)
    gamma_prior <- list(
      mean = rep(0, q),
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

  # Prior list validation constraints for importance sampling
  is_importance <- tolower(method) %in% c("is", "importance sampling")
  if (is_importance && (typeof(beta_prior) != "list" ||
                        typeof(s2_prior) != "list" ||
                        typeof(gamma_prior) != "list")) {
    stop("importance sampling methods require list input of prior parameters,
         prior functions cannot be changed.")
  }

  hetdat <- list(y=y, X=X, W_tilde=W_tilde, network=network, method=method,
                 p=p, q=q, n=n, I=diag(n), beta_prior=beta_prior,
                 s2_prior=s2_prior, gamma_prior=gamma_prior, compile = compile,
                 ...)

  # Dispatch system
  if(tolower(method) == "normal"){
    hetout <- hs_nem_normal(hetdat)
  }else if(tolower(method) == "stan"){
    hetout <- hs_nem_stan(hetdat)
  }else{
    stop("Non-valid estimation method.
         Heterogeneous models support 'normal' or 'stan'.")
  }

  # Return string immediately if compile is FALSE
  if(tolower(method) == "stan" && !compile){
    return(hetout)
  }

  final_out <- list(
    estimates = hetout,
    formula   = formula,
    W_formula = W_formula,
    model     = "heterogeneous_effects",
    method    = tolower(method),
    data      = list(y = y, X = X, W_tilde = W_tilde, network = network),
    priors    = list(beta = beta_prior, s2 = s2_prior, gamma = gamma_prior),
    call      = match.call()
  )

  class(final_out) <- c("heterogeneous_effects", "netreg")
  return(final_out)
}






# Log-Likelihood for Heterogeneous Effects Model
hs_nem_likelihood <- function(theta, hetdat, log = TRUE){
  p <- hetdat$p
  q <- hetdat$q
  n <- hetdat$n

  # Unpack parameters
  beta  <- theta[1:p]
  s2    <- theta[p + 1]
  gamma <- theta[(p + 2):(p + 1 + q)]

  # Calculate individual rho vector elements: rho = W_tilde %*% gamma
  rho <- as.numeric(hetdat$W_tilde %*% gamma)

  lambda_1 <- max(Re(eigen(hetdat$network, only.values = TRUE)$values))

  # Structural boundary checks: check variance scale and EVERY individual rho_i
  if(s2 <= 0 || any(rho <= -1/lambda_1) || any(rho >= 1/lambda_1)){
    if(log){
      return(-1e10)
    }else{
      return(0)
    }
  }


  # Construct Launtif
  L <- hetdat$I - diag(rho) %*% hetdat$network

  # Verify matrix structural stability via log determinant tracking
  det_L <- det(L)
  if(det_L <= 0){
    if(log){
      return(-1e10)
    }else{
      return(0)
    }
  }


  # Transform y directly to isolate independent error innovations
  # System: L %*% y = X %*% beta + error
  y_trans <- as.numeric(L %*% hetdat$y)
  mu      <- as.numeric(hetdat$X %*% beta)

  # Evaluate Multivariate Normal log density plus spatial Jacobian adjustment
  log_lik <- sum(dnorm(y_trans, mean = mu, sd = sqrt(s2), log = TRUE)) +
    log(det_L)

  if(log){
    return(log_lik)
  }else{
    return(exp(log_lik))
  }
}

# Log-Posterior for Heterogeneous Effects Model
hs_nem_posterior <- function(theta, hetdat, log = TRUE,
                             bprior, s2prior, gprior){
  p <- hetdat$p
  q <- hetdat$q

  # Parameters
  beta  <- theta[1:p]
  s2    <- theta[p + 1]
  gamma <- theta[(p + 2):(p + 1 + q)]

  # Base parameter boundary checks
  if(s2 <= 0){
    if(log){
      return(-Inf)
    }else{
      return(0)
    }
  }

  # Calculate log-likelihood
  log_lik <- hs_nem_likelihood(theta, hetdat, log = TRUE)
  if(!is.finite(log_lik) || log_lik == -1e10){
    return(if(log) -1e10 else 0)
  }

  # Calculate Log Posterior
  log_post <- log_lik + bprior(beta) + s2prior(s2) + gprior(gamma)

  if(log){
    return(log_post)
  }else{
    return(exp(log_post))
  }
}


hs_nem_normal <- function(hetdat){
  p <- hetdat$p
  q <- hetdat$q

  # Beta Prior
  if(typeof(hetdat$beta_prior) == "list"){
    if(length(hetdat$beta_prior$mean) == 1){
      bmean <- rep(hetdat$beta_prior$mean, p)
    }else{
      bmean <- hetdat$beta_prior$mean
    }
    if(length(hetdat$beta_prior$sigma) == 1){
      bsigma <- diag(p) * hetdat$beta_prior$sigma
    }else{
      bsigma <- hetdat$beta_prior$sigma
    }
    bprior <- function(beta){
      mvtnorm::dmvnorm(beta, bmean, bsigma, log = TRUE)
    }
  }else{
    bprior <- hetdat$beta_prior
  }

  # s2 Prior
  if(typeof(hetdat$s2_prior) == "list"){
    s2prior <- function(s2){
      invgamma::dinvgamma(s2, hetdat$s2_prior$shape, hetdat$s2_prior$scale,
                          log = TRUE)
    }
  }else{
    s2prior <- hetdat$s2_prior
  }

  # Gamma Prior
  if(typeof(hetdat$gamma_prior) == "list"){
    if(length(hetdat$gamma_prior$mean) == 1){
      gmean <- rep(hetdat$gamma_prior$mean, q)
    }else{
      gmean <- hetdat$gamma_prior$mean
    }
    if(length(hetdat$gamma_prior$sigma) == 1){
      gsigma <- diag(q) * hetdat$gamma_prior$sigma
    }else{
      gsigma <- hetdat$gamma_prior$sigma
    }
    gprior <- function(gamma){
      mvtnorm::dmvnorm(gamma, gmean, gsigma, log = TRUE)
    }
  }else{
    gprior <- hetdat$gamma_prior
  }


  # Initial value calculation
  init_betas <- as.numeric(lm.fit(hetdat$X, hetdat$y)$coefficients)
  init_gamma <- rep(0, q)
  init_par   <- c(init_betas, 1, init_gamma)

  # Checking if additional Arguments were input
  if(!is.null(hetdat$gr)){
    gr=hetdat$gr
  }else{
    gr=NULL
  }

  if(!is.null(hetdat$lower)){
    lower <- hetdat$lower
  }else{
    lower <- c(rep(-Inf, p), 1e-6, rep(-Inf, q)) # Goes c(beta, sigma, gamma)
  }

  if(!is.null(hetdat$upper)){
    upper <- hetdat$upper
  }else{
    upper <- c(rep(Inf, p), Inf, rep(Inf, q))
  }

  if(!is.null(hetdat$control)){
    control=hetdat$control
  }else{
    control=list(fnscale = -1)
  }

  # Execute optimization loop
  myopt <- optim(par = init_par, fn = hs_nem_posterior, gr=gr,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(fnscale = -1),
    hessian = TRUE,

    hetdat = hetdat,
    log = TRUE,
    bprior = bprior,
    s2prior = s2prior,
    gprior = gprior
  )

  beta_ests <- myopt$par[1:p]
  names(beta_ests) <- colnames(hetdat$X)

  gamma_ests <- myopt$par[(p + 2):(p + 1 + q)]
  names(gamma_ests) <- colnames(hetdat$W_tilde)

  list(
    beta      = beta_ests,
    s2        = myopt$par[p + 1],
    gamma     = gamma_ests,
    rho       = as.numeric(hetdat$W_tilde %*% gamma_ests),
    hessian   = myopt$hessian,
    loglik    = hs_nem_likelihood(myopt$par, hetdat, log = TRUE),
    logpost   = myopt$value,
    counts    = myopt$counts,
    message   = myopt$message,
    converged = myopt$convergence == 0,
    q         = q,
    hetdat    = hetdat
  )
}


hs_nem_stan <- function(hetdat){
  if(!requireNamespace("cmdstanr", quietly = TRUE)){
    stop("Package 'cmdstanr' is needed for this method. Please install it.",
         call. = FALSE)
  }

  if(is.null(hetdat$iterations)){
    iterations <- 2000
  }else{
    iterations <- hetdat$iterations
  }
  if(is.null(hetdat$chains)){
    chains <- 4
  }else{
    chains <- hetdat$chains
  }
  K <- hetdat$p
  Q <- hetdat$q

  # Format beta prior
  if(is.list(hetdat$beta_prior)){
    if(length(hetdat$beta_prior$mean) == 1){
      b_mean_vec <- rep(hetdat$beta_prior$mean, K)
    }else{
      b_mean_vec <- as.vector(hetdat$beta_prior$mean)
    }
    if(length(hetdat$beta_prior$sigma) == 1){
      b_sigma_mat <- diag(K) * hetdat$beta_prior$sigma
    }else{
      b_sigma_mat <- as.matrix(hetdat$beta_prior$sigma)
    }

    b_mean_str <- paste0("[", paste(b_mean_vec, collapse = ", "), "]'")
    b_sigma_rows <- apply(b_sigma_mat, 1,
                          function(row) paste0("[",
                                               paste(row, collapse = ", "),
                                               "]"))
    b_sigma_str <- paste0("[", paste(b_sigma_rows, collapse = ", "), "]")

    beta_prior_line <- paste0("beta ~ multi_normal(", b_mean_str, ", ",
                              b_sigma_str, ");")
  }else if(is.character(hetdat$beta_prior)){
    beta_prior_line <- hetdat$beta_prior
  }else{
    stop("beta_prior must be a list of parameters
         or a valid Stan character string.")
  }

  # Format s2 prior line
  if(is.list(hetdat$s2_prior)){
    s2_prior_line <- paste0("s2 ~ inv_gamma(", hetdat$s2_prior$shape, ", ",
                            hetdat$s2_prior$scale, ");")
  }else if(is.character(hetdat$s2_prior)){
    s2_prior_line <- hetdat$s2_prior
  }else{
    stop("s2_prior must be a list of parameters
         or a valid Stan character string.")
  }

  # Format gamma prior line (Multivariate)
  if(is.list(hetdat$gamma_prior)){
    if(length(hetdat$gamma_prior$mean) == 1){
      g_mean_vec <- rep(hetdat$gamma_prior$mean, Q)
    }else{
      g_mean_vec <- as.vector(hetdat$gamma_prior$mean)
    }
    if(length(hetdat$gamma_prior$sigma) == 1){
      g_sigma_mat <- diag(Q) * hetdat$gamma_prior$sigma
    }else{
      g_sigma_mat <- as.matrix(hetdat$gamma_prior$sigma)
    }

    g_mean_str <- paste0("[", paste(g_mean_vec, collapse = ", "), "]'")
    g_sigma_rows <- apply(g_sigma_mat, 1,
                          function(row) paste0("[",
                                               paste(row, collapse = ", "),
                                               "]"))
    g_sigma_str <- paste0("[", paste(g_sigma_rows, collapse = ", "), "]")

    gamma_prior_line <- paste0("gamma ~ multi_normal(", g_mean_str, ", ",
                               g_sigma_str, ");")
  }else if(is.character(hetdat$gamma_prior)){
    gamma_prior_line <- hetdat$gamma_prior
  }else{
    stop("gamma_prior must be a list of parameters
         or a valid Stan character string.")
  }

  # Uncollapsed template
  stan_template <- paste0("
data {
  int<lower=0> N;                // Number of observations
  int<lower=0> K;                // Number of predictors for X
  int<lower=0> Q;                // Number of predictors for W_tilde
  vector[N] y;                   // Response variable
  matrix[N, K] X;                // Predictor matrix X
  matrix[N, Q] W_tilde;          // Susceptibility matrix W (with intercept)
  matrix[N, N] A;                // Adjacency/Network weights matrix
}

parameters {
  vector[K] beta;
  real<lower=0> s2;
  vector[Q] gamma;      // Estimates elements determining spatial susceptibility
}

transformed parameters {
  real sigma = sqrt(s2);
  vector[N] rho = W_tilde * gamma; //Node-specific spatial autoregressive vector
}

model {
  matrix[N, N] G;

  // Enforce boundary constraints using a scalar loop and negative infinity
  for (i in 1:N) {
    if (rho[i] <= -1.0 || rho[i] >= 1.0) {
      target += negative_infinity();
    }
  }

  // Construct G = I - diag(rho) * A row by row
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        G[i, j] = 1.0 - rho[i] * A[i, j];
      } else {
        G[i, j] = -rho[i] * A[i, j];
      }
    }
  }

  // NEM Form: Transform y directly
  vector[N] y_transformed = G * y;
  vector[N] mu = X * beta;

  // Priors dynamically inserted via R
  ", beta_prior_line, "
  ", s2_prior_line, "
  ", gamma_prior_line, "

  // Likelihood with spatial Jacobian correction
  target += log_determinant(G);
  y_transformed ~ normal(mu, sigma);
}
")

  if(!is.null(hetdat$compile) && !hetdat$compile){
    return(stan_template)
  }

  stan_file_tmp <- cmdstanr::write_stan_file(stan_template)
  compiled_model <- cmdstanr::cmdstan_model(stan_file_tmp)

  stan_data <- list(
    N        = hetdat$n,
    K        = K,
    Q        = Q,
    y        = hetdat$y,
    X        = hetdat$X,
    W_tilde  = hetdat$W_tilde,
    A        = hetdat$network
  )

  # Extract stan arguments
  all_args <- names(formals(cmdstanr:::CmdStanModel$public_methods$sample))
  valid_sample_args <- setdiff(all_args, "data")
  user_sampling_params <- dat[names(hetdat) %in% valid_sample_args]

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



#' Summary method for heterogeneous susceptibilities effects objects
#' @method summary heterogeneous_effects
#' @param object An object of class "heterogeneous_effects"
#' @param probs A numeric vector of probabilities for the interval estimates.
#' @param ... Additional arguments
#' @export
summary.heterogeneous_effects <- function(object, probs = c(0.025, 0.975), ...){
  p <- ncol(object$data$X)
  q <- ncol(object$data$W_tilde)
  n <- length(object$data$y)
  prob_names <- paste0(probs * 100, "%")

  if(object$method == "stan"){
    draws <- object$estimates

    beta_cols  <- grep("^beta\\[", colnames(draws))
    s2_col     <- honest_idx <- which(colnames(draws) == "s2")
    gamma_cols <- grep("^gamma\\[", colnames(draws))
    rho_cols   <- grep("^rho\\[", colnames(draws))
    lp_col     <- which(colnames(draws) == "lp__")

    target_cols <- c(beta_cols, s2_col, gamma_cols)

    means <- colMeans(draws[, target_cols, drop = FALSE])
    quants <- apply(draws[, target_cols, drop = FALSE], 2,
                    quantile, probs = probs, na.rm = TRUE)

    if(length(target_cols) == 1){
      coef_table <- matrix(c(means, quants), nrow = 1)
    }else{
      coef_table <- cbind(means, t(quants))
    }

    colnames(coef_table) <- c("Estimate", prob_names)
    rownames(coef_table)[1:p] <- colnames(object$data$X)
    rownames(coef_table)[(p + 2):(p + 1 + q)] <- paste0("gamma:",
                                                 colnames(object$data$W_tilde))
    rownames(coef_table)[p + 1] <- "sigma2"

    # Process rho metrics out of full sampling chains
    rho_draws <- draws[, rho_cols, drop = FALSE]
    rho_means <- colMeans(rho_draws)
    rho_summary <- matrix(c(mean(rho_means), sd(rho_means), min(rho_means),
                            max(rho_means)), nrow = 1)
    colnames(rho_summary) <- c("Mean", "SD", "Min", "Max")
    rownames(rho_summary) <- "Individual rho"

    logpost_val <- mean(draws[, lp_col])
    theta_star  <- as.numeric(means)
    loglik_val  <- hs_nem_likelihood(theta_star, object$estimates, log=TRUE)

  }else{
    # Normal Approximation / Optimization method
    theta_star <- c(object$estimates$beta, object$estimates$s2,
                    object$estimates$gamma)
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
    rownames(coef_table) <- c(colnames(object$data$X), "sigma2",
                              paste0("gamma:", colnames(object$data$W_tilde)))

    # Process rho vector distribution from point estimates
    rho_vec <- object$estimates$rho
    rho_summary <- matrix(c(mean(rho_vec), sd(rho_vec), min(rho_vec),
                            max(rho_vec)), nrow = 1)
    colnames(rho_summary) <- c("Mean", "SD", "Min", "Max")
    rownames(rho_summary) <- "Individual rho"

    loglik_val  <- object$estimates$loglik
    logpost_val <- object$estimates$logpost
  }

  res <- list(
    call         = object$call,
    coefficients = coef_table,
    rho_summary  = rho_summary,
    loglik       = loglik_val,
    logpost      = logpost_val,
    marg_lik     = NA,
    model        = "heterogeneous_effects",
    method       = object$method,
    converged = if(object$method == "normal") object$estimates$converged else T
  )

  class(res) <- "summary.netreg"
  return(res)
}


#' Print method for heterogeneous_effects objects
#' @method print heterogeneous_effects
#' @export
print.heterogeneous_effects <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  if (x$method == "stan") {
    sm    <- summary(x)
    coefs <- sm$coefficients[, "Estimate"]
  } else { # normal
    coefs <- c(x$estimates$beta, sigma2 = x$estimates$s2, x$estimates$gamma)
    names(coefs) <- c(colnames(x$data$X), "sigma2", paste0("gamma:", colnames(x$data$W_tilde)))
  }

  cat("Coefficients & Parameters:\n")
  print(round(coefs, 5))
  cat("\n")

  # Inject the summary profile of individual rho variations
  sm_obj <- summary(x)
  cat("Summary of Individual Network Influences (rho):\n")
  print(round(sm_obj$rho_summary, 5))
  cat("\n")

  invisible(x)
}


# Internal utility to calculate fitted values and residuals for HS-NEM.
# Keeps NetReg::: completely clean of S3 clutter.
calc_heterogeneous_diags <- function(object) {
  # 1. Extract beta and gamma parameters based on the estimation method
  if (object$method == "stan") {
    sm    <- summary(object)
    coefs <- sm$coefficients
    beta  <- coefs[colnames(object$data$X), "Estimate"]

    gamma_names <- paste0("gamma:", colnames(object$data$W_tilde))
    gamma <- coefs[gamma_names, "Estimate"]
  } else { # normal
    beta  <- object$estimates$beta
    gamma <- object$estimates$gamma
  }

  X       <- object$data$X
  y       <- object$data$y
  W_tilde <- object$data$W_tilde
  W       <- object$data$network

  # 2. Reconstruct individual network influence parameter vector (rho_i)
  rho_vec <- 1 / (1 + exp(-as.vector(W_tilde %*% gamma)))

  # 3. Compute structural fitted values and raw residuals
  fitted_vals <- as.vector(X %*% beta + rho_vec * (W %*% y))
  resids      <- as.vector(y - fitted_vals)

  return(list(fitted = fitted_vals, residuals = resids))
}

#' Plot method for heterogeneous_effects objects
#' @method plot heterogeneous_effects
#' @export
plot.heterogeneous_effects <- function(x, ...) {
  diags    <- calc_heterogeneous_diags(x)
  fit_vals <- diags$fitted
  res_vals <- diags$residuals

  # Minimal diagnostic scatter plot inheriting base margins
  plot(fit_vals, res_vals,
       xlab = "Fitted values",
       ylab = "Residuals",
       main = "Residuals vs Fitted",
       panel.first = grid(),
       pch = 20, col = "darkgray", ...)

  # Horizontal reference line at 0
  abline(h = 0, lty = 3, col = "gray40", lwd = 1.5)

  # Solid black local-regression trend line
  if (length(res_vals) > 3) {
    lines(lowess(fit_vals, res_vals), col = "black", lwd = 1.5)
  }

  invisible(x)
}



