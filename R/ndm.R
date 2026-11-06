#' ndm - A function to perform the Network Disturbances Model
#' @param formula an object of class "formula" (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible
#' by as.data.frame to a data frame) containing the variables in the model.
#' If not found in data, the variables are taken from environment(formula),
#' typically the environment from which nam is called. You can also input an
#' igraph object and the dataset and adjacency matrix will be extracted for
#' you (in which case the network argument will be ignored).
#' @param network An adjacency matrix describing your network. Can be of the class:
#' c("matrix", "dgCMatrix", "igraph", "network") or an edgelist.
#' @param method A character indicating the method to use for simulation.
#' Options vary by model but all potential options are c("normal", "IS",
#' "stan", "lnam"). Run ?methods.nam() for model-specific designations.
#' @param beta_prior The prior for the beta vector of the model. Requires list
#' or function inputs. Default is a multivariate normal with mean=0 and sigma=1
#' If you supply a list, you must have one element named 'mean' which is the
#' mean of the multivariate normal and one named 'sigma' which is the
#' variance-covariance matrix. If mean is a single integer,
#' that mean will be used for every element in the mean vector and if sigma is
#' a single integer, it will multiply that integer by the identity matrix.
#' Otherwise, both mean and sigma will be what you specify them to be. Or the
#' user can specify their own prior function with which to input. If
#' method="stan" only list input will be accepted.
#' @param s2_prior Prior for the variance parameter. Requires list or function
#' inputs. Defaults to a gamma(2,0.5). If Input is a list, the list must have
#' two elements named 'shape' and 'scale' to represent the shape and scale of
#' the gamma distribution respectively. Or the user can specify their own prior
#' function with which to input. If method="stan" only list input
#' will be accepted.
#' @param rho_prior The prior for the correlation of the adjacency matrix of
#' the model. Requires list or function inputs. Default is a N(0.36,0.49)
#' (Dittrich et al. 2017). If input is a list, the list must have two elements
#' named 'mean' and 'sd' which are the mean and standard deviation of the normal
#' distribution respectively. Or the user can specify their own prior
#' function with which to input. If method="stan" only list input
#' will be accepted.
#' @param rownorm Default to TRUE. Will row-normalize the adjacency matrix.
#' @param ... Additional arguments to be passed into the optim() function for
#' when method="normal".
#' @return A list with two classes, the first being "nam" and the second is the
#' class of "disturbances" that contains values:
#' @return \code{estimates} The raw estimate output of the chosen method.
#' @return \code{formula} The inputted formula
#' @return \code{model} The chosen model ("disturbances").
#' @return \code{method} The chosen method.
#' @return \code{data} The inputted data of the form list(y,X,network).
#' @return \code{priors} The inputted prior settings/functions.
#' @return \code{call} The match.call() of what was inputted.
#'
#' @import mvtnorm
#' @importFrom truncnorm dtruncnorm
#' @importFrom invgamma dinvgamma
#' @importFrom sna lnam
#' @importFrom cmdstanr cmdstan_model write_stan_file
#' @importFrom stats rgamma rnorm
#' @export


# Model-specific wrapper function for Network Disturbances Model
ndm <- function(formula, data=NULL, network, method,
                beta_prior,
                s2_prior = list(shape=0.0001, scale=0.0001),
                rho_prior = list(mean=0.36, sd=0.7),
                rownorm=TRUE, compile = TRUE, ...){

  # Extract all individual variables needed by the formula
  vars_in_formula <- all.vars(formula)

  # Initialize a list to hold the scavenged vectors
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

  # Consolidate all columns into one data frame
  data <- as.data.frame(scavenged_list)

  X = model.matrix(formula, data=data)
  y = model.frame(formula, data)[,1]
  n = NROW(y)
  p = NCOL(X)

  if(rownorm){
    network = row_norm(network)
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

  # Check importance sampling inputs here.
  is_importance <- tolower(method) %in% c("is", "importance sampling")
  if(is_importance && (typeof(beta_prior) != "list" ||
                        typeof(s2_prior) != "list" ||
                        typeof(rho_prior) != "list")){
    stop("importance sampling methods require list input of prior parameters,
         prior functions cannot be changed.")
  }


  ndmdat <- list(y=y, X=X, network=network, method=method, p=p, n=n,
                 I=diag(n), beta_prior=beta_prior, s2_prior = s2_prior,
                 rho_prior = rho_prior, compile = compile, ...)

  # Dispatch to NDM-specific estimation functions
  if(tolower(method) == "normal"){
    ndmout <- ndm_normal(ndmdat)
  }else if(tolower(method) == "stan"){
    ndmout <- ndm_stan(ndmdat)
  }else if(tolower(method) == "lnam"){
    ndmout <- ndm_lnam(ndmdat)
  }else if(tolower(method) %in% c("is", "importance sampling")){
    ndmout <- ndm_IS(ndmdat)
  }else{
    stop("Unknown estimation method.")
  }

  # Return string immediately if compile is FALSE
  if(tolower(method) == "stan" && !compile){
    return(ndmout)
  }

  final_out <- list(
    estimates = ndmout,
    formula   = formula,
    model     = "disturbances",
    method    = tolower(method),
    data      = list(y = y, X = X, network = network),
    priors    = list(beta = beta_prior, s2 = s2_prior, rho = rho_prior),
    call      = match.call()
  )

  # Assign class
  class(final_out) <- c("disturbances", "netreg")
  return(final_out)
}



# Likelihood of the disturbances model.
ndm_likelihood <- function(theta, A, y, X, log=TRUE){
  # Extract parameters
  beta <- theta[1:(length(theta)-2)]
  s2 <- theta[length(theta)-1]
  rho <- theta[length(theta)]

  n <- NROW(y)
  I <- diag(n)

  # Matrix (I - rho*A)
  L <- I - rho * A

  mu <- X %*% beta

  L_inv <- qr.solve(L)
  V <- s2 * tcrossprod(L_inv)

  mvtnorm::dmvnorm(y, mean = mu, sigma = V, log = log)
}


# Posterior of the disturbances model.
ndm_posterior <- function(theta, A, y, X, log=TRUE,
                          bprior, s2prior, rhoprior){
  # Extract parameters
  beta <- theta[1:(length(theta)-2)]
  s2 <- theta[length(theta)-1]
  rho <- theta[length(theta)]

  # Parameter boundary checks to prevent optimization crashes
  if(s2 <= 0 || rho <= -1 || rho >= 1){
    if(log){
      return(-Inf)
    }else{
      return(0)
    }
  }

  # Calculate Likelihood using the NDM logic (mean = Xb)
  log_lik <- ndm_likelihood(theta, A, y, X, log=TRUE)

  # Calculate log-posterior by adding log-prior densities
  log_post <- log_lik + bprior(beta) + s2prior(s2) + rhoprior(rho)

  if(log){
    return(log_post)
  }else{
    return(exp(log_post))
  }
}



# Code for when method="normal"
ndm_normal <- function(ndmdat){

  # Checking Prior Inputs
  if(typeof(ndmdat$beta_prior) == "list"){
    if(length(ndmdat$beta_prior$mean) == 1){
      bmean <- rep(ndmdat$beta_prior$mean, ndmdat$p)
    }else{
      bmean <- ndmdat$beta_prior$mean
    }
    if(length(ndmdat$beta_prior$sigma) == 1){
      bsigma <- diag(ndmdat$p) * ndmdat$beta_prior$sigma
    }else{
      bsigma <- ndmdat$beta_prior$sigma
    }
    bprior <- function(beta){
      mvtnorm::dmvnorm(beta, bmean, bsigma, log=TRUE)
    }
  }else{
    bprior <- ndmdat$beta_prior
  }

  if(typeof(ndmdat$s2_prior) == "list"){
    s2prior <- function(s2){
      invgamma::dinvgamma(s2, ndmdat$s2_prior$shape, ndmdat$s2_prior$scale,
                          log=TRUE)
    }
  }else{
    s2prior <- ndmdat$s2_prior
  }

  if(typeof(ndmdat$rho_prior) == "list"){
    rhoprior <- function(rho){
      log(truncnorm::dtruncnorm(rho, a=-1, b=1, mean=ndmdat$rho_prior$mean,
                                sd=ndmdat$rho_prior$sd))
    }
  }else{
    rhoprior <- ndmdat$rho_prior
  }

  # Checking if additional Arguments were input
  if(!is.null(ndmdat$gr)){
    gr=ndmdat$gr
  }else{
    gr=NULL
  }

  if(!is.null(ndmdat$lower)){
    lower <- ndmdat$lower
  }else{
    lower <- c(rep(-Inf, ndmdat$p), 1e-6, -0.999) # Safeguards s2 and rho
  }

  if(!is.null(ndmdat$upper)){
    upper <- ndmdat$upper
  }else{
    upper <- c(rep(Inf, ndmdat$p), Inf, 0.999)
  }

  if(!is.null(ndmdat$control)){
    control=ndmdat$control
  }else{
    control=list(fnscale = -1)
  }

  # Finding good initial beta values
  init.betas <- as.numeric(lm.fit(ndmdat$X, ndmdat$y)$coefficients)

  # Finding posterior estimates using the normal approximation of the posterior
  myopt <- optim(par = c(init.betas, 1, 0.36), fn = ndm_posterior, gr=gr,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = control,
    hessian = TRUE,

    A = ndmdat$network,
    y = ndmdat$y,
    X = ndmdat$X,
    log = TRUE,
    bprior = bprior,
    s2prior = s2prior,
    rhoprior = rhoprior
  )

  beta_ests <- myopt$par[1:ndmdat$p]
  names(beta_ests) <- colnames(ndmdat$X)

  # Outputting the results
  list(
    beta      = beta_ests,
    s2        = myopt$par[ndmdat$p + 1],
    rho       = myopt$par[ndmdat$p + 2],
    hessian   = myopt$hessian,
    loglik    = ndm_likelihood(myopt$par, ndmdat$network,
                               ndmdat$y, ndmdat$X, log=TRUE),
    logpost   = myopt$value,
    counts    = myopt$counts,
    message   = myopt$message,
    converged = myopt$convergence == 0
  )
}



# Code for when method="IS"
ndm_IS <- function(ndmdat){
  if(is.null(ndmdat$samples)){
    ndmdat$samples = 10000
  }
  samples <- ndmdat$samples
  p <- ndmdat$p

  # Run fast Normal optimizer first to find the true posterior peak
  fit_norm <- ndm_normal(ndmdat)

  mu_proposal <- c(fit_norm$beta, s2 = fit_norm$s2, rho = fit_norm$rho)

  # Extract the covariance matrix from the inverted negative Hessian
  # Add ridge inflation (* 2) to ensure the proposal is wider than the target
  cov_proposal <- solve(-fit_norm$hessian) * 2

  # Draw proposals from the calibrated Multivariate Normal distribution
  proposal_samples <- mvtnorm::rmvnorm(samples, mean = mu_proposal,
                                       sigma = cov_proposal)

  beta_samples <- proposal_samples[, 1:p]
  s2_samples   <- proposal_samples[, p + 1]
  rho_samples  <- proposal_samples[, p + 2]

  # Setup standard prior closures for weighting calculations
  if(length(ndmdat$beta_prior$mean) == 1){
    b_mean <- rep(ndmdat$beta_prior$mean, p)
  }else{
    b_mean <- ndmdat$beta_prior$mean
  }
  if(length(ndmdat$beta_prior$sigma) == 1){
    b_sigma <- diag(p) * ndmdat$beta_prior$sigma
  }else{
    b_sigma <- ndmdat$beta_prior$sigma
  }

  bprior <- function(b){
    mvtnorm::dmvnorm(b, b_mean, b_sigma, log = TRUE)
  }
  s2prior <- function(s){
    invgamma::dinvgamma(s, ndmdat$s2_prior$shape, ndmdat$s2_prior$scale,
                        log = TRUE)
  }
  rhoprior <- function(r){
    if(r <= -1 || r >= 1){
      -Inf
    }else{
      log(truncnorm::dtruncnorm(r, a=-1, b=1, mean=ndmdat$rho_prior$mean,
                                sd=ndmdat$rho_prior$sd))
    }
  }
  # Calculate Normalized Importance Weights
  # Weight = [ Likelihood * Prior ] / [ Proposal_Density ]
  log_weights <- numeric(samples)

  for(i in 1:samples){
    if(s2_samples[i] <= 0 || rho_samples[i] <= -1 || rho_samples[i] >= 1){
      log_weights[i] <- -Inf
      next
    }

    theta <- c(beta_samples[i,], s2_samples[i], rho_samples[i])

    log_target <- ndm_likelihood(theta, ndmdat$network, ndmdat$y, ndmdat$X,
                                 log = TRUE) + bprior(beta_samples[i,]) +
                  s2prior(s2_samples[i]) + rhoprior(rho_samples[i])

    # Proposal denominator log-density
    log_proposal_dens <- mvtnorm::dmvnorm(theta, mean = mu_proposal,
                                          sigma = cov_proposal, log = TRUE)

    log_weights[i] <- log_target - log_proposal_dens
  }

  # Normalize using stable log-sum-exp
  max_log_w <- max(log_weights[is.finite(log_weights)])
  weights <- exp(log_weights - max_log_w)
  weights[is.na(weights)] <- 0
  weights <- weights / sum(weights)

  # Resample to form the empirical posterior distribution
  indices <- sample(1:samples, size = samples, replace = TRUE, prob = weights)

  ess <- 1 / sum(weights^2)

  # Compare ESS to 5% of samples
  ess_threshold <- 0.05 * samples
  is_converged <- ess >= ess_threshold

  if (!is_converged){
    warning(paste0("Effective Sample Size (ESS = ", round(ess, 1), ") is below
                   5% of nominal samples. Consider increasing nominal samples or
                   checking your model bounds."))
  }

  resampled_beta <- beta_samples[indices, , drop = FALSE]
  resampled_s2   <- s2_samples[indices]
  resampled_rho  <- rho_samples[indices]

  theta_mean <- c(colMeans(resampled_beta), mean(resampled_s2),
                  mean(resampled_rho))

  list(
    beta = colMeans(resampled_beta),
    s2 = as.numeric(mean(resampled_s2)),
    rho = as.numeric(mean(resampled_rho)),
    samples = list(beta = resampled_beta, s2 = resampled_s2,
                   rho = resampled_rho),
    logpost = ndm_posterior(theta_mean, ndmdat$network, ndmdat$y, ndmdat$X,
                            bprior = bprior, s2prior = s2prior,
                            rhoprior = rhoprior),
    loglik = ndm_likelihood(theta_mean, ndmdat$network, ndmdat$y, ndmdat$X,
                            log = TRUE),
    log_marginal_lik = max_log_w + log(mean(exp(log_weights - max_log_w))),
    ess = ess,
    converged = is_converged
  )

}



# Code for when method="stan"
ndm_stan <- function(ndmdat){
  # Checking cmdstanr
  if(!requireNamespace("cmdstanr", quietly = TRUE)){
    stop("Package 'cmdstanr' is needed for this method. Please install it.",
         call. = FALSE)
  }

  # Checking for additional inputs
  if(is.null(ndmdat$iterations)){
    iterations <- 2000
  }else{
    iterations <- ndmdat$iterations
  }

  if(is.null(ndmdat$chains)){
    chains <- 4
  }else{
    chains <- ndmdat$chains
  }

  K <- ndmdat$p

  # Format beta prior line
  if(is.list(ndmdat$beta_prior)){
    if(length(ndmdat$beta_prior$mean) == 1){
      b_mean_vec <- rep(ndmdat$beta_prior$mean, K)
    }else{
      b_mean_vec <- as.vector(ndmdat$beta_prior$mean)
    }
    if(length(ndmdat$beta_prior$sigma) == 1){
      b_sigma_mat <- diag(K) * ndmdat$beta_prior$sigma
    }else{
      b_sigma_mat <- as.matrix(ndmdat$beta_prior$sigma)
    }

    b_mean_str <- paste0("[", paste(b_mean_vec, collapse = ", "), "]'")
    b_sigma_rows <- apply(b_sigma_mat, 1,
                          function(row) paste0("[",
                                               paste(row, collapse = ", "),
                                               "]"))
    b_sigma_str <- paste0("[", paste(b_sigma_rows, collapse = ", "), "]")

    beta_prior_line <- paste0("beta ~ multi_normal(", b_mean_str, ", ",
                              b_sigma_str, ");")
  }else if(is.character(ndmdat$beta_prior)){
    beta_prior_line <- ndmdat$beta_prior
  }else{
    stop("beta_prior must be a list of parameters
         or a valid Stan character string.")
  }

  # Format s2 prior line
  if(is.list(ndmdat$s2_prior)){
    s2_prior_line <- paste0("s2 ~ inv_gamma(", ndmdat$s2_prior$shape, ", ",
                            ndmdat$s2_prior$scale, ");")
  }else if(is.character(ndmdat$s2_prior)){
    s2_prior_line <- ndmdat$s2_prior
  }else{
    stop("s2_prior must be a list of parameters
         or a valid Stan character string.")
  }

  # Format rho prior line
  if(is.list(ndmdat$rho_prior)){
    rho_prior_line <- paste0("rho ~ normal(", ndmdat$rho_prior$mean, ", ",
                             ndmdat$rho_prior$sd, ");")
  }else if(is.character(ndmdat$rho_prior)){
    rho_prior_line <- ndmdat$rho_prior
  }else{
    stop("rho_prior must be a list of parameters
         or a valid Stan character string.")
  }

  # Clean template block with exact equations and structure preserved
  stan_template <- paste0("
data {
  int<lower=0> N;
  int<lower=0> K;
  vector[N] y;
  matrix[N, K] X;
  matrix[N, N] A;
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

  // Priors dynamically inserted via R
  ", beta_prior_line, "
  ", s2_prior_line, "
  ", rho_prior_line, "

  // Likelihood
  target += log_determinant(W);
  epsilon ~ normal(0, sigma);
}
")

  if(!is.null(ndmdat$compile) && !ndmdat$compile){
    return(stan_template)
  }

  stan_file_tmp <- cmdstanr::write_stan_file(stan_template)
  compiled_model <- cmdstanr::cmdstan_model(stan_file_tmp)

  stan_data <- list(
    N = length(ndmdat$y),
    K = ncol(ndmdat$X),
    y = ndmdat$y,
    X = ndmdat$X,
    A = ndmdat$network
  )

  # Extract stan arguments
  all_args <- names(formals(cmdstanr:::CmdStanModel$public_methods$sample))
  valid_sample_args <- setdiff(all_args, "data")
  user_sampling_params <- dat[names(nemdat) %in% valid_sample_args]

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



# Code for when method="lnam"
ndm_lnam <- function(ndmdat){

  # Checking if additional arguments were input
  if(!is.null(nemdat$theta.seed)){
    theta.seed=nemdat$theta.seed
  }else{
    theta.seed=NULL
  }

  if(!is.null(nemdat$null.model)){
    null.model = nemdat$null.model
  }else{
    null.model="meanstd"
  }

  if(!is.null(nemdat$control)){
    control=nemdat$control
  }else{
    control=list()
  }

  if(!is.null(nemdat$tol)){
    tol=nemdat$tol
  }else{
    tol=1e-10
  }

  # Running lnam for the Disturbances model
  res <- sna::lnam(
    y = ndmdat$y,
    x = ndmdat$X,
    W1 = NULL,
    W2 = ndmdat$network,
    theta.seed = theta.seed,
    null.model = null.model,
    method = "BFGS",
    control = control,
    tol = tol
  )

  list(
    beta = res$beta,
    s2   = res$sigma.sq,
    rho  = res$rho2,
    raw  = res
  )
}

#' Summary method for disturbances objects
#' @method summary disturbances
#' @param object An object of class "disturbances"
#' @param probs A numeric vector of probabilities for the interval estimates.
#' @param ... Additional arguments
#' @export
summary.disturbances <- function(object, probs = c(0.025, 0.975), ...){
  p <- ncol(object$data$X)
  loglik_val <- NA
  logpost_val <- NA
  marg_lik <- NA
  prob_names <- paste0(probs * 100, "%")

  if(object$method == "stan"){
    draws <- object$estimates

    beta_cols <- grep("^beta\\[", colnames(draws))
    s2_col    <- which(colnames(draws) == "s2")
    rho_col   <- which(colnames(draws) == "rho")
    lp_col    <- which(colnames(draws) == "lp__")

    target_cols <- c(beta_cols, s2_col, rho_col)

    means <- colMeans(draws[, target_cols, drop = FALSE])
    quants <- apply(draws[, target_cols, drop = FALSE], 2,
                    quantile, probs = probs, na.rm = TRUE)

    if(length(target_cols) == 1){
      coef_table <- matrix(c(means, quants), nrow = 1)
    }else{
      coef_table <- cbind(means, t(quants))
    }

    colnames(coef_table) <- c("Estimate", prob_names)
    rownames(coef_table) <- c(colnames(object$data$X), "sigma2", "rho")

    logpost_val <- mean(draws[, lp_col])
    theta_star <- as.numeric(means)
    loglik_val <- ndm_likelihood(theta_star, object$data$network,
                                 object$data$y, object$data$X, log=TRUE)

  }else if(object$method %in% c("importance", "is", "importance sampling")){
    samples <- object$estimates$samples
    theta_samples <- cbind(samples$beta, sigma2 = samples$s2, rho = samples$rho)

    coef_table <- t(apply(theta_samples, 2, function(x){
      c(Estimate = mean(x), quantile(x, probs = probs))
    }))
    rownames(coef_table)[1:p] <- colnames(object$data$X)

    loglik_val <- object$estimates$loglik
    logpost_val <- object$estimates$logpost
    marg_lik <- object$estimates$log_marginal_lik

  }else{
    if(object$method == "normal"){
      theta_star <- c(object$estimates$beta, object$estimates$s2,
                      object$estimates$rho)
      hinv <- tryCatch(solve(-object$estimates$hessian),
                       error = function(e) matrix(NA, length(theta_star),
                                                  length(theta_star)))
      std_errs <- sqrt(diag(hinv))
      loglik_val <- object$estimates$loglik
      logpost_val <- object$estimates$logpost
    }else{
      est <- object$estimates$raw
      theta_star <- c(est$beta, est$sigma.sq, est$rho2)
      std_errs <- c(est$se.beta, est$se.sigma.sq, est$se.rho2)
      loglik_val <- est$loglik
    }

    z_crit <- qnorm(probs)
    intervals <- matrix(NA, nrow = length(theta_star), ncol = length(probs))
    for(i in 1:length(probs)) {
      intervals[,i] <- theta_star + z_crit[i] * std_errs
    }
    coef_table <- cbind(Estimate = theta_star, intervals)
    colnames(coef_table) <- c("Estimate", prob_names)
    rownames(coef_table) <- c(colnames(object$data$X), "sigma2", "rho")
  }

  res <- list(
    call = object$call,
    coefficients = coef_table,
    loglik = loglik_val,
    logpost = logpost_val,
    marg_lik = marg_lik,
    model = object$model,
    method = object$method,
    converged = if(object$method == "normal") object$estimates$converged else T
  )

  class(res) <- "summary.netreg"
  return(res)
}


#' Print method for disturbances objects
#' @method print disturbances
#' @export
print.disturbances <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  if (x$method %in% c("stan", "importance", "is", "importance sampling")) {
    sm <- summary(x)
    coefs <- sm$coefficients[, "Estimate"]
  } else if (x$method == "normal") {
    coefs <- c(x$estimates$beta, sigma2 = x$estimates$s2, rho = x$estimates$rho)
    names(coefs)[1:ncol(x$data$X)] <- colnames(x$data$X)
  } else { # lnam
    est <- x$estimates$raw
    coefs <- c(est$beta, sigma2 = est$sigma.sq, rho = est$rho2)
    names(coefs)[1:ncol(x$data$X)] <- colnames(x$data$X)
  }

  cat("Coefficients & Parameters:\n")
  print(round(coefs, 5))
  cat("\n")
  invisible(x)
}


# Internal utility to calculate fitted values and residuals for NDM.
# This is a standard hidden function—not an S3 method—to keep NetReg::: clean.
calc_disturbances_diags <- function(object) {
  # 1. Extract beta and rho estimates based on the method
  if (object$method %in% c("stan", "importance", "is", "importance sampling")) {
    sm    <- summary(object)
    coefs <- sm$coefficients
    beta  <- coefs[colnames(object$data$X), "Estimate"]
  } else if (object$method == "normal") {
    beta  <- object$estimates$beta
  } else { # lnam
    beta  <- object$estimates$raw$beta
  }

  # 2. Compute diagnostics for a Network Disturbance structure
  X           <- object$data$X
  y           <- object$data$y

  fitted_vals <- as.vector(X %*% beta)
  resids      <- as.vector(y - fitted_vals)

  return(list(fitted = fitted_vals, residuals = resids))
}

#' Plot method for disturbances objects
#' @method plot disturbances
#' @export
plot.disturbances <- function(x, ...) {
  diags    <- calc_disturbances_diags(x)
  fit_vals <- diags$fitted
  res_vals <- diags$residuals

  # Clean, minimal plot inheriting active user/base margins
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

