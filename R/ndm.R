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
#' @import rstan
#' @import Rcpp
#' @importFrom truncnorm dtruncnorm
#' @importFrom invgamma dinvgamma
#' @importFrom sna lnam
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom stats rgamma rnorm
#' @useDynLib nafm, .registration = TRUE
#' @export


# Model-specific wrapper function for Network Disturbances Model
ndm <- function(formula, data=NULL, network, method,
                beta_prior = list(mean=0, sigma = 1),
                s2_prior = list(shape=2, scale=0.5),
                rho_prior = list(mean=0.36, sd=0.7^2),
                rownorm=TRUE, ...){

  X = model.matrix(formula, data=data)
  y = model.frame(formula, data)[,1]
  n = NROW(y)
  p = NCOL(X)

  network = convert_adjacency(network, "matrix")$A
  if(rownorm){
    network = row_norm(network)
  }

  # Check stan prior inputs
  if(tolower(method) == "stan" &&
     (!is.list(beta_prior) | !is.list(s2_prior) | !is.list(rho_prior))){
    stop("stan method requires list input of prior parameters.")
  }

  ndmdat <- list(y=y, X=X, network=network, method=method, p=p, n=n,
                 I=diag(n), beta_prior=beta_prior, s2_prior = s2_prior,
                 rho_prior = rho_prior, ...)

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

  final_out <- list(
    estimates = ndmout,
    formula   = formula,
    model     = "disturbances",
    method    = tolower(method),
    data      = list(y = y, X = X, network = network),
    priors    = list(beta = beta_prior, s2 = s2_prior, rho = rho_prior),
    call      = match.call()
  )

  # Assign S3 classes for disturbances
  class(final_out) <- c("nam", "disturbances")
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


# Posterior of the effects model.
ndm_posterior <- function(theta, A, y, X, log=TRUE, bprior, s2prior, rhoprior){
  # Extract parameters
  beta <- theta[1:(length(theta)-2)]
  s2 <- theta[length(theta)-1]
  rho <- theta[length(theta)]

  # Calculate Likelihood using the NDM logic (mean = Xb)
  lik <- ndm_likelihood(theta, A, y, X, log=FALSE)

  # Calculate Posterior
  out <- lik * bprior(beta) * s2prior(s2) * rhoprior(rho)

  if(log){
    return(log(out))
  }

  return(out)
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
      mvtnorm::dmvnorm(beta, bmean, bsigma)
    }
  }else{
    bprior <- ndmdat$beta_prior
  }

  if(typeof(ndmdat$s2_prior) == "list"){
    s2prior <- function(s2){
      invgamma::dinvgamma(s2, ndmdat$s2_prior$shape, ndmdat$s2_prior$scale)
    }
  }else{
    s2prior <- ndmdat$s2_prior
  }

  if(typeof(ndmdat$rho_prior) == "list"){
    rhoprior <- function(rho){
      truncnorm::dtruncnorm(rho, a=-1, b=1, mean=ndmdat$rho_prior$mean, sd=ndmdat$rho_prior$sd)
    }
  }else{
    rhoprior <- ndmdat$rho_prior
  }

  # Checking if additional Arguments were input
  if(!is.null(nemdat$gr)){
    gr=nemdat$gr
  }else{
    gr=NULL
  }

  if(!is.null(nemdat$lower)){
    lower = nemdat$lower
  }else{
    lower=-Inf
  }

  if(!is.null(nemdat$upper)){
    upper = nemdat$upper
  }else{
    upper=Inf
  }

  if(!is.null(nemdat$control)){
    control=nemdat$control
  }else{
    control=list(fnscale = -1)
  }

  # Finding good initial beta values
  init.betas <- as.numeric(lm.fit(ndmdat$X, ndmdat$y)$coefficients)

  # Finding posterior estimates using the normal approximation
  myopt <- optim(
    par = c(init.betas, 1, 0.36),
    fn = ndm_posterior,
    method = "BFGS",
    lower = lower,
    upper = upper,
    control = control,
    hessian = TRUE,
    y = ndmdat$y,
    X = ndmdat$X,
    A = ndmdat$network,
    bprior = bprior,
    s2prior = s2prior,
    rhoprior = rhoprior
  )

  beta_ests <- myopt$par[1:ndmdat$p]
  names(beta_ests) <- colnames(ndmdat$X)

  # Outputting the results
  list(
    beta = beta_ests,
    s2 = myopt$par[ndmdat$p + 1],
    rho = myopt$par[ndmdat$p + 2],
    hessian = myopt$hessian,
    loglik = ndm_likelihood(myopt$par, ndmdat$network, ndmdat$y, ndmdat$X),
    logpost = myopt$value,
    counts = myopt$counts,
    message = myopt$message,
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

  # 1. Sample from the Priors (The Proposal Distribution)
  # Beta
  if(is.list(ndmdat$beta_prior)){
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
    beta_samples <- mvtnorm::rmvnorm(samples, b_mean, b_sigma)
    bprior <- function(beta){ mvtnorm::dmvnorm(beta, b_mean, b_sigma) }
  }else{
    stop("Importance sampling doesn't support prior distribution changes.")
  }

  # S2 (Variance)
  if(is.list(ndmdat$s2_prior)){
    s2_samples <- invgamma::rinvgamma(samples, ndmdat$s2_prior$shape, ndmdat$s2_prior$scale)
    s2prior <- function(s2){ invgamma::dinvgamma(s2, ndmdat$s2_prior$shape, ndmdat$s2_prior$scale) }
  }else{
    stop("Importance sampling doesn't support prior distribution changes.")
  }

  # Rho
  if(is.list(ndmdat$rho_prior)){
    rho_samples <- truncnorm::rtruncnorm(samples, a = -1, b = 1, mean = ndmdat$rho_prior$mean, sd = ndmdat$rho_prior$sd)
    rhoprior <- function(rho){ truncnorm::dtruncnorm(rho, a=-1, b=1, mean=ndmdat$rho_prior$mean, sd=ndmdat$rho_prior$sd) }
  }else{
    stop("Importance sampling doesn't support prior distribution changes.")
  }

  # 2. Calculate Importance Weights (NDM Likelihood)
  log_weights <- numeric(samples)
  for(i in 1:samples) {
    theta <- c(beta_samples[i,], s2_samples[i], rho_samples[i])
    log_weights[i] <- ndm_likelihood(theta, ndmdat$network, ndmdat$y, ndmdat$X, log = TRUE)
  }

  # Normalize weights
  max_log_w <- max(log_weights)
  weights <- exp(log_weights - max_log_w)
  weights <- weights / sum(weights)

  # 3. Resample to get the "Posterior"
  indices <- sample(1:samples, size = samples, replace = TRUE, prob = weights)

  # Calculate means for output
  resampled_beta <- beta_samples[indices, ]
  resampled_s2   <- s2_samples[indices]
  resampled_rho  <- rho_samples[indices]

  theta_mean <- c(colMeans(resampled_beta), mean(resampled_s2), mean(resampled_rho))

  # Output
  list(
    beta = colMeans(resampled_beta),
    s2 = as.numeric(mean(resampled_s2)),
    rho = as.numeric(mean(resampled_rho)),
    samples = list(beta = resampled_beta, s2 = resampled_s2, rho = resampled_rho),
    logpost = ndm_posterior(theta_mean, ndmdat$network, ndmdat$y, ndmdat$X, bprior = bprior, s2prior = s2prior, rhoprior = rhoprior),
    loglik = ndm_likelihood(theta_mean, ndmdat$network, ndmdat$y, ndmdat$X, log=TRUE),
    log_marginal_lik = max_log_w + log(mean(exp(log_weights - max_log_w))),
    converged = TRUE
  )
}


# Code for when method="stan"
ndm_stan <- function(ndmdat){
  # Checking stan
  if(!requireNamespace("rstan", quietly = TRUE)){
    stop("Package 'rstan' is needed for this method. Please install it.", call. = FALSE)
  }

  # Checking for additional inputs
  iterations <- if(is.null(ndmdat$iterations)) 2000 else ndmdat$iterations
  chains <- if(is.null(ndmdat$chains)) 4 else ndmdat$chains

  # Ensure priors match the number of predictors (K)
  K <- ndmdat$p

  # 1. Expand b_mean if it's a single value
  b_mean_vec <- if(length(ndmdat$beta_prior$mean) == 1) {
    rep(ndmdat$beta_prior$mean, K)
  }else{
    as.vector(ndmdat$beta_prior$mean)
  }

  # 2. Expand b_sigma if it's a single value
  b_sigma_mat <- if(length(ndmdat$beta_prior$sigma) == 1) {
    diag(K) * ndmdat$beta_prior$sigma
  }else{
    as.matrix(ndmdat$beta_prior$sigma)
  }

  # Prepare data list for Stan
  stan_data <- list(
    N = length(ndmdat$y),
    K = ncol(ndmdat$X),
    y = ndmdat$y,
    X = ndmdat$X,
    A = ndmdat$network,
    b_mean = b_mean_vec,
    b_sigma = b_sigma_mat,
    s2_shape = ndmdat$s2_prior$shape,
    s2_scale = ndmdat$s2_prior$scale,
    rho_mean = ndmdat$rho_prior$mean,
    rho_sd = ndmdat$rho_prior$sd
  )

  # Fit the model
  fit <- rstan::sampling(
    stanmodels$ndm_disturbances,
    data = stan_data,
    iter = iterations,
    chains = chains
  )

  return(fit)
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

#' Summary method for ndm objects
#' @param object An object of secondary class "disturbances"
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
    s_mat <- rstan::summary(object$estimates, probs = probs)$summary
    beta_idx <- grep("^beta\\[", rownames(s_mat))
    s2_idx <- which(rownames(s_mat) == "s2")
    rho_idx <- which(rownames(s_mat) == "rho")

    coef_table <- s_mat[c(beta_idx, s2_idx, rho_idx), c("mean", prob_names)]
    rownames(coef_table)[1:p] <- colnames(object$data$X)

    logpost_val <- mean(rstan::extract(object$estimates, "lp__")$lp__)
    theta_star <- as.numeric(s_mat[c(beta_idx, s2_idx, rho_idx), "mean"])
    loglik_val <- ndm_likelihood(theta_star, object$data$network, object$data$y, object$data$X, log=TRUE)

  }else if(object$method %in% c("importance", "is")){
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
    # Normal/LNAM
    if(object$method == "normal"){
      theta_star <- c(object$estimates$beta, object$estimates$s2, object$estimates$rho)
      hinv <- tryCatch(solve(-object$estimates$hessian),
                       error = function(e) matrix(NA, length(theta_star), length(theta_star)))
      std_errs <- sqrt(diag(hinv))
      loglik_val <- object$estimates$loglik
      logpost_val <- object$estimates$logpost
    } else {
      # LNAM Specifics for Disturbances
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
    model = object$model,
    method = object$method,
    converged = if(object$method == "normal") object$estimates$converged else TRUE
  )

  class(res) <- "summary.nam"
  return(res)
}

#' @export
print.summary.nam <- function(x, ...) {
  cat("\nSummary of Network Analysis Model\n")
  cat("Model Type:", x$model, "| Estimation Method:", x$method, "\n")

  cat("\nCall:\n")
  print(x$call)

  cat("\nParameter Estimates and Intervals:\n")
  print(round(x$coefficients, 4))

  cat("\nModel Fit Statistics:\n")
  cat("Log-Likelihood: ", round(x$loglik, 4))

  if(!is.null(x$logpost) && !is.na(x$logpost)) {
    cat("\nLog-Posterior:  ", round(x$logpost, 4))
  }

  # Display Marginal Likelihood if available (from Importance Sampling)
  if(!is.null(x$marg_lik) && !is.na(x$marg_lik)) {
    cat("\nLog-Marginal-Lik:", round(x$marg_lik, 4))
  }

  if(!is.null(x$converged) && isFALSE(x$converged)) {
    cat("\n\n*** WARNING: Optimization did not converge! ***\n")
  }

  cat("\n")
  invisible(x)
}


