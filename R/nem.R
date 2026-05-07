#' nem - A function to perform the Network Effects Model
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
#' class of "effects" that contains values:
#' @return \code{estimates} The raw estimate output of the chosen method.
#' @return \code{formula} The inputted formula
#' @return \code{model} The chosen model ("effects").
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


# Model-specific wrapper function to call different methods
nem <- function(formula, data=NULL, network, method,
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
    network=row_norm(network)
  }


  # Check stan prior inputs are list only.
  if((tolower(method)=="stan" && typeof(beta_prior) != "list") |
     (tolower(method)=="stan" && typeof(s2_prior) != "list") |
     (tolower(method)=="stan" && typeof(rho_prior) != "list")){
    stop("stan method requires list input of prior parameters,
         prior functions cannot be changed.")
  }


  nemdat <- list(y=y, X=X, network=network, method=method, p=p, n=n, I=diag(n),
                 beta_prior=beta_prior,
                 s2_prior = s2_prior,
                 rho_prior = rho_prior,
                 ...)

  if(tolower(method)=="normal"){
    nemout <- nem_normal(nemdat)
  }else if(tolower(method)=="stan"){
    nemout <- nem_stan(nemdat)
  }else if(tolower(method)=="lnam"){
    nemout <- nem_lnam(nemdat)
  }else if((tolower(method)=="is") | (tolower(method)=="importance sampling")){
    nemout <- nem_IS(nemdat)
  }
  final_out <- list(
    estimates = nemout,        # The raw results from the specific method
    formula   = formula,       # The symbolic formula used
    model     = "effects",     # String identifying the model type
    method    = tolower(method),    # String identifying the estimation method
    data      = list(y = y, X = X, network = network), # Original data for summaries/plots
    priors    = list(beta = beta_prior, s2 = s2_prior, rho = rho_prior), # Prior settings
    call      = match.call()
  )

  # Assign the S3 classes
  class(final_out) <- c("nam", "effects")

  return(final_out)

}

# Likelihood of the effects model.
nem_likelihood <- function(theta, A, y, X, log=TRUE){
  beta <- theta[1:(NROW(theta)-2)]
  s2 <- theta[NROW(theta)-1]
  rho <- theta[NROW(theta)]

  n <- NROW(y)
  p <- length(beta)
  I <- diag(n)

  Launtif <- qr.solve(I - rho*A)

  mvtnorm::dmvnorm(y, Launtif %*% X %*% beta, s2 * tcrossprod(Launtif), log = log)

}

# Posterior of the effects model.
nem_posterior <- function(theta, A, y, X, log=TRUE,
                          bprior, s2prior, rhoprior){
  beta <- theta[1:(NROW(theta)-2)]
  s2 <- theta[NROW(theta)-1]
  rho <- theta[NROW(theta)]

  n <- NROW(y)
  p <- length(beta)
  I <- diag(n)

  Launtif <- qr.solve(I - rho*A)

  out <- nem_likelihood(theta, A, y, X, log=FALSE) *
    bprior(beta) * s2prior(s2) * rhoprior(rho)
  if(log){
    out <- log(out)
  }
  out
}

# Code for when method="normal"
nem_normal <- function(nemdat){

  # Checking Prior Inputs
  if(typeof(nemdat$beta_prior) == "list"){
    if(length(nemdat$beta_prior$mean) == 1){
      bmean <- rep(nemdat$beta_prior$mean,nemdat$p)
    }else{
      bmean <- nemdat$beta_prior$mean
    }

    if(length(nemdat$beta_prior$sigma) == 1){
      bsigma <- diag(nemdat$p) * nemdat$beta_prior$sigma
    }else{
      bsigma <- nemdat$beta_prior$sigma
    }

    bprior <- function(beta){
      mvtnorm::dmvnorm(beta, bmean, bsigma)
    }
  }else{
    bprior <- nemdat$beta_prior
  }

  if(typeof(nemdat$s2_prior) == "list"){
    s2prior <- function(s2){
      invgamma::dinvgamma(s2, nemdat$s2_prior$shape, nemdat$s2_prior$scale)
    }
  }else{
    s2prior <- nemdat$s2_prior
  }

  if(typeof(nemdat$rho_prior) == "list"){
    rhoprior <- function(rho){
      truncnorm::dtruncnorm(rho, a=-1, b=1, mean=nemdat$rho_prior$mean,
                            sd=nemdat$rho_prior$sd)
    }
  }else{
    rhoprior <- nemdat$rho_prior
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
  init.betas <- lm(nemdat$y ~ nemdat$X[,-1])$coefficients

  # Finding posterior estimates using the normal approximation of the posterior
  myopt <- optim(par=c(init.betas, 1, 0.36), fn=nem_posterior,
        method = "BFGS",
        lower=lower,
        upper=upper,
        control=control,
        hessian=TRUE,
        y=nemdat$y,
        X=nemdat$X,
        A=nemdat$network,
        bprior=bprior,
        s2prior=s2prior,
        rhoprior=rhoprior
        )

  beta_ests <- myopt$par[1:nemdat$p]
  names(beta_ests) <- colnames(nemdat$X)
  # Outputting the results
  list(
    beta     = beta_ests,
    s2       = myopt$par[nemdat$p + 1],
    rho      = myopt$par[nemdat$p + 2],
    hessian  = myopt$hessian,
    loglik   = nem_likelihood(myopt$par, nemdat$network, nemdat$y, nemdat$X),
    logpost  = myopt$value,
    counts   = myopt$counts,
    message  = myopt$message,
    converged = myopt$convergence == 0
  )
}

# Code for when method="IS"
nem_IS <- function(nemdat){
  if(is.null(nemdat$samples)){
    nemdat$samples = 10000
  }
  samples <- nemdat$samples
  p <- nemdat$p

  # 1. Sample from the Priors (The Proposal Distribution)
  # Beta
  if(is.list(nemdat$beta_prior)){
    if(length(nemdat$beta_prior$mean) == 1){
      b_mean <- rep(nemdat$beta_prior$mean, p)
    }else{
      b_mean <- nemdat$beta_prior$mean
    }

    if(length(nemdat$beta_prior$sigma) == 1){
      b_sigma <- diag(p) * nemdat$beta_prior$sigma
    }else{
      b_sigma <- nemdat$beta_prior$sigma
    }
    beta_samples <- mvtnorm::rmvnorm(samples, b_mean, b_sigma)
    bprior <- function(beta){
      mvtnorm::dmvnorm(beta, b_mean, b_sigma)
    }
  }else{
    stop("Importance sampling doesn't support prior distribution changes.")
  }

  # S2 (Variance) - Sampling from Inverse Gamma
  if(is.list(nemdat$s2_prior)){
    s2_samples <- invgamma::rinvgamma(samples, nemdat$s2_prior$shape,
                                      nemdat$s2_prior$scale)
    s2prior <- function(s2){
      invgamma::dinvgamma(s2, nemdat$s2_prior$shape, nemdat$s2_prior$scale)
    }
  }else{
    stop("Importance sampling doesn't support prior distribution changes.")
  }

  # Rho - Sampling from Truncated Normal
  if(is.list(nemdat$rho_prior)){
    rho_samples <- truncnorm::rtruncnorm(samples, a = -1, b = 1,
                                       mean = nemdat$rho_prior$mean,
                                       sd = nemdat$rho_prior$sd)
    rhoprior <- function(rho){
      truncnorm::dtruncnorm(rho, a=-1, b=1, mean=nemdat$rho_prior$mean,
                            sd=nemdat$rho_prior$sd)
    }
  }else{
    stop("Importance sampling doesn't support prior distribution changes.")
  }

  # 2. Calculate Importance Weights (Likelihood)
  log_weights <- numeric(samples)
  for(i in 1:samples) {
    theta <- c(beta_samples[i,], s2_samples[i], rho_samples[i])
    log_weights[i] <- nem_likelihood(theta, nemdat$network, nemdat$y, nemdat$X, log = TRUE)
  }

  # Normalize weights using the log-sum-exp trick for stability
  max_log_w <- max(log_weights)
  weights <- exp(log_weights - max_log_w)
  weights <- weights / sum(weights)

  # 3. Resample to get the "Posterior"
  indices <- sample(1:samples, size = samples, replace = TRUE, prob = weights)

  theta_mean <- c(colMeans(beta_samples[indices, ]),
                  mean(s2_samples[indices]),
                  mean(rho_samples[indices]))

  # Output formatted for your existing print/summary methods
  list(
    beta = colMeans(beta_samples[indices, ]),
    s2 = as.numeric(mean(s2_samples[indices])),
    rho = as.numeric(mean(rho_samples[indices])),
    samples = list(beta = beta_samples[indices, ], s2 = s2_samples[indices], rho = rho_samples[indices]),
    logpost = nem_posterior(theta_mean, nemdat$network, nemdat$y, nemdat$X,
                            bprior = bprior, s2prior = s2prior,
                            rhoprior = rhoprior),
    loglik = nem_likelihood(theta_mean, nemdat$network, nemdat$y, nemdat$X,
                            log=TRUE),
    log_marginal_lik  = max_log_w + log(mean(exp(log_weights - max_log_w))), # Marginal Likelihood estimate
    converged = TRUE
  )
}


# Code for when method="stan"
nem_stan <- function(nemdat){

  #Checking stan
  if(!requireNamespace("rstan", quietly = TRUE)){
    stop("Package 'rstan' is needed for this method. Please install it.", call. = FALSE)
  }

  #Checking for additional inputs
  if(is.null(nemdat$iterations)){
    iterations=2000
  }else{
      iterations=nemdat$iterations
    }
  if(is.null(nemdat$chains)){
    chains=4
  }else{
    chains=nemdat$chains
  }

  # Ensure priors match the number of predictors (K)
  K <- nemdat$p

  # 1. Expand b_mean if it's a single value
  b_mean_vec <- if(length(nemdat$beta_prior$mean) == 1) {
    rep(nemdat$beta_prior$mean, K)
  } else {
    as.vector(nemdat$beta_prior$mean)
  }

  # 2. Expand b_sigma if it's a single value (Identity * sigma)
  b_sigma_mat <- if(length(nemdat$beta_prior$sigma) == 1) {
    diag(K) * nemdat$beta_prior$sigma
  } else {
    as.matrix(nemdat$beta_prior$sigma)
  }

  # Prepare data list for Stan
  stan_data <- list(
    N = length(nemdat$y),
    K = ncol(nemdat$X),
    y = nemdat$y,
    X = nemdat$X,
    A = nemdat$network,
    # Mapping the list priors to Stan inputs
    b_mean  = b_mean_vec,
    b_sigma = b_sigma_mat,
    s2_shape = nemdat$s2_prior$shape,
    s2_scale = nemdat$s2_prior$scale,
    rho_mean = nemdat$rho_prior$mean,
    rho_sd   = nemdat$rho_prior$sd
  )

  # Fit the model
  fit <- rstan::sampling(
    stanmodels$nem_effects,
    data = stan_data,
    iter = iterations,
    chains = chains
  )

  return(fit)
}

# Code for when method="lnam"
nem_lnam <- function(nemdat){

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

  # Running lnam
  res <- sna::lnam(y=nemdat$y, x=nemdat$X, W1=nemdat$network, theta.seed=theta.seed,
            null.model=null.model, method="BFGS", control=control, tol=tol)

  list(
    beta = res$beta,
    s2   = res$sigma.sq,
    rho  = res$rho1,
    raw  = res
  )
}






#' Summary method for nem objects
#' @param object An object of secondary class "effects"
#' @param probs A numeric vector of probabilities for the interval estimates.
#' Default is c(0.025, 0.975).
#' @param ... Additional arguments (ignored)
#' @export
summary.effects <- function(object, probs = c(0.025, 0.975), ...){
  p <- ncol(object$data$X)
  loglik_val <- NA
  logpost_val <- NA
  marg_lik <- NA

  # Ensure probs are named for the table columns
  prob_names <- paste0(probs * 100, "%")

  if(object$method == "stan"){

    s_mat <- rstan::summary(object$estimates, probs = probs)$summary

    beta_idx <- grep("^beta\\[", rownames(s_mat))
    s2_idx   <- which(rownames(s_mat) == "s2")
    rho_idx  <- which(rownames(s_mat) == "rho")

    # Extract mean and the specific quantile columns requested
    coef_table <- s_mat[c(beta_idx, s2_idx, rho_idx), c("mean", prob_names)]
    rownames(coef_table)[1:p] <- colnames(object$data$X)

    logpost_val <- mean(rstan::extract(object$estimates, "lp__")$lp__)
    theta_star  <- as.numeric(s_mat[c(beta_idx, s2_idx, rho_idx), "mean"])
    loglik_val  <- nem_likelihood(theta_star, object$data$network, object$data$y, object$data$X, log=TRUE)

  }else if(object$method %in% c("importance", "is")){
    samples <- object$estimates$samples
    theta_samples <- cbind(samples$beta, sigma2 = samples$s2, rho = samples$rho)

    # Calculate quantiles from the resampled distribution
    coef_table <- t(apply(theta_samples, 2, function(x){
      c(Estimate = mean(x), quantile(x, probs = probs))
    }))

    rownames(coef_table)[1:p] <- colnames(object$data$X)

    # Pull values directly from nem_IS output
    loglik_val  <- object$estimates$loglik
    logpost_val <- object$estimates$logpost
    marg_lik    <- object$estimates$log_marginal_lik

  }else{
    # 2. Normal/LNAM: Use Normal Approximation (Estimate +/- z * SE)
    if(object$method == "normal"){
      theta_star <- c(object$estimates$beta, object$estimates$s2, object$estimates$rho)
      hinv <- tryCatch(solve(-object$estimates$hessian),
                       error = function(e) matrix(NA, length(theta_star), length(theta_star)))
      std_errs <- sqrt(diag(hinv))
      loglik_val <- object$estimates$loglik
      logpost_val <- object$estimates$logpost
    } else {
      est <- object$estimates$raw
      theta_star <- c(est$beta, est$sigma.sq, est$rho1)
      std_errs <- c(est$se.beta, est$se.sigma.sq, est$se.rho1)
      loglik_val <- est$loglik
    }

    # Calculate intervals based on Normal Distribution
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






