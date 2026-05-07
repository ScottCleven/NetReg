#' hand - A function to perform the Homophily Adjusted Network Disturbances Model
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
#' Options vary by model but all potential options are c("normal",
#' " ", "stan", "lnam"). Run ?methods.nam() for model-specific designations.
#' @param beta_prior The prior for the beta vector of the model. Requires list
#' or function inputs. Default is a multivariate normal with mean=0 and sigma=1
#' If you supply a list, you must have one element named 'mean' which is the
#' mean of the multivariate normal and one named 'sigma' which is the
#' variance-covariance matrix. If mean is a single integer,
#' that mean will be used for every element in the mean vector and if sigma is
#' a single integer, it will multiply that integer by the identity matrix.
#' Otherwise, both mean and sigma will be what you specify them to be. Or the
#' user can specify their own prior function with which to input. If
#' method=" ", the user-specified function must be a function from the
#'   package.
#' @param s2_prior Prior for the variance parameter. Requires list or function
#' inputs. Defaults to a gamma(2,0.5). If Input is a list, the list must have
#' two elements named 'shape' and 'scale' to represent the shape and scale of
#' the gamma distribution respectively. Or the user can specify their own prior
#' function with which to input. If method=" ", the user-specified function
#' must be a function from the   package.
#' @param rho_prior The prior for the correlation of the adjacency matrix of
#' the model. Requires list or function inputs. Default is a N(0.36,0.49)
#' (Dittrich et al. 2017). If input is a list, the list must have two elements
#' named 'mean' and 'sd' which are the mean and standard deviation of the normal
#' distribution respectively. Or the user can specify their own prior
#' function with which to input. If method=" ", the user-specified function
#' must be a function from the   package. If method="stan" only list input
#' will be accepted.
#' @param ... Additional arguments to be passed into the optim() function for
#' when method="normal".
#' @return A list with two classes, the first being "nam" and the second is the
#' class of the method that contains values:
#' @return \code{parameters} A list containing the parameters of the chosen method.
#' @return \code{A} The adjacency matrix. Will be row-normalized if rownorm=TRUE.
#' @return \code{loglik} The log of the likelihood of the model.
#' @return \code{AIC} The Akaike Information Criterion of the model.
#' @return \code{BIC} The Bayesian Information Criterion of the model.
#'
#' @import mvtnorm
#' @import rstan
#' @import Rcpp
#' @importFrom truncnorm dtruncnorm
#' @importFrom invgamma dinvgamma
#' @importFrom sna lnam
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @useDynLib nafm, .registration = TRUE
#' @export


# Model-specific wrapper function to call different methods
hand <- function(formula, data, network, method,
                beta_prior = list(mean=0, sigma = 1),
                s2_prior = list(shape=2, scale=0.5),
                rho_prior = list(mean=0.36, sd=0.7^2),
                ...){
  # Cleaning inputs
  X = model.matrix(formula, data=data)
  y = model.frame(formula, data)[,1]
  n = NROW(y)
  p = NCOL(X)

  if(tolower(method)=="stan" && typeof(beta_prior) != "list" |
     tolower(method)=="stan" && typeof(s2_prior) != "list" |
     tolower(method)=="stan" && typeof(rho_prior) != "list"){
    stop("stan method requires list input priors")
  }

  # Creating dataset input for individual methods
  handdat <- list(y=y, X=X, network=network, method=method, p=p, n=n, I=diag(n),
                 beta_prior=beta_prior,
                 s2_prior = s2_prior,
                 rho_prior = rho_prior,
                 ...)

  if(tolower(method)=="normal"){
    handout <- hand_normal(handdat)
  }else if(tolower(method)=="stan"){
    handout <- hand_stan(handdat)
  }else if(tolower(method)=="lnam"){
    handout <- hand_lnam(handdat)
  }
  handout
}

# Likelihood of the HAND model
hand_likelihood <- function(theta, A, y, X, Lambda, Omega, Psi, log = TRUE){
  p <- NCOL(X)
  D <- NCOL(Lambda)
  n <- NROW(y)

  beta <- theta[1:p]
  gamma <- theta[(p + 1):(p + D)]
  s2 <- theta[p + D + 1]
  rho <- theta[p + D + 2]

  # Mean: X*beta + Lambda*gamma (Simple structure)
  mu <- X %*% beta + Lambda %*% gamma

  # Covariance structure
  # Sigma = (gamma' * Psi * gamma) * Omega + M %*% (s2 * I) %*% M'
  M <- qr.solve(diag(n) - rho * A)
  gPg <- as.numeric(t(gamma) %*% Psi %*% gamma)

  V_hand <- gPg * Omega + tcrossprod(M * s2, M)

  mvtnorm::dmvnorm(y, mean = as.vector(mu), sigma = V_hand, log = log)
}


# Posterior of the HAND model
hand_posterior <- function(theta, A, y, X, Lambda, Omega, Psi,
                           bprior, gprior, s2prior, rhoprior, log = TRUE){
  ll <- hand_likelihood(theta, A, y, X, Lambda, Omega, Psi, log = TRUE)

  # Extract for priors
  p <- NCOL(X); D <- NCOL(Lambda)
  beta <- theta[1:p]; gamma <- theta[(p+1):(p+D)]
  s2 <- theta[p+D+1]; rho <- theta[p+D+2]

  lp <- ll + bprior(beta, log=T) + gprior(gamma, log=T) +
    s2prior(s2, log=T) + rhoprior(rho, log=T)

  return(if(log) lp else exp(lp))
}


# Code for when method="normal"
hand_normal <- function(handdat){
  # 1. Prior Setup
  if(is.list(handdat$beta_prior)){
    if(length(handdat$beta_prior$mean) == 1){
      b_mean <- rep(handdat$beta_prior$mean, handdat$p)
    }else{
      b_mean <- handdat$beta_prior$mean
    }
    b_sigma <- if(length(handdat$beta_prior$sigma) == 1){
      b_sigma <- diag(handdat$p) * handdat$beta_prior$sigma
    }else{
      b_sigma <- handdat$beta_prior$sigma
    }
    bprior <- function(beta, log=TRUE){
      mvtnorm::dmvnorm(beta, b_mean, b_sigma, log=log)
    }
  }else{
    bprior <- handdat$beta_prior
  }

  # Gamma Prior (Latent Space Coefficients)
  if(is.list(handdat$gamma_prior)){
    if(length(handdat$gamma_prior$mean) == 1){
      g_mean <- rep(handdat$gamma_prior$mean, handdat$D)
    }else{
      g_mean <- handdat$gamma_prior$mean
    }
    if(length(handdat$gamma_prior$sigma) == 1){
      g_sigma <- diag(handdat$D) * handdat$gamma_prior$sigma
    }else{
      g_sigma <- handdat$gamma_prior$sigma
    }
    gprior <- function(gamma, log=TRUE){
      mvtnorm::dmvnorm(gamma, g_mean, g_sigma, log=log)
    }
  }else{
    gprior <- handdat$gamma_prior
  }

  # S2 Prior
  if(is.list(handdat$s2_prior)){
    s2prior <- function(s2, log=TRUE){
      out <- invgamma::dinvgamma(s2, handdat$s2_prior$shape, handdat$s2_prior$scale)
      if(log){
        return(log(out))
      }else{
        return(out)
      }
    }
  }else{
    s2prior <- handdat$s2_prior
  }

  # Rho Prior
  if(is.list(handdat$rho_prior)){
    rhoprior <- function(rho, log=TRUE){
      truncnorm::dtruncnorm(rho, a=-1, b=1, mean=handdat$rho_prior$mean,
                            sd=handdat$rho_prior$sd)
    }
  }else{
    rhoprior <- handdat$rho_prior
  }

  # 2. Initialization (2SLS to get the optimizer in the ballpark)
  fit2sls <- lnam2sls_effect(handdat$y, cbind(handdat$X, handdat$Lambda),
                             handdat$A)
  # par structure: [beta, gamma, s2, rho]
  init_vec <- c(fit2sls$coefs[1:(handdat$p + handdat$D)], fit2sls$s2,
                fit2sls$coefs[handdat$p + handdat$D + 1])

  # 3. Optimization
  myopt <- optim(par = init_vec,
                 fn = hand_posterior,
                 method = "L-BFGS-B",
                 lower = c(rep(-Inf, handdat$p + handdat$D), 1e-5, -0.999),
                 upper = c(rep(Inf, handdat$p + handdat$D + 1), 0.999),
                 control = list(fnscale = -1),
                 hessian = TRUE,
                 y = handdat$y, X = handdat$X, A = handdat$A,
                 Lambda = handdat$Lambda, Omega = handdat$Omega, Psi = handdat$Psi,
                 bprior = bprior, gprior = gprior, s2prior = s2prior, rhoprior = rhoprior)

  list(beta = myopt$par[1:handdat$p],
       gamma = myopt$par[(handdat$p + 1):(handdat$p + handdat$D)],
       s2 = myopt$par[handdat$p + handdat$D + 1],
       rho = myopt$par[handdat$p + handdat$D + 2],
       hessian = myopt$hessian, loglik = myopt$value, converged = (myopt$convergence == 0))
}


# Code for when method="IS"
hand_IS <- function(handdat) {
  if(is.null(handdat$samples)) handdat$samples <- 10000
  S <- handdat$samples
  p <- handdat$p
  D <- handdat$D

  # 1. Sample from Priors (Proposal Distribution)
  # Beta
  if(is.list(handdat$beta_prior)) {
    if(length(handdat$beta_prior$mean) == 1){
      b_mu <- rep(handdat$beta_prior$mean, p)
    }else{
      b_mu <- handdat$beta_prior$mean
    }
    if(length(handdat$beta_prior$sigma) == 1){
      b_sig <- diag(p) * handdat$beta_prior$sigma
    }else{
      b_sig <- handdat$beta_prior$sigma
    }
    beta_samples <- mvtnorm::rmvnorm(S, b_mu, b_sig)
    bprior <- function(beta, log=TRUE){
      mvtnorm::dmvnorm(beta, b_mu, b_sig, log=log)
    }
  }else{
    stop("Importance sampling requires list-based priors for HAND.")
  }

  # Gamma
  if(is.list(handdat$gamma_prior)){
    if(length(handdat$gamma_prior$mean) == 1){
      g_mu <- rep(handdat$gamma_prior$mean, D)
    }else{
      g_mu <- handdat$gamma_prior$mean
    }
    if(length(handdat$gamma_prior$sigma) == 1){
      g_sig <- diag(D) * handdat$gamma_prior$sigma
    }else{
      g_sig <- handdat$gamma_prior$sigma
    }
    gamma_samples <- mvtnorm::rmvnorm(S, g_mu, g_sig)
    gprior <- function(gamma, log=TRUE){
      mvtnorm::dmvnorm(gamma, g_mu, g_sig, log=log)
    }
  }else{
    stop("Importance sampling requires list-based priors for HAND.")
  }

  # S2 & Rho
  s2_samples <- invgamma::rinvgamma(S, handdat$s2_prior$shape,
                                    handdat$s2_prior$scale)
  s2prior <- function(s2, log=TRUE){
    out <- invgamma::dinvgamma(s2, handdat$s2_prior$shape, handdat$s2_prior$scale)
    if(log){
      log(out)
    }else{
      out
    }
  }
  rho_samples <- truncnorm::rtruncnorm(S, a = -1, b = 1,
                                       mean = handdat$rho_prior$mean,
                                       sd = handdat$rho_prior$sd)
  rhoprior <- function(rho, log=TRUE){
    if(log){
      log(truncnorm::dtruncnorm(rho, a=-1, b=1,
                            mean=handdat$rho_prior$mean,
                            sd=handdat$rho_prior$sd))
    }else{
      truncnorm::dtruncnorm(rho, a=-1, b=1,
                            mean=handdat$rho_prior$mean,
                            sd=handdat$rho_prior$sd)
    }
  }

  # 2. Calculate Importance Weights (Log-Likelihood)
  log_weights <- numeric(S)
  for(i in 1:S) {
    theta <- c(beta_samples[i,], gamma_samples[i,], s2_samples[i],
               rho_samples[i])
    log_weights[i] <- hand_likelihood(theta, handdat$A, handdat$y, handdat$X,
                                      handdat$Lambda, handdat$Omega,
                                      handdat$Psi, log = TRUE)
  }

  # 3. Resample (Log-Sum-Exp Trick)
  max_log_w <- max(log_weights)
  weights <- exp(log_weights - max_log_w)
  weights <- weights / sum(weights)
  indices <- sample(1:S, size = S, replace = TRUE, prob = weights)

  # 4. Format Output
  theta_mean <- c(colMeans(beta_samples[indices,]),
                  colMeans(gamma_samples[indices,]),
                  mean(s2_samples[indices]), mean(rho_samples[indices]))

  list(
    beta = colMeans(beta_samples[indices,]),
    gamma = colMeans(gamma_samples[indices,]),
    s2 = mean(s2_samples[indices]),
    rho = mean(rho_samples[indices]),
    samples = list(beta = beta_samples[indices,],
                   gamma = gamma_samples[indices,],
                   s2 = s2_samples[indices], rho = rho_samples[indices]),
    logpost = hand_posterior(theta_mean, handdat$A, handdat$y, handdat$X,
                             handdat$Lambda, handdat$Omega, handdat$Psi,
                             bprior, gprior, s2prior, rhoprior),
    loglik = hand_likelihood(theta_mean, handdat$A, handdat$y, handdat$X,
                             handdat$Lambda, handdat$Omega, handdat$Psi,
                             log=TRUE),
    log_marginal_lik = max_log_w + log(mean(exp(log_weights - max_log_w))),
    converged = TRUE
  )
}



# Code for when method="stan"
hand_stan <- function(handdat){
  if(!requireNamespace("rstan", quietly = TRUE)){
    stop("rstan is required.")
  }

  # Set defaults for MCMC
  if(is.null(handdat$iterations)){
    iters <- 2000
  }else{
    iters <- handdat$iterations
  }
  if(is.null(handdat$chains)){
    chains <- 4
  }else{
    chains <- handdat$chains
  }

  # Prior Vectors
  if(length(handdat$beta_prior$mean) == 1){
    b_mu <- rep(handdat$beta_prior$mean, handdat$p)
  }else{
    b_mu <- as.vector(handdat$beta_prior$mean)
  }
  if(length(handdat$beta_prior$sigma) == 1){
    b_sig <- diag(handdat$p) * handdat$beta_prior$sigma
  }else{
    b_sig <- as.matrix(handdat$beta_prior$sigma)
  }

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
  stan_data <- list(
    N = handdat$n, K = handdat$p, D = handdat$D,
    y = handdat$y, X = handdat$X, A = handdat$A,
    Lambda = handdat$Lambda, Omega = handdat$Omega, Psi = handdat$Psi,
    b_mean = b_mu, b_sigma = b_sig,
    g_mean = g_mu, g_sigma = g_sig,
    s2_shape = handdat$s2_prior$shape, s2_scale = handdat$s2_prior$scale,
    rho_mean = handdat$rho_prior$mean, rho_sd = handdat$rho_prior$sd
  )

  # Assuming the HAND specific stan model is named 'hand_model'
  fit <- rstan::sampling(stanmodels$hand_model, data = stan_data, iter = iters, chains = chains)
  return(fit)
}


#' Summary method for HAND and HANE objects
#' @export
summary.hand <- function(object, probs = c(0.025, 0.975), ...){
  # HAND uses the same summary structure as HANE
  p <- ncol(object$data$X)
  D <- ncol(object$data$Lambda)
  prob_names <- paste0(probs * 100, "%")

  if(object$method == "stan"){
    s_mat <- rstan::summary(object$estimates, probs = probs)$summary

    # Identify indices
    beta_idx  <- grep("^beta\\[", rownames(s_mat))
    gamma_idx <- grep("^gamma\\[", rownames(s_mat))
    s2_idx    <- which(rownames(s_mat) == "s2")
    rho_idx   <- which(rownames(s_mat) == "rho")

    coef_table <- s_mat[c(beta_idx, gamma_idx, s2_idx, rho_idx),
                        c("mean", prob_names)]

    # Labeling
    rownames(coef_table)[1:p] <- colnames(object$data$X)
    rownames(coef_table)[(p+1):(p+D)] <- paste0("gamma_", 1:D)

    logpost_val <- mean(rstan::extract(object$estimates, "lp__")$lp__)
    theta_star <- as.numeric(s_mat[c(beta_idx, gamma_idx, s2_idx, rho_idx),
                                   "mean"])

    # Dynamic likelihood call based on model type
    if(object$model == "hand"){
      loglik_val <- hand_likelihood(theta_star, object$data$network,
                                    object$data$y,
                                    object$data$X, object$data$Lambda,
                                    object$data$Omega, object$data$Psi)
    } else {
      loglik_val <- hane_likelihood(theta_star, object$data$network,
                                    object$data$y,
                                    object$data$X, object$data$Lambda,
                                    object$data$Omega, object$data$Psi)
    }

  } else if(object$method %in% c("importance", "is")){
    samples <- object$estimates$samples
    theta_samples <- cbind(samples$beta, samples$gamma, sigma2 = samples$s2,
                           rho = samples$rho)

    coef_table <- t(apply(theta_samples, 2, function(x){
      c(Estimate = mean(x), quantile(x, probs = probs))
    }))

    rownames(coef_table)[1:p] <- colnames(object$data$X)
    rownames(coef_table)[(p+1):(p+D)] <- paste0("gamma_", 1:D)

    loglik_val <- object$estimates$loglik
    logpost_val <- object$estimates$logpost
    marg_lik <- object$estimates$log_marginal_lik

  } else {
    # Normal Approximation
    theta_star <- c(object$estimates$beta, object$estimates$gamma,
                    object$estimates$s2, object$estimates$rho)

    hinv <- tryCatch(solve(-object$estimates$hessian),
                     error = function(e) matrix(NA, length(theta_star),
                                                length(theta_star)))

    std_errs <- sqrt(pmax(0, diag(hinv)))
    z_crit <- qnorm(probs)
    intervals <- matrix(NA, nrow = length(theta_star), ncol = length(probs))
    for(i in 1:length(probs)) {
      intervals[,i] <- theta_star + z_crit[i] * std_errs
    }

    coef_table <- cbind(Estimate = theta_star, intervals)
    colnames(coef_table) <- c("Estimate", prob_names)
    rownames(coef_table) <- c(colnames(object$data$X), paste0("gamma_", 1:D),
                              "sigma2", "rho")

    loglik_val <- object$estimates$loglik
    logpost_val <- object$estimates$logpost
  }

  res <- list(
    call = object$call,
    coefficients = coef_table,
    loglik = loglik_val,
    logpost = logpost_val,
    model = ifelse(object$model == "hand", "HAND (Disturbances)", "HANE (Effects)"),
    method = object$method,
    converged = if(object$method == "normal") object$estimates$converged else TRUE
  )
  class(res) <- "summary.nam"
  return(res)
}
