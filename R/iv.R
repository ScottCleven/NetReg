#' IV - A function to perform the Kelejian and Prucha (SARAR) Model
#' @param formula symbolic description of the model to be fitted.
#' @param data optional data frame or environment.
#' @param network Adjacency matrix.
#' @param method Estimation method: c("normal", "IS", "stan").
#' @param beta_prior List or function for beta vector.
#' @param s2_prior List or function for variance.
#' @param rho_prior List or function for the outcome lag correlation.
#' @param lambda_prior List or function for the error lag correlation.
#' @param rownorm Default to TRUE.
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

IV <- function(formula, data=NULL, network, method,
               beta_prior = list(mean=0, sigma = 1),
               s2_prior = list(shape=2, scale=0.5),
               rho_prior = list(mean=0.36, sd=0.7^2),
               lambda_prior = list(mean=0.36, sd=0.7^2),
               rownorm=TRUE, ...){

  X = model.matrix(formula, data=data)
  y = model.frame(formula, data)[,1]
  n = NROW(y)
  p = NCOL(X)

  network = convert_adjacency(network, "matrix")$A
  if(rownorm){
    network = row_norm(network)
  }

  if(tolower(method) == "stan" | tolower(method) == "is"){
    check_list <- list(beta_prior, s2_prior, rho_prior, lambda_prior)
    if(any(sapply(check_list, function(x) typeof(x) != "list"))){
      stop("stan method requires list input of prior parameters;
           custom prior functions are not supported for this method.")
    }
  }

  ivdat <- list(
    y = y, X = X, network = network, method = tolower(method),
    p = p, n = n, I = diag(n),
    beta_prior = beta_prior, s2_prior = s2_prior,
    rho_prior = rho_prior, lambda_prior = lambda_prior,
    ...
  )

  if(tolower(method) == "normal"){
    ivout <- iv_normal(ivdat)
  } else if(tolower(method) == "stan"){
    ivout <- iv_stan(ivdat)
  } else if(tolower(method) %in% c("is", "importance sampling")){
    ivout <- iv_IS(ivdat)
  } else {
    stop("Unknown method. Please choose 'normal', 'stan', or 'IS'.")
  }


  final_out <- list(
    estimates = ivout,
    formula = formula,
    model = "iv",
    method = tolower(method),
    data = list(y = y, X = X, network = network),
    priors = list(beta = beta_prior, s2 = s2_prior, rho = rho_prior,
                  lambda = lambda_prior),
    call = match.call()
  )

  class(final_out) <- c("nam", "iv")
  return(final_out)
}




# Likelihood of the IV (SARAR 1,1) model
iv_likelihood <- function(theta, A, y, X, log = TRUE) {
  p <- NCOL(X)
  n <- NROW(y)

  # Extract parameters
  beta <- theta[1:p]
  s2 <- theta[p + 1]
  rho <- theta[p + 2]
  lambda <- theta[p + 3]

  I <- diag(n)

  # 1. Structural Matrices
  W_rho <- I - rho * A
  W_lambda <- I - lambda * A

  # 2. Calculate the 'innovations' (epsilon)
  epsilon <- W_lambda %*% (W_rho %*% y - X %*% beta)

  # 3. Log-Determinants (The Jacobian adjustments)
  log_det_rho <- as.numeric(determinant(W_rho, logarithm = TRUE)$modulus)
  log_det_lambda <- as.numeric(determinant(W_lambda, logarithm = TRUE)$modulus)

  # 4. Calculate the Log-Likelihood
  ll <- log_det_rho + log_det_lambda - (n/2) * log(2 * pi * s2) -
    (1 / (2 * s2)) * sum(epsilon^2)

  if (log) {
    return(ll)
  } else {
    return(exp(ll))
  }
}



# Posterior of the IV (SARAR 1,1) model
iv_posterior <- function(theta, A, y, X, log=TRUE, bprior, s2prior, rhoprior,
                         lprior){
  p <- NCOL(X)

  # Extract parameters for prior functions
  beta <- theta[1:p]
  s2 <- theta[p + 1]
  rho <- theta[p + 2]
  lambda <- theta[p + 3]

  # 1. Calculate Log-Likelihood
  ll <- iv_likelihood(theta, A, y, X, log = TRUE)

  # 2. Calculate Log-Priors
  lp_beta   <- bprior(beta, log = TRUE)
  lp_s2     <- s2prior(s2, log = TRUE)
  lp_rho    <- rhoprior(rho, log = TRUE)
  lp_lambda <- lprior(lambda, log = TRUE)

  # 3. Aggregate
  out <- ll + lp_beta + lp_s2 + lp_rho + lp_lambda

  if(log){
    return(out)
  } else {
    return(exp(out))
  }
}


# Code for when method="normal" for IV model
iv_normal <- function(ivdat){

  # 1. Prior Setup (List-to-Function Conversion)
  # Beta
  if(is.list(ivdat$beta_prior)){
    if(length(ivdat$beta_prior$mean) == 1){
      b_mean <- rep(ivdat$beta_prior$mean, ivdat$p)
    }else{
      b_mean <- ivdat$beta_prior$mean
    }
    if(length(ivdat$beta_prior$sigma) == 1){
      b_sigma <- diag(ivdat$p) * ivdat$beta_prior$sigma
    }else{
      b_sigma <- ivdat$beta_prior$sigma
    }
    bprior <- function(beta, log=TRUE){
      mvtnorm::dmvnorm(beta, b_mean, b_sigma, log=log)
      }
  }else{
    bprior <- ivdat$beta_prior
    }

  # S2
  if(is.list(ivdat$s2_prior)){
    s2prior <- function(s2, log=TRUE){
      out <- invgamma::dinvgamma(s2, ivdat$s2_prior$shape, ivdat$s2_prior$scale)
      if(log){
        return(log(out))
      }else{
        return(out)
      }
    }
  }else{
    s2prior <- ivdat$s2_prior
  }

  # Rho (Outcome Lag)
  if(is.list(ivdat$rho_prior)){
    rhoprior <- function(rho, log=TRUE){
      if(log){
        log(truncnorm::dtruncnorm(rho, a=-1, b=1, mean=ivdat$rho_prior$mean,
                                  sd=ivdat$rho_prior$sd))
      }else{
        truncnorm::dtruncnorm(rho, a=-1, b=1, mean=ivdat$rho_prior$mean,
                              sd=ivdat$rho_prior$sd)
      }
    }
  }else{
    rhoprior <- ivdat$rho_prior
  }

  # Lambda (Error Lag)
  if(is.list(ivdat$lambda_prior)){
    lprior <- function(lambda, log=TRUE){
      if(log){
        log(truncnorm::dtruncnorm(lambda, a=-1, b=1,
                                  mean=ivdat$lambda_prior$mean,
                                  sd=ivdat$lambda_prior$sd))
      }else{
        truncnorm::dtruncnorm(lambda, a=-1, b=1, mean=ivdat$lambda_prior$mean,
                              sd=ivdat$lambda_prior$sd)
      }
    }
  }else{
    lprior <- ivdat$lambda_prior
  }

  # 2. Initialization using IV (2SLS) Logic
  # Instruments: X, AX, A^2X
  AX <- ivdat$network %*% ivdat$X
  AAX <- ivdat$network %*% AX
  H <- cbind(ivdat$X, AX, AAX)

  # Stage 1: Project Ay onto instruments H
  Ay <- ivdat$network %*% ivdat$y
  Ay_hat <- H %*% solve(crossprod(H), crossprod(H, Ay))

  # Stage 2: Estimate Rho and Beta
  X_combined <- cbind(Ay_hat, ivdat$X)
  fit_2sls <- solve(crossprod(X_combined), crossprod(X_combined, ivdat$y))

  # Initial Theta: [Beta, S2, Rho, Lambda]
  init_vec <- c(fit_2sls[-1], var(ivdat$y), fit_2sls[1], 0.1)

  # 3. Optimization
  myopt <- optim(par = init_vec,
                 fn = iv_posterior,
                 method = "L-BFGS-B",
                 lower = c(rep(-Inf, ivdat$p), 1e-5, -0.999, -0.999),
                 upper = c(rep(Inf, ivdat$p + 1), 0.999, 0.999),
                 control = list(fnscale = -1),
                 hessian = TRUE,
                 y = ivdat$y, X = ivdat$X, A = ivdat$network,
                 bprior = bprior, s2prior = s2prior, rhoprior = rhoprior,
                 lprior = lprior)

  # 4. Return formatted results
  list(
    beta = myopt$par[1:ivdat$p],
    s2 = myopt$par[ivdat$p + 1],
    rho = myopt$par[ivdat$p + 2],
    lambda = myopt$par[ivdat$p + 3],
    hessian = myopt$hessian,
    loglik = iv_likelihood(myopt$par, ivdat$network, ivdat$y, ivdat$X),
    logpost = myopt$value,
    converged = (myopt$convergence == 0)
  )
}


# Code for IV method="IS"
iv_IS <- function(ivdat) {
  if(is.null(ivdat$samples)) ivdat$samples <- 10000
  S <- ivdat$samples
  p <- ivdat$p

  # 1. Sample from Priors (Proposal Distribution)

  # Beta
  if(is.list(ivdat$beta_prior)){
    if(length(ivdat$beta_prior$mean) == 1){
      b_mu <- rep(ivdat$beta_prior$mean, p)
    }else{
      b_mu <- ivdat$beta_prior$mean
    }
    if(length(ivdat$beta_prior$sigma) == 1){
      b_sig <- diag(p) * ivdat$beta_prior$sigma
    }else{
      b_sig <- ivdat$beta_prior$sigma
    }
    beta_samples <- mvtnorm::rmvnorm(S, b_mu, b_sig)
    bprior <- function(beta, log=TRUE){
      mvtnorm::dmvnorm(beta, b_mu, b_sig, log=log)
    }
  }else{
    stop("Importance sampling requires list-based priors.")
  }

  # S2 (Variance)
  s2_samples <- invgamma::rinvgamma(S, ivdat$s2_prior$shape,
                                    ivdat$s2_prior$scale)
  s2prior <- function(s2, log=TRUE){
    out <- invgamma::dinvgamma(s2, ivdat$s2_prior$shape, ivdat$s2_prior$scale)
    if(log){
      log(out)
    }else{
      out
    }
  }

  # Rho (Outcome Lag)
  rho_samples <- truncnorm::rtruncnorm(S, a = -1, b = 1,
                                       mean = ivdat$rho_prior$mean,
                                       sd = ivdat$rho_prior$sd)
  rhoprior <- function(rho, log=TRUE){
    if(log){
      log(truncnorm::dtruncnorm(rho, a=-1, b=1, mean=ivdat$rho_prior$mean,
                                sd=ivdat$rho_prior$sd))
    }else{
      truncnorm::dtruncnorm(rho, a=-1, b=1, mean=ivdat$rho_prior$mean,
                            sd=ivdat$rho_prior$sd)
    }
  }

  # Lambda (Error Lag)
  lambda_samples <- truncnorm::rtruncnorm(S, a = -1, b = 1,
                                          mean = ivdat$lambda_prior$mean,
                                          sd = ivdat$lambda_prior$sd)
  lprior <- function(lambda, log=TRUE){
    if(log){
      log(truncnorm::dtruncnorm(lambda, a=-1, b=1, mean=ivdat$lambda_prior$mean,
                                sd=ivdat$lambda_prior$sd))
    }else{
      truncnorm::dtruncnorm(lambda, a=-1, b=1, mean=ivdat$lambda_prior$mean,
                            sd=ivdat$lambda_prior$sd)
    }
  }
  # 2. Calculate Importance Weights (Log-Likelihood)
  log_weights <- numeric(S)
  for(i in 1:S){
    # theta structure: [beta, s2, rho, lambda]
    theta <- c(beta_samples[i,], s2_samples[i], rho_samples[i],
               lambda_samples[i])
    log_weights[i] <- iv_likelihood(theta, ivdat$network, ivdat$y, ivdat$X,
                                    log = TRUE)
  }

  # 3. Resample using Log-Sum-Exp Trick for stability
  max_log_w <- max(log_weights)
  weights <- exp(log_weights - max_log_w)
  weights <- weights / sum(weights)

  indices <- sample(1:S, size = S, replace = TRUE, prob = weights)

  # 4. Extract resampled "Posterior" parameters
  resampled_beta <- beta_samples[indices, ]
  resampled_s2 <- s2_samples[indices]
  resampled_rho <- rho_samples[indices]
  resampled_lambda <- lambda_samples[indices]

  theta_mean <- c(colMeans(resampled_beta), mean(resampled_s2),
                  mean(resampled_rho), mean(resampled_lambda))

  # 5. Format Output
  list(
    beta = colMeans(resampled_beta),
    s2 = mean(resampled_s2),
    rho = mean(resampled_rho),
    lambda = mean(resampled_lambda),
    samples = list(beta = resampled_beta, s2 = resampled_s2,
                   rho = resampled_rho, lambda = resampled_lambda),
    logpost = iv_posterior(theta_mean, ivdat$network, ivdat$y, ivdat$X,
                           bprior=bprior, s2prior=s2prior, rhoprior=rhoprior,
                           lprior=lprior),
    loglik = iv_likelihood(theta_mean, ivdat$network, ivdat$y, ivdat$X,
                           log=TRUE),
    log_marginal_lik = max_log_w + log(mean(exp(log_weights - max_log_w))),
    converged = TRUE
  )
}



# Code for when method="stan" for IV (SARAR 1,1)
iv_stan <- function(ivdat) {
  if(!requireNamespace("rstan", quietly = TRUE)){
    stop("rstan is required.")
  }

  # Set defaults for MCMC
  if(is.null(hanedat$iterations)){
    iters <- 2000
  }else{
    iters <- hanedat$iterations
  }
  if(is.null(hanedat$chains)){
    chains <- 4
  }else{
    chains <- hanedat$chains
  }

  # Expand Priors
  if(length(ivdat$beta_prior$mean) == 1){
    b_mu <- rep(ivdat$beta_prior$mean, ivdat$p)
  }else{
    b_mu <- as.vector(ivdat$beta_prior$mean)
  }
  if(length(ivdat$beta_prior$sigma) == 1){
    b_sig <- diag(ivdat$p) * ivdat$beta_prior$sigma
  }else{
    b_sig <- as.matrix(ivdat$beta_prior$sigma)
  }

  stan_data <- list(
    N = ivdat$n, K = ivdat$p,
    y = ivdat$y, X = ivdat$X, A = ivdat$network,
    b_mean = b_mu, b_sigma = b_sig,
    s2_shape = ivdat$s2_prior$shape, s2_scale = ivdat$s2_prior$scale,
    rho_mean = ivdat$rho_prior$mean, rho_sd = ivdat$rho_prior$sd,
    lambda_mean = ivdat$lambda_prior$mean, lambda_sd = ivdat$lambda_prior$sd
  )

  # Call the pre-compiled SARAR model
  fit <- rstan::sampling(stanmodels$iv_model, data = stan_data, iter = iters,
                         chains = chains)
  return(fit)
}


#' Summary method for IV objects
#' @export
summary.iv <- function(object, probs = c(0.025, 0.975), ...){
  p <- ncol(object$data$X)
  prob_names <- paste0(probs * 100, "%")

  if(object$method == "stan"){
    s_mat <- rstan::summary(object$estimates, probs = probs)$summary

    # Identify indices for IV parameters
    beta_idx   <- grep("^beta\\[", rownames(s_mat))
    s2_idx     <- which(rownames(s_mat) == "s2")
    rho_idx    <- which(rownames(s_mat) == "rho")
    lambda_idx <- which(rownames(s_mat) == "lambda")

    coef_table <- s_mat[c(beta_idx, s2_idx, rho_idx, lambda_idx),
                        c("mean", prob_names)]

    # Labeling
    rownames(coef_table)[1:p] <- colnames(object$data$X)

    logpost_val <- mean(rstan::extract(object$estimates, "lp__")$lp__)
    theta_star <- as.numeric(s_mat[c(beta_idx, s2_idx, rho_idx, lambda_idx),
                                   "mean"])
    loglik_val <- iv_likelihood(theta_star, object$data$network, object$data$y,
                                object$data$X)

  }else if(object$method %in% c("importance", "is")){
    samples <- object$estimates$samples
    # Combine: Beta, S2, Rho, Lambda
    theta_samples <- cbind(samples$beta, sigma2 = samples$s2, rho = samples$rho,
                           lambda = samples$lambda)

    coef_table <- t(apply(theta_samples, 2, function(x){
      c(Estimate = mean(x), quantile(x, probs = probs))
    }))

    rownames(coef_table)[1:p] <- colnames(object$data$X)

    loglik_val <- object$estimates$loglik
    logpost_val <- object$estimates$logpost
    marg_lik <- object$estimates$log_marginal_lik

  }else{
    # Normal Approximation
    theta_star <- c(object$estimates$beta, object$estimates$s2,
                    object$estimates$rho, object$estimates$lambda)

    # Robust Hessian inversion
    hinv <- tryCatch(solve(-object$estimates$hessian),
                     error = function(e) matrix(NA, length(theta_star),
                                                length(theta_star)))

    std_errs <- sqrt(pmax(0, diag(hinv)))
    z_crit <- qnorm(probs)
    intervals <- matrix(NA, nrow = length(theta_star), ncol = length(probs))
    for(i in 1:length(probs)){
      intervals[,i] <- theta_star + z_crit[i] * std_errs
    }

    coef_table <- cbind(Estimate = theta_star, intervals)
    colnames(coef_table) <- c("Estimate", prob_names)
    rownames(coef_table) <- c(colnames(object$data$X),
                              "sigma2", "rho", "lambda")

    loglik_val <- object$estimates$loglik
    logpost_val <- object$estimates$logpost
  }

  res <- list(
    call = object$call,
    coefficients = coef_table,
    loglik = loglik_val,
    logpost = logpost_val,
    model = "IV (SARAR 1,1)",
    method = object$method,
    converged = if(object$method == "normal") object$estimates$converged else TRUE
  )
  class(res) <- "summary.nam"
  return(res)
}
















