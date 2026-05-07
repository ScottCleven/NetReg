#' hane - A function to perform the Homophily Adjusted Network Effects Model
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
#' @param Usample a sample of latent locations to use as priors; an array of
#' size (size, n, D), where size is the number of sample, n is the number of
#' nodes in the network, and D is the number of latent dimension.
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
#' @param gamma_prior The prior for the gamma vector of the model. Requires list
#' or function inputs. Default is a multivariate normal with mean=0 and
#' sigma=2.25^2. If you supply a list, you must have one element named
#' mean' which is the mean of the multivariate normal and one named 'sigma'
#' which is the variance-covariance matrix. If mean is a single integer,
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
#' @param Usample.maxiter a integer value indicating the maximum iteration used
#' in approximating Usample with a matrix normal distribution.
#' @param Usample.eps a numeric value indicating the precision used in
#' approximating Usample with a matrix normal distribution.
#' @param init_vec an initialization vector for \eqn{(\beta, \gamma,}
#' \ifelse{html}{\out{&sigma;<sup>2</sup>}}{\eqn{\sigma^2}}, \eqn{\rho)} in
#' \code{\link{optim}}. If not provided, the initial values are estimated using
#' the two-stage least squares (TSLS) estimator proposed by Kelejian & Prucha
#' (1998).
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

hane <- function(formula, data=NULL, network, Usample, method,
                 beta_prior = list(mean=0, sigma=2.25^2),
                 gamma_prior = list(mean=0, sigma=2.25^2),
                 s2_prior = list(shape=2, scale=2),
                 rho_prior = list(mean=0.36, sd=0.7),
                 rownorm=TRUE, ...){

  # 1. Data Prep
  X <- model.matrix(formula, data=data)
  y <- model.frame(formula, data)[,1]
  n <- NROW(y)
  p <- NCOL(X)
  A <- convert_adjacency(network, "matrix")$A
  if(rownorm) A <- row_norm(A)

  # 2. U-Matrix Approximation (Required for HANE)
  # This uses the Uprior_appox logic from Pham & Sewell
  Uapprox <- Uprior_appox(Usample)
  D <- dim(Usample)[3]

  # 3. Consolidate Data
  hanedat <- list(y=y, X=X, A=A, p=p, n=n, D=D,
                  Lambda=Uapprox$Lambda, Omega=Uapprox$Omega, Psi=Uapprox$Psi,
                  beta_prior=beta_prior, gamma_prior=gamma_prior,
                  s2_prior=s2_prior, rho_prior=rho_prior, ...)

  # 4. Dispatch
  if(tolower(method) == "normal"){
    haneout <- hane_normal(hanedat)
  }else if(tolower(method)=="stan"){
    haneout <- hane_stan(nemdat)
  }else if((tolower(method)=="is") | (tolower(method)=="importance sampling")){
    haneout <- hane_IS(nemdat)
  }

  # 5. Format Output (S3 classes)
  final_out <- list(
    estimates = haneout,
    formula = formula,
    model = "hane",
    method = tolower(method),
    data = list(y=y, X=X, network=A, Usample=Usample),
    priors = list(beta=beta_prior, gamma=gamma_prior, s2=s2_prior, rho=rho_prior),
    call = match.call()
  )
  class(final_out) <- c("nam", "hane")
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

  # The 'M' matrix (I - rho*A)^{-1}
  M <- qr.solve(I - rho*A)

  # Calculate the combined variance structure from the latent space approximation
  gPg <- as.numeric(t(gamma) %*% Psi %*% gamma)

  # Total internal variance before the network transformation
  Sigma_internal <- gPg * Omega + s2 * I

  # The mean structure: M %*% (X %*% beta + Lambda %*% gamma)
  mu_hane <- M %*% (X %*% beta + Lambda %*% gamma)

  # The total covariance: M %*% Sigma_internal %*% M'
  # We use tcrossprod(M %*% Sigma_internal, M) for efficiency
  V_hane <- tcrossprod(M %*% Sigma_internal, M)

  # Use mvtnorm for the likelihood calculation
  mvtnorm::dmvnorm(y, mean = as.vector(mu_hane), sigma = V_hane, log = log)
}




# Posterior of the HANE model.
# Posterior of the HANE model
hane_posterior <- function(theta, A, y, X, Lambda, Omega, Psi,
                           bprior, gprior, s2prior, rhoprior, log=TRUE){

  p <- NCOL(X)
  D <- NCOL(Lambda)

  # Extract parameters for prior functions
  beta <- theta[1:p]
  gamma <- theta[(p+1):(p+D)]
  s2 <- theta[p+D+1]
  rho <- theta[p+D+2]

  # 1. Calculate Log-Likelihood
  # We force log=TRUE here to perform stable addition
  ll <- hane_likelihood(theta, A, y, X, Lambda, Omega, Psi, log=TRUE)

  # 2. Calculate Log-Priors
  # These are the functions created in the engine (normal or IS)
  lp_beta  <- bprior(beta, log=TRUE)
  lp_gamma <- gprior(gamma, log=TRUE)
  lp_s2    <- s2prior(s2, log=TRUE)
  lp_rho   <- rhoprior(rho, log=TRUE)

  # 3. Combine
  out <- ll + lp_beta + lp_gamma + lp_s2 + lp_rho

  # If the user explicitly wants the raw probability (not recommended for HANE)
  if(!log){
    out <- exp(out)
  }

  return(out)
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

# Code for when method="normal"
hane_normal <- function(hanedat){

  # 1. Prior Setup
  if(is.list(hanedat$beta_prior)){
    if(length(hanedat$beta_prior$mean) == 1){
      b_mean <- rep(hanedat$beta_prior$mean, hanedat$p)
    }else{
      b_mean <- hanedat$beta_prior$mean
    }
    b_sigma <- if(length(hanedat$beta_prior$sigma) == 1){
      b_sigma <- diag(hanedat$p) * hanedat$beta_prior$sigma
    }else{
      b_sigma <- hanedat$beta_prior$sigma
    }
    bprior <- function(beta, log=TRUE){
      mvtnorm::dmvnorm(beta, b_mean, b_sigma, log=log)
      }
  }else{
    bprior <- hanedat$beta_prior
  }

  # Gamma Prior (Latent Space Coefficients)
  if(is.list(hanedat$gamma_prior)){
    if(length(hanedat$gamma_prior$mean) == 1){
      g_mean <- rep(hanedat$gamma_prior$mean, hanedat$D)
    }else{
      g_mean <- hanedat$gamma_prior$mean
    }
    if(length(hanedat$gamma_prior$sigma) == 1){
      g_sigma <- diag(hanedat$D) * hanedat$gamma_prior$sigma
    }else{
      g_sigma <- hanedat$gamma_prior$sigma
    }
    gprior <- function(gamma, log=TRUE){
      mvtnorm::dmvnorm(gamma, g_mean, g_sigma, log=log)
      }
  }else{
    gprior <- hanedat$gamma_prior
  }

  # S2 Prior
  if(is.list(hanedat$s2_prior)){
    s2prior <- function(s2, log=TRUE){
      out <- invgamma::dinvgamma(s2, hanedat$s2_prior$shape, hanedat$s2_prior$scale)
      if(log){
        return(log(out))
      }else{
        return(out)
      }
    }
  }else{
    s2prior <- hanedat$s2_prior
  }

  # Rho Prior
  if(is.list(hanedat$rho_prior)){
    rhoprior <- function(rho, log=TRUE){
      truncnorm::dtruncnorm(rho, a=-1, b=1, mean=hanedat$rho_prior$mean,
                            sd=hanedat$rho_prior$sd)
    }
  }else{
    rhoprior <- hanedat$rho_prior
  }

  # 2. Initialization (2SLS to get the optimizer in the ballpark)
  fit2sls <- lnam2sls_effect(hanedat$y, cbind(hanedat$X, hanedat$Lambda),
                             hanedat$A)
  # par structure: [beta, gamma, s2, rho]
  init_vec <- c(fit2sls$coefs[1:(hanedat$p + hanedat$D)], fit2sls$s2,
                fit2sls$coefs[hanedat$p + hanedat$D + 1])

  # 3. Optimization
  # We use L-BFGS-B to enforce the -1 to 1 constraint on rho and >0 on s2
  myopt <- optim(par = init_vec,
                 fn = hane_posterior,
                 method = "L-BFGS-B",
                 lower = c(rep(-Inf, hanedat$p + hanedat$D), 1e-5, -0.999),
                 upper = c(rep(Inf, hanedat$p + hanedat$D + 1), 0.999),
                 control = list(fnscale = -1),
                 hessian = TRUE,
                 y = hanedat$y, X = hanedat$X, A = hanedat$A,
                 Lambda = hanedat$Lambda, Omega = hanedat$Omega, Psi = hanedat$Psi,
                 bprior = bprior, gprior = gprior, s2prior = s2prior, rhoprior = rhoprior)

  # 4. Output formatting
  list(
    beta = myopt$par[1:hanedat$p],
    gamma = myopt$par[(hanedat$p + 1):(hanedat$p + hanedat$D)],
    s2 = myopt$par[hanedat$p + hanedat$D + 1],
    rho = myopt$par[hanedat$p + hanedat$D + 2],
    hessian = myopt$hessian,
    loglik = hane_likelihood(myopt$par, hanedat$A, hanedat$y, hanedat$X,
                             hanedat$Lambda, hanedat$Omega, hanedat$Psi),
    logpost = myopt$value,
    converged = (myopt$convergence == 0)
  )
}

# Code for when method="IS"
hane_IS <- function(hanedat) {
  if(is.null(hanedat$samples)) hanedat$samples <- 10000
  S <- hanedat$samples
  p <- hanedat$p
  D <- hanedat$D

  # 1. Sample from Priors (Proposal Distribution)

  # Beta
  if(is.list(hanedat$beta_prior)) {
    if(length(hanedat$beta_prior$mean) == 1){
      b_mu <- rep(hanedat$beta_prior$mean, p)
    }else{
      b_mu <- hanedat$beta_prior$mean
    }
    if(length(hanedat$beta_prior$sigma) == 1){
      b_sig <- diag(p) * hanedat$beta_prior$sigma
    }else{
      b_sig <- hanedat$beta_prior$sigma
    }
    beta_samples <- mvtnorm::rmvnorm(S, b_mu, b_sig)
    bprior <- function(beta, log=TRUE){
      mvtnorm::dmvnorm(beta, b_mu, b_sig, log=log)
    }
  }else{
    stop("Importance sampling requires list-based priors for HANE.")
  }

  # Gamma (Latent space coefficients)
  if(is.list(hanedat$gamma_prior)){
    if(length(hanedat$gamma_prior$mean) == 1){
      g_mu <- rep(hanedat$gamma_prior$mean, D)
    }else{
      g_mu <- hanedat$gamma_prior$mean
    }
    if(length(hanedat$gamma_prior$sigma) == 1){
      g_sig <- diag(D) * hanedat$gamma_prior$sigma
    }else{
      g_sig <- hanedat$gamma_prior$sigma
    }
    gamma_samples <- mvtnorm::rmvnorm(S, g_mu, g_sig)
    gprior <- function(gamma, log=TRUE){
      mvtnorm::dmvnorm(gamma, g_mu, g_sig, log=log)
    }
  }else{
    stop("Importance sampling requires list-based priors for HANE.")
  }

  # S2 (Variance)
  s2_samples <- invgamma::rinvgamma(S, hanedat$s2_prior$shape,
                                    hanedat$s2_prior$scale)
  s2prior <- function(s2, log=TRUE){
    out <- invgamma::dinvgamma(s2, hanedat$s2_prior$shape,
                               hanedat$s2_prior$scale)
    if(log){
      log(out)
    }else{
      out
    }
  }

  # Rho (Network Correlation)
  rho_samples <- truncnorm::rtruncnorm(S, a = -1, b = 1,
                                       mean = hanedat$rho_prior$mean,
                                       sd = hanedat$rho_prior$sd)
  rhoprior <- function(rho, log=TRUE){
    truncnorm::dtruncnorm(rho, a=-1, b=1, mean=hanedat$rho_prior$mean,
                          sd=hanedat$rho_prior$sd, log=log)
  }

  # 2. Calculate Importance Weights (Log-Likelihood)
  log_weights <- numeric(S)
  for(i in 1:S) {
    theta <- c(beta_samples[i,], gamma_samples[i,], s2_samples[i],
               rho_samples[i])
    log_weights[i] <- hane_likelihood(theta, hanedat$A, hanedat$y, hanedat$X,
                                      hanedat$Lambda, hanedat$Omega,
                                      hanedat$Psi, log = TRUE)
  }

  # 3. Resample using Log-Sum-Exp Trick
  max_log_w <- max(log_weights)
  weights <- exp(log_weights - max_log_w)
  weights <- weights / sum(weights)

  indices <- sample(1:S, size = S, replace = TRUE, prob = weights)

  # Extract resampled parameters
  resampled_beta <- beta_samples[indices, ]
  resampled_gamma <- gamma_samples[indices, ]
  resampled_s2 <- s2_samples[indices]
  resampled_rho <- rho_samples[indices]

  theta_mean <- c(colMeans(resampled_beta), colMeans(resampled_gamma),
                  mean(resampled_s2), mean(resampled_rho))

  # 4. Format Output
  list(
    beta = colMeans(resampled_beta),
    gamma = colMeans(resampled_gamma),
    s2 = mean(resampled_s2),
    rho = mean(resampled_rho),
    samples = list(beta = resampled_beta, gamma = resampled_gamma,
                   s2 = resampled_s2, rho = resampled_rho),
    logpost = hane_posterior(theta_mean, hanedat$A, hanedat$y, hanedat$X,
                             hanedat$Lambda, hanedat$Omega, hanedat$Psi,
                             bprior, gprior, s2prior, rhoprior),
    loglik = hane_likelihood(theta_mean, hanedat$A, hanedat$y, hanedat$X,
                             hanedat$Lambda, hanedat$Omega, hanedat$Psi,
                             log=TRUE),
    log_marginal_lik = max_log_w + log(mean(exp(log_weights - max_log_w))),
    converged = TRUE
  )
}




# Code for when method="stan"
hane_stan <- function(hanedat) {
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

  # Prepare prior vectors/matrices
  if(length(hanedat$beta_prior$mean) == 1){
    b_mu <- rep(hanedat$beta_prior$mean, hanedat$p)
  }else{
    b_mu <- as.vector(hanedat$beta_prior$mean)
  }
  if(length(hanedat$beta_prior$sigma) == 1){
    b_sig <- diag(hanedat$p) * hanedat$beta_prior$sigma
  }else{
    b_sig <- as.matrix(hanedat$beta_prior$sigma)
  }

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

  stan_data <- list(
    N = hanedat$n, K = hanedat$p, D = hanedat$D,
    y = hanedat$y, X = hanedat$X, A = hanedat$A,
    Lambda = hanedat$Lambda, Omega = hanedat$Omega, Psi = hanedat$Psi,
    b_mean = b_mu, b_sigma = b_sig,
    g_mean = g_mu, g_sigma = g_sig,
    s2_shape = hanedat$s2_prior$shape, s2_scale = hanedat$s2_prior$scale,
    rho_mean = hanedat$rho_prior$mean, rho_sd = hanedat$rho_prior$sd
  )

  fit <- rstan::sampling(stanmodels$hane_model, data = stan_data,
                         iter = iters, chains = chains)
  return(fit)
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



#' Summary method for HANE objects
#' @export
summary.hane <- function(object, probs = c(0.025, 0.975), ...){
  p <- ncol(object$data$X)
  D <- ncol(object$data$Usample) # Dimensions of latent space
  prob_names <- paste0(probs * 100, "%")

  if(object$method == "stan"){
    s_mat <- rstan::summary(object$estimates, probs = probs)$summary

    # Identify indices for HANE parameters
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
    loglik_val <- hane_likelihood(theta_star, object$data$network,
                                  object$data$y,
                                  object$data$X, object$estimates$Lambda,
                                  object$estimates$Omega, object$estimates$Psi)

  }else if(object$method %in% c("importance", "is")){
    samples <- object$estimates$samples
    # Combine samples: Beta, Gamma, Sigma2, Rho
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

  }else{
    # Normal Approximation / L-BFGS-B
    theta_star <- c(object$estimates$beta, object$estimates$gamma,
                    object$estimates$s2, object$estimates$rho)

    # Robust Hessian inversion using your established tryCatch pattern
    hinv <- tryCatch(solve(-object$estimates$hessian),
                     error = function(e) matrix(NA, length(theta_star),
                                                length(theta_star)))

    # Ensure no negative variances before sqrt
    std_errs <- sqrt(pmax(0, diag(hinv)))

    z_crit <- qnorm(probs)
    intervals <- matrix(NA, nrow = length(theta_star), ncol = length(probs))
    for(i in 1:length(probs)) {
      intervals[,i] <- theta_star + z_crit[i] * std_errs
    }

    coef_table <- cbind(Estimate = theta_star, intervals)
    colnames(coef_table) <- c("Estimate", prob_names)

    # Labeling
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
    model = "HANE (Homophily-Adjusted)",
    method = object$method,
    converged = if(object$method == "normal") object$estimates$converged else TRUE
  )
  class(res) <- "summary.nam"
  return(res)
}
















