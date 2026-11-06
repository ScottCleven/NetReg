#' nmr - A function to perform Network Moderated Regression
#'
#' @description This function performs a Network Moderated Regression (NMR),
#' which uses community detection to assign nodes to groups and then fits a
#' hierarchical mixed-effects model where the influence of covariates is
#' moderated by community membership.
#'
#' @param formula an object of class "formula" (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible
#' by as.data.frame to a data frame) containing the variables in the model.
#' @param network An adjacency matrix describing your network. Can be of the
#' class: c("matrix", "dgCMatrix", "igraph", "network") or an edgelist. Used for
#' automatic community detection if \code{communities} is NULL.
#' @param family A character string or family object specifying the error distribution
#' and link function. Supported options include \code{"gaussian"}, \code{"binomial"},
#' and \code{"poisson"}. Defaults to \code{"gaussian"}.
#' @param method A character indicating the method to use for estimation.
#' Options for NMR include \code{"stan"}, \code{"normal"}, and \code{"IS"}.
#' @param moderated_formula A half formula line to designate which of the
#' variables in the model to be estimated as moderated e.g. =~ x1 + x3. If
#' NULL (default) will assume all variables in formula input are moderated. It
#' is our suggestion to leave it as NULL and let the Bayes factors inform which
#' variables are being moderated on or not.
#' @param communities An optional vector of community assignments. If NULL
#' (default), communities are detected using \code{igraph::cluster_label_prop}.
#' @param test Logical. If TRUE and \code{method="stan"}, the function calculates
#' Bayes Factors for each moderated variable to test the significance of the
#' community-level moderation using bridge sampling. Defaults to TRUE.
#' @param beta_prior The prior for the fixed effects (beta). Requires list
#' or function inputs. Default is a multivariate normal with mean=0 and sigma=5.
#' If input is a list, must contain elements 'mean' and 'sigma'. If
#' \code{method="stan"}, or \code{method="IS"}, only list input is accepted.
#' @param psi_prior The prior for the log-relative scale of random effects
#' (psi). Requires list or function inputs. Defaults to a normal distribution
#' with mean=0 and sd=1. If input is a list, must contain elements 'mean' and
#' 'sd'. If \code{method="stan"}, or \code{method="IS"}, only list input is
#' accepted.
#' @param re_corr_prior The prior for the off diagonal elements of the
#' correlation matrix of random effects. Requires list or function inputs.
#' Defaults to an LKJ with shape=1 (identity/uniform). If input is a list, must
#' contain element 'shape'. If \code{method="stan"},or \code{method="IS"}, only
#' list input is accepted.
#' @param s2_prior Prior for the variance parameter. Requires list or function
#' inputs. Defaults to a gamma distribution with shape=2 and scale=0.5.
#' If Input is a list, the list must have two elements named 'shape' and 'scale'
#' to represent the shape and scale of the gamma distribution respectively.
#' (rate won't be accepted, it's just 1/scale if you want to use it) Or the user
#' can specify their own prior function with which to input. If
#' \code{method="stan"}, or \code{method="IS"}, only list input is accepted.
#' @param test Logical of whether to compute Bayes factors when method="stan".
#' other methods this will be ignored. TRUE (default) performs the computations.
#' Useful for computation speed.
#' @param ... Additional arguments to be passed to the specific estimation
#' function (e.g., \code{iter}, \code{chains} for Stan, or \code{control}
#' for \code{optim}).
#'
#' @return A list with two classes, \code{"nmr"} and \code{"netreg"}, containing:
#' @return \code{estimates} The raw estimate output (Stan fit, optim list, or IS samples).
#' @return \code{formula} The symbolic formula used.
#' @return \code{model} String identifying the model type ("nmr").
#' @return \code{method} String identifying the estimation method.
#' @return \code{data} List containing the processed Y, X, Z matrices and community vector.
#' @return \code{call} The matched call.
#'
#' @import mvtnorm
#' @import igraph
#' @import lme4
#' @importFrom bridgesampling bridge_sampler bf
#' @importFrom stats model.matrix model.frame as.formula rgamma rnorm
#' @importFrom truncnorm dtruncnorm
#' @importFrom invgamma dinvgamma
#' @importFrom sna lnam
#' @importFrom cmdstanr cmdstan_model write_stan_file
#' @export

nmr <- function(formula, data = NULL, network, family = "gaussian",
                method = "stan", moderated_formula = NULL, communities = NULL,
                beta_prior,
                psi_prior,
                re_corr_prior = list(shape=1),
                s2_prior = list(shape=0.0001, scale=0.0001),
                test=TRUE, compile = TRUE, ...){

  # Extract all individual variables needed across both formulas
  main_vars <- all.vars(formula)
  if(!is.null(moderated_formula)){
    mod_vars <- all.vars(moderated_formula)
  }else{
    mod_vars <- c()
  }
  vars_in_formula <- unique(c(main_vars, mod_vars))

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

  # Resolve family object
  if(is.character(family)){
    family_obj <- get(family, mode = "function", envir = parent.frame())()
  }else if (is.function(family)){
    family_obj <- family()
  }else{
    family_obj <- family
  }

  X <- model.matrix(formula, data = data)
  mf_Y <- model.frame(formula, data = data)[, 1]

  if(family_obj$family == "binomial"){
    if(is.matrix(mf_Y) && ncol(mf_Y) == 2){
      Y <- mf_Y
    }else{
      Y <- cbind(as.vector(mf_Y), 1 - as.vector(mf_Y))
      colnames(Y) <- c("successes", "failures")
    }
  }else{
    Y <- as.vector(mf_Y)
  }

  # Community logic
  if(is.null(communities)){
    A_norm <- row_norm(row_norm(network), reverse = TRUE)
    g <- convert_network(A_norm, out="igraph")$A
    communities <- igraph::cluster_label_prop(g)$membership
  }

  # Moderator logic
  if(is.null(moderated_formula)){
    X2 <- X
  }else{
    X2 <- model.matrix(moderated_formula, data = data)
  }

  # Build Adaptive Scale-Invariant Prior for Psi
  if(missing(psi_prior)){
    message("\nThe psi_prior parameter is not specified. It will be set
            automatically based on log-data scales.")
    has_intercept <- all(X2[, 1] == 1)
    X2_slopes <- if(has_intercept) X2[, -1, drop = FALSE] else X2
    s_X2_slopes <- apply(X2_slopes, 2, sd)
    s_X2_slopes[s_X2_slopes == 0] <- 1

    slope_sds_psi <- log(1 + s_X2_slopes)

    if(has_intercept){
      intercept_sd_psi <- 1.5
      psi_sd_vec <- c(intercept_sd_psi, slope_sds_psi)
    } else {
      psi_sd_vec <- slope_sds_psi
    }
    psi_prior <- list(
      mean = rep(0, ncol(X2)),
      sd   = psi_sd_vec
    )
  }


  # Default Prior for Beta
  if(missing(beta_prior)){
    message("\nThe beta_prior parameter is not specified. It will be set
            automatically based on link-data scales.")
    s_j <- apply(X[, -1, drop = FALSE], 2, sd)
    s_j[s_j == 0] <- 1
    if(ncol(X) > 1){
      if (family_obj$family == "poisson"){
        s_y_link <- sd(log(Y + 0.5))
        intercept_var <- max(2.5^2, (2.5 * s_y_link)^2)
      }else if(family_obj$family == "binomial"){
        emp_log_odds <- log((Y[, 1] + 0.5) / (Y[, 2] + 0.5))
        s_y_link     <- sd(emp_log_odds)
        intercept_var <- max(2.5^2, (2.5 * s_y_link)^2)
      }else{
        s_y <- sd(Y)
        intercept_var <- max(2.5^2, (2.5 * s_y)^2)
      }
      slope_vars    <- 2.5^2 * s_j^2
      Sigma_beta_mat <- diag(c(intercept_var, slope_vars), nrow = ncol(X))
    }else{
      intercept_var <- 2.5^2
      Sigma_beta_mat <- matrix(intercept_var, nrow = 1, ncol = 1)
    }
    beta_prior <- list(
      mean  = rep(0, ncol(X)),
      sigma = Sigma_beta_mat
    )
  }


  # list limits for importance sampling
  is_importance <- tolower(method) %in% c("is", "importance sampling")
  if (is_importance && (typeof(beta_prior) != "list" ||
                        typeof(psi_prior) != "list" ||
                        typeof(re_corr_prior) != "list" ||
                        typeof(s2_prior) != "list")) {
    stop("importance sampling methods require list input of prior parameters,
         prior functions cannot be changed.")
  }

  nmrdat <- list(Y=Y, X=X, X2=X2, communities=communities,
                 family=family_obj, moderated_formula=moderated_formula,
                 beta_prior=beta_prior, psi_prior=psi_prior, s2_prior=s2_prior,
                 re_corr_prior=re_corr_prior, p=ncol(X), m=ncol(X2),
                 test=test, compile = compile, ...)

  if(tolower(method) == "stan"){
    nmrout <- nmr_stan(nmrdat, ...)
  }else if(tolower(method) == "normal"){
    nmrout <- nmr_normal(nmrdat, ...)
  }else if(tolower(method) %in% c("importance", "is", "importance sampling")){
    if(family_obj$family != "gaussian"){
      stop("Importance sampling is only supported for family = 'gaussian'")
    }else{
      nmrout <- nmr_IS(nmrdat, ...)
    }
  }else{
    stop("Unknown estimation method.")
  }

  # Return string immediately if compile is FALSE
  if (tolower(method) == "stan" && !compile) {
    return(nmrout)
  }

  final_out <- list(
    estimates = nmrout,
    formula = formula,
    moderated_formula = moderated_formula,
    family = family_obj,
    model = "nmr",
    method = tolower(method),
    data = list(Y = Y, X = X, X2 = X2, communities = communities,
                network = convert_network(network, "Matrix")$A),
    call = match.call()
  )
  class(final_out) <- c("nmr", "netreg")
  return(final_out)
}




# Log-Likelihood for Network Moderated Regression
nmr_likelihood <- function(theta, Y, X, X2, communities,
                           family_obj, log = TRUE){
  # Parameters
  p <- ncol(X)
  m <- ncol(X2)
  J <- as.integer(as.factor(communities))
  K_comm <- length(unique(J))

  beta <- theta[1:p]
  psi  <- theta[(p + 1):(p + m)]

  has_sigma <- (family_obj$family %in% c("gaussian", "Gamma",
                                         "inverse.gaussian"))

  if(has_sigma){
    sigma <- theta[length(theta)]
    if(sigma <= 0) return(if(log) -Inf else 0)
    b_vec <- theta[(p + m + 1):(length(theta) - 1)]
  }else{
    b_vec <- theta[(p + m + 1):length(theta)]
  }

  # Reshape unscaled random variations into community matrix lines
  b_matrix <- matrix(b_vec, nrow = K_comm, ncol = m)

  # Compute Linear Predictor (eta) vectorization
  eta <- as.vector(X %*% beta) + rowSums(b_matrix[J, , drop = FALSE] * X2)

  # Map to conditional mean via inverse link function
  mu <- family_obj$linkinv(eta)

  # Calculate Density based on GLM Family
  if(family_obj$family == "gaussian"){
    out <- dnorm(Y, mean = mu, sd = sigma, log = TRUE)
  }else if(family_obj$family == "binomial"){
    trials <- Y[, 1] + Y[, 2]
    out <- dbinom(Y[, 1], size = trials, prob = mu, log = TRUE)
  }else if(family_obj$family == "poisson"){
    out <- dpois(Y, lambda = mu, log = TRUE)
  }else{
    stop(paste("Family", family_obj$family, "is not supported."))
  }

  if(log){
    return(sum(out))
  }else{
    return(exp(sum(out)))
  }
}

# Log-Posterior for Network Moderated Regression (NMR)
nmr_posterior <- function(theta, Y, X, X2, communities, family_obj, log = TRUE,
                          bprior, psiprior, reprior, sigmaprior){
  p <- ncol(X)
  m <- ncol(X2)
  J <- as.integer(as.factor(communities))
  K_comm <- length(unique(J))

  has_sigma <- (family_obj$family %in% c("gaussian", "Gamma",
                                         "inverse.gaussian"))

  beta <- theta[1:p]
  psi  <- theta[(p + 1):(p + m)]

  if(has_sigma){
    sigma <- theta[length(theta)]
    if(sigma <= 0) return(if(log) -Inf else 0)
    b_vec <- theta[(p + m + 1):(length(theta) - 1)]
  }else{
    b_vec <- theta[(p + m + 1):length(theta)]
  }
  b_matrix <- matrix(b_vec, nrow = K_comm, ncol = m)

  # Log-Likelihood
  log_lik <- nmr_likelihood(theta, Y, X, X2, communities, family_obj,
                            log = TRUE)
  if(!is.finite(log_lik)){
    if(log){
      return(-Inf)
    }else{
      return(0)
    }
  }

  # Penalty (Independent Random Effects)
  log_prior_re_hier <- 0
  for(j in 1:m){
    col_sd <- if (has_sigma) sigma * exp(psi[j]) else exp(psi[j])
    log_prior_re_hier <- log_prior_re_hier + sum(dnorm(b_matrix[, j], 0,
                                                       col_sd, log = TRUE))
  }

  # Fixed Parameter Log-Priors
  log_prior_beta  <- bprior(beta)
  log_prior_psi   <- psiprior(psi)
  log_prior_re    <- reprior(b_matrix)

  log_prior_sigma <- 0
  if(has_sigma){
    log_prior_sigma <- sigmaprior(sigma)
  }

  total <- log_lik + log_prior_re_hier + log_prior_beta + log_prior_psi +
    log_prior_re + log_prior_sigma

  if(log){
    return(total)
  }else{
    return(exp(total))
  }
}


# Code for when method="normal"
nmr_normal <- function(nmrdat, ...){
  # Parameter Layout
  p <- nmrdat$p
  m <- nmrdat$m
  J <- as.integer(as.factor(nmrdat$communities))
  K_comm <- length(unique(J))

  has_sigma <- (nmrdat$family$family %in% c("gaussian", "Gamma",
                                            "inverse.gaussian"))

  # Priors
  if(typeof(nmrdat$beta_prior) == "list"){
    if(length(nmrdat$beta_prior$mean) == 1){
      bmean <- rep(nmrdat$beta_prior$mean, p)
    }else{
      bmean <- nmrdat$beta_prior$mean
    }
    if(length(nmrdat$beta_prior$sigma) == 1){
      bsigma <- diag(p) * nmrdat$beta_prior$sigma
    }else{
      bsigma <- nmrdat$beta_prior$sigma
    }
    bprior <- function(beta){
      mvtnorm::dmvnorm(beta, bmean, bsigma, log = TRUE)
    }
  }else{
    bprior <- nmrdat$beta_prior
  }

  if(typeof(nmrdat$psi_prior) == "list"){
    psiprior <- function(psi){
      sum(dnorm(psi, nmrdat$psi_prior$mean, nmrdat$psi_prior$sd, log = TRUE))
    }
  }else{
    psiprior <- nmrdat$psi_prior
  }

  # Regularizing prior for the unscaled random variations matrix b_matrix
  if(typeof(nmrdat$re_corr_prior) == "list"){
    reprior <- function(b_matrix){
      sum(dnorm(b_matrix, 0, 1, log = TRUE))
    }
  }else{
    reprior <- nmrdat$re_corr_prior
  }

  if(typeof(nmrdat$s2_prior) == "list"){
    sigmaprior <- function(sigma){
      dgamma(sigma, shape = nmrdat$s2_prior$shape,
             scale = nmrdat$s2_prior$scale, log = TRUE)
    }
  }else{
    sigmaprior <- nmrdat$s2_prior
  }

  # Initialization via lme4
  lmer_df <- data.frame(G = as.factor(nmrdat$communities))
  lmer_df$Y <- nmrdat$Y
  for(i in 1:p){
    lmer_df[[paste0("X", i)]] <- nmrdat$X[, i]
  }
  for(i in 1:m){
    lmer_df[[paste0("X2_", i)]] <- nmrdat$X2[, i]
  }

  x_vars <- paste0("X", 1:p, collapse = " + ")
  z_vars <- paste0("(0 + X2_", 1:m, " | G)", collapse = " + ")
  lmer_form <- as.formula(paste0("Y ~ 0 + ", x_vars, " + ", z_vars))

  if(has_sigma){
    lmer_fit <- lme4::lmer(lmer_form, data = lmer_df)
    sig_init <- sigma(lmer_fit)
    psi_init <- log(as.data.frame(lme4::VarCorr(lmer_fit))$sdcor[1:m] /
                      sig_init + 1e-15)
    init_theta <- c(
      lme4::fixef(lmer_fit),
      psi_init,
      as.vector(as.matrix(lme4::ranef(lmer_fit)$G)),
      sig_init
    )
    if(!is.null(nmrdat$lower)){
      lower <- nmrdat$lower
    }else{
      lower <- c(rep(-Inf, p), rep(-Inf, m), rep(-Inf, K_comm * m), 1e-6)
      # Goes c(beta, psi, b, sigma)
    }

    if(!is.null(nmrdat$upper)){
      upper <- nmrdat$upper
    }else{
      upper <- c(rep(Inf, p), rep(Inf, m), rep(Inf, K_comm * m), Inf)
    }
  }else{
    lmer_fit <- lme4::glmer(lmer_form, data = lmer_df, family = nmrdat$family)
    psi_init <- log(as.data.frame(lme4::VarCorr(lmer_fit))$sdcor[1:m] + 1e-15)
    init_theta <- c(
      lme4::fixef(lmer_fit),
      psi_init,
      as.vector(as.matrix(lme4::ranef(lmer_fit)$G))
    )
    if(!is.null(nmrdat$lower)){
      lower <- nmrdat$lower
    }else{
      lower <- c(rep(-Inf, p), rep(-Inf, m), rep(-Inf, K_comm * m))
      # Goes c(beta, psi, b)
    }

    if(!is.null(nmrdat$upper)){
      upper <- nmrdat$upper
    }else{
      upper <- c(rep(Inf, p), rep(Inf, m), rep(Inf, K_comm * m))
    }
  }

  # Checking if additional Arguments were input
  if(!is.null(nmrdat$gr)){
    gr=nmrdat$gr
  }else{
    gr=NULL
  }
  if(!is.null(nmrdat$control)){
    control=nmrdat$control
  }else{
    control=list(fnscale = -1)
  }

  # Optimization using L-BFGS-B
  suppressWarnings(
    myopt <- optim(par = init_theta,
                   fn = nmr_posterior,
                   gr = gr,
                   method = "L-BFGS-B",
                   lower = lower,
                   upper = upper,
                   control = control,
                   hessian = TRUE,
                   Y = nmrdat$Y,
                   X = nmrdat$X,
                   X2 = nmrdat$X2,
                   communities = nmrdat$communities,
                   family_obj = nmrdat$family,
                   bprior = bprior,
                   psiprior = psiprior,
                   reprior = reprior,
                   sigmaprior = sigmaprior)
  )

  # Format Output
  beta_ests <- myopt$par[1:p]
  names(beta_ests) <- colnames(nmrdat$X)

  psi_ests <- myopt$par[(p + 1):(p + m)]
  names(psi_ests) <- colnames(nmrdat$X2)

  sigma_val <- if (has_sigma) myopt$par[length(myopt$par)] else NULL

  list(
    beta      = beta_ests,
    psi       = psi_ests,
    b_matrix  = matrix(myopt$par[(p + m + 1):(p + m + (K_comm * m))],
                       K_comm, m),
    sigma     = sigma_val,
    loglik    = nmr_likelihood(myopt$par, nmrdat$Y, nmrdat$X, nmrdat$X2,
                               nmrdat$communities, family_obj = nmrdat$family,
                               log = TRUE),
    logpost   = myopt$value,
    converged = myopt$convergence == 0,
    raw_opt   = myopt
  )
}


# Code for when method="IS"
nmr_IS <- function(nmrdat, ...){
  if(is.null(nmrdat$samples)){
    nmrdat$samples <- 5000
  }
  S <- nmrdat$samples
  p <- nmrdat$p
  m <- nmrdat$m
  J <- as.integer(as.factor(nmrdat$communities))
  K_comm <- length(unique(J))

  # Run fast Normal optimizer first to find posterior peak
  fit_norm <- nmr_normal(nmrdat)
  mu <- fit_norm$raw_opt$par

  # Extract the covariance matrix from the inverted negative Hessian
  # Add a ridge inflation (* 2) to ensure the proposal is wider than the target
  Sigma <- solve(-fit_norm$raw_opt$hessian) * 2

  # Draw proposals from Multivariate Normal distribution
  cat("Generating", S, "samples from the proposal distribution...\n")
  theta_samples <- mvtnorm::rmvnorm(S, mean = mu, sigma = Sigma)

  # Setup standard prior closures for weighting calculations
  if(typeof(nmrdat$beta_prior) == "list"){
    if(length(nmrdat$beta_prior$mean) == 1){
      bmean <- rep(nmrdat$beta_prior$mean, p)
    }else{
      bmean <- nmrdat$beta_prior$mean
    }
    if(length(nmrdat$beta_prior$sigma) == 1){
      bsigma <- diag(p) * nmrdat$beta_prior$sigma
    }else{
      bsigma <- nmrdat$beta_prior$sigma
    }
    bprior <- function(beta){
      mvtnorm::dmvnorm(beta, bmean, bsigma, log = TRUE)
    }
  }else{
    bprior <- nmrdat$beta_prior
  }

  if(typeof(nmrdat$psi_prior) == "list"){
    psiprior <- function(psi){
      sum(dnorm(psi, nmrdat$psi_prior$mean, nmrdat$psi_prior$sd, log = TRUE))
    }
  }else{
    psiprior <- nmrdat$psi_prior
  }

  if(typeof(nmrdat$re_corr_prior) == "list"){
    reprior <- function(b_matrix){
      sum(dnorm(b_matrix, 0, 1, log = TRUE))
    }
  }else{
    reprior <- nmrdat$re_corr_prior
  }

  if(typeof(nmrdat$s2_prior) == "list"){
    sigmaprior <- function(sigma) dgamma(sigma, shape = nmrdat$s2_prior$shape,
                                         scale = nmrdat$s2_prior$scale,
                                         log = TRUE)
  }else{
    sigmaprior <- nmrdat$s2_prior
  }

  # Calculate Normalized Importance Weights
  # Weight = [ Likelihood * Priors ] / [ Proposal_Density ]
  log_weights <- numeric(S)

  cat("Calculating weights...\n")
  for (i in 1:S){
    sigma_curr <- theta_samples[i, ncol(theta_samples)]
    if (sigma_curr <= 0) {
      log_weights[i] <- -Inf
      next
    }

    # Target numerator log-density
    log_post <- nmr_posterior(theta_samples[i, ], nmrdat$Y, nmrdat$X, nmrdat$X2,
                              nmrdat$communities, family_obj = nmrdat$family,
                              log = TRUE, bprior = bprior, psiprior = psiprior,
                              reprior = reprior, sigmaprior = sigmaprior)

    # Proposal denominator log-density
    log_prop <- mvtnorm::dmvnorm(theta_samples[i, ], mean = mu, sigma = Sigma,
                                 log = TRUE)

    log_weights[i] <- log_post - log_prop
  }

  # Normalize using log-sum-exp
  max_log_w <- max(log_weights[is.finite(log_weights)])
  weights <- exp(log_weights - max_log_w)
  weights[is.na(weights)] <- 0
  weights <- weights / sum(weights)

  # Resample to form the empirical posterior distribution
  indices <- sample(1:S, size = S, replace = TRUE, prob = weights)
  resampled_theta <- theta_samples[indices, , drop = FALSE]

  # Calculate the effective sample size
  ess <- 1 / sum(weights^2)
  ess_threshold <- 0.05 * S
  is_converged <- ess >= ess_threshold

  if (!is_converged){
    warning(paste0("Importance sampler exhibits high weight degeneracy! ",
                   "Effective Sample Size (ESS = ", round(ess, 1), ") is below
                   5% of nominal samples."))
  }

  # Extract means from resampled distribution
  beta_means <- colMeans(resampled_theta[, 1:p, drop = FALSE])
  names(beta_means) <- colnames(nmrdat$X)

  psi_means <- colMeans(resampled_theta[, (p + 1):(p + m), drop = FALSE])
  names(psi_means) <- colnames(nmrdat$X2)

  list(
    beta      = beta_means,
    psi       = psi_means,
    sigma     = mean(resampled_theta[, ncol(resampled_theta)]),
    samples   = list(
      beta    = resampled_theta[, 1:p, drop = FALSE],
      psi     = resampled_theta[, (p + 1):(p + m), drop = FALSE],
      b_matrix = resampled_theta[, (p + m + 1):(p + m + (K_comm * m)),
                                 drop = FALSE],
      sigma   = resampled_theta[, ncol(resampled_theta)]
    ),
    log_marginal_lik = max_log_w + log(mean(exp(log_weights - max_log_w))),
    ess       = ess,
    converged = is_converged
  )
}



# Code for when method="stan"
nmr_stan <- function(nmrdat, ...){
  test <- nmrdat$test
  if(!requireNamespace("cmdstanr", quietly = TRUE)){
    stop("Package 'cmdstanr' is needed for this method. Please install it.",
         call. = FALSE)
  }
  if(!requireNamespace("bridgesampling", quietly = TRUE)){
    stop("Package 'bridgesampling' is needed. Please install it.",
         call. = FALSE)
  }

  # Checking for additional inputs
  if(is.null(nmrdat$iterations)){
    iterations=2000
  }else{
    iterations=nmrdat$iterations
  }
  if(is.null(nmrdat$chains)){
    chains=4
  }else{
    chains=nmrdat$chains
  }

  p_dim <- nmrdat$p
  family_name <- nmrdat$family$family

  # Format beta prior
  if(is.list(nmrdat$beta_prior)){
    if(length(nmrdat$beta_prior$mean) == 1){
      b_mean <- rep(nmrdat$beta_prior$mean, p_dim)
    }else{
      b_mean <- nmrdat$beta_prior$mean
    }
    if(length(nmrdat$beta_prior$sigma) == 1){
      b_sigma <- diag(p_dim) * nmrdat$beta_prior$sigma
    }else{
      b_sigma <- as.matrix(nmrdat$beta_prior$sigma)
    }

    b_mean_str <- paste0("[", paste(b_mean, collapse = ", "), "]'")
    b_sigma_rows <- apply(b_sigma, 1,
                          function(row) paste0("[",
                                               paste(row, collapse = ", "),
                                               "]"))
    b_sigma_str <- paste0("[", paste(b_sigma_rows, collapse = ", "), "]")

    beta_prior_line <- paste0("target += multi_normal_lpdf(beta | ",
                              b_mean_str, ", ", b_sigma_str, ");")
  }else if(is.character(nmrdat$beta_prior)){
    beta_prior_line <- nmrdat$beta_prior
  }else{
    stop("beta_prior must be a list or a character string.")
  }

  # Format psi prior cleanly as independent univariate components passed via data block
  if(is.list(nmrdat$psi_prior)){
    psi_prior_line <- "target += normal_lpdf(psi | psi_mean, psi_sd);"
  }else if(is.character(nmrdat$psi_prior)){
    psi_prior_line <- nmrdat$psi_prior
  }else{
    stop("psi_prior must be a list or a character string.")
  }

  # Format re_corr prior
  if(is.list(nmrdat$re_corr_prior)){
    re_corr_prior_line <- paste0("
      if(m == 1){
        target += 0.0;
      }else{
        target += lkj_corr_cholesky_lpdf(L_b | ", nmrdat$re_corr_prior$shape, ");
      }")
  }else if(is.character(nmrdat$re_corr_prior)){
    re_corr_prior_line <- nmrdat$re_corr_prior
  }else{
    stop("re_corr_prior must be a list or a character string.")
  }

  # Format s2 prior
  if(is.list(nmrdat$s2_prior)){
    s2_prior_line <- "target += student_t_lpdf(sigma | 3, 0, 2.5) -
    student_t_lccdf(0 | 3, 0, 2.5);"
  }else if(is.character(nmrdat$s2_prior)){
    s2_prior_line <- nmrdat$s2_prior
  }else{
    stop("s2_prior must be a list or a character string.")
  }

  # Select family-specific structural template
  if(family_name == "gaussian"){
    stan_template <- paste0("
    data {
      int<lower=1> N;
      vector[N] Y;
      int<lower=1> K;
      matrix[N, K] X;
      int<lower=1> K_comm;
      int<lower=1> m;
      array[N] int<lower=1> J;
      matrix[N, m] X2;
      vector[m] psi_mean;
      vector[m] psi_sd;
    }

    parameters {
      vector[K] beta;                  // Fixed effects
      real<lower=0> sigma;             // Residual SD (global scale)
      vector[m] psi;                   // Log-transformed relative scale
      matrix[m, K_comm] z_b;           // Standardized random variation variations
      cholesky_factor_corr[m] L_b;     // Cholesky factor of correlations
    }

    transformed parameters {
      vector[m] lambda = exp(psi);
      matrix[K_comm, m] b_matrix;

      b_matrix = transpose(diag_pre_multiply(lambda * sigma, L_b) * z_b);
    }

    model {
      vector[N] mu = X * beta;

      for (n in 1:N) {
        mu[n] += dot_product(b_matrix[J[n]], X2[n]);
      }

      // Priors dynamically inserted via R
      ", beta_prior_line, "
      ", s2_prior_line, "
      ", psi_prior_line, "
      ", re_corr_prior_line, "
      target += std_normal_lpdf(to_vector(z_b));

      // Likelihood
      target += normal_lpdf(Y | mu, sigma);
    }
    ")
  }else if(family_name %in% c("binomial", "poisson")){
    if(family_name == "binomial"){
      y_decl_line <- "array[N, 2] int<lower=0> Y_int;"
      likelihood_line <- "
        for (n in 1:N){
          target += binomial_logit_lpmf(
          Y_int[n, 1] | Y_int[n, 1] + Y_int[n, 2], eta[n]);
        }"
    }else if(family_name == "poisson"){
      y_decl_line <- "array[N] int<lower=0> Y_int;"
      likelihood_line <- "target += poisson_log_lpmf(Y_int | eta);"
    }
    stan_template <- paste0("
      data {
        int<lower=1> N;
        ", y_decl_line, "
        int<lower=1> K;
        matrix[N, K] X;
        int<lower=1> K_comm;
        int<lower=1> m;
        array[N] int<lower=1> J;
        matrix[N, m] X2;
        vector[m] psi_mean;
        vector[m] psi_sd;
      }

      parameters {
        vector[K] beta;                  // Fixed effects
        vector[m] psi;                   // Log-transformed relative scale
        matrix[m, K_comm] z_b;           // Standardized random variations
        cholesky_factor_corr[m] L_b;     // Cholesky factor of correlations
      }

      transformed parameters {
        vector[m] lambda = exp(psi);
        matrix[K_comm, m] b_matrix;

        b_matrix = transpose(diag_pre_multiply(lambda, L_b) * z_b);
      }

      model {
        vector[N] eta = X * beta;
        for (n in 1:N) {
          eta[n] += dot_product(b_matrix[J[n]], X2[n]);
        }

        // Priors dynamically inserted via R
        ", beta_prior_line, "
        ", psi_prior_line, "
        ", re_corr_prior_line, "
        target += std_normal_lpdf(to_vector(z_b));

        // Discrete Likelihood
        ", likelihood_line, "
      }
      ")
  }else{
    stop(paste("Family", family_name, "is not supported by the Stan engine."))
  }

  if(!is.null(nmrdat$compile) && !nmrdat$compile){
    return(stan_template)
  }

  # Compile the template file exactly ONCE at runtime launch
  stan_file_tmp <- cmdstanr::write_stan_file(stan_template,
                                             force_overwrite = TRUE)
  compiled_model <- cmdstanr::cmdstan_model(stan_file_tmp, force=TRUE)

  fit_and_bridge <- function(X2_matrix){
    current_m <- ncol(X2_matrix)

    if(current_m == nmrdat$m){
      p_mu  <- as.vector(nmrdat$psi_prior$mean)
      if(!is.null(nmrdat$psi_prior$sd)){
        p_sd  <- as.vector(nmrdat$psi_prior$sd)
      }else{
        p_sd  <- sqrt(diag(as.matrix(nmrdat$psi_prior$sigma)))
      }
    }else{
      p_mu  <- rep(0, current_m)
      if(!is.null(nmrdat$psi_prior$sd)){
        all_sds <- as.vector(nmrdat$psi_prior$sd)
      }else{
        all_sds <- sqrt(diag(as.matrix(nmrdat$psi_prior$sigma)))
      }
      p_sd  <- all_sds[1:current_m]
    }

    stan_data <- list(
      N        = NROW(nmrdat$Y),
      K        = p_dim,
      X        = nmrdat$X,
      K_comm   = length(unique(nmrdat$communities)),
      m        = current_m,
      J        = as.integer(as.factor(nmrdat$communities)),
      X2       = X2_matrix,
      psi_mean = p_mu,
      psi_sd   = p_sd
    )

    if(family_name == "gaussian"){
      stan_data$Y <- as.numeric(nmrdat$Y)
    }else{
      if(family_name == "binomial"){
        storage.mode(nmrdat$Y) <- "integer"
        stan_data$Y_int <- nmrdat$Y
      }else{
        stan_data$Y_int <- as.integer(nmrdat$Y)
      }
    }

    if(test){
      if(is.null(nmrdat$iter_sampling)){
        iters <- 10000
      }else{
        iters <- nmrdat$iter_sampling
      }
      if(is.null(nmrdat$iter_warmup)){
        warmup_iters <- 2000
      }else{
        warmup_iters <- nmrdat$iter_warmup
      }
    }else{
      if(is.null(nmrdat$iter_sampling)){
        iters <- 2000
      }else{
        iters <- nmrdat$iter_sampling
      }
      if(is.null(nmrdat$iter_warmup)){
        warmup_iters <- 2000
      }else{
        warmup_iters <- nmrdat$iter_warmup
      }
    }

    init_fun <- function(){
      init_list <- list(
        beta = rep(0, p_dim),
        psi  = rep(0, current_m),
        z_b  = matrix(0, nrow = current_m, ncol = stan_data$K_comm)
      )
      if (current_m >= 2) {
        init_list$L_b <- diag(current_m)
      }
      return(init_list)
    }

    # Extract stan arguments
    all_args <- names(formals(compiled_model$sample))
    valid_sample_args <- setdiff(all_args, "data")
    user_sampling_params <- nmrdat[names(nmrdat) %in% valid_sample_args]

    # Sample
    final_args <- c(
      list(
        data = stan_data,
        iter_sampling = iters,
        iter_warmup = warmup_iters,
        chains = chains,
        refresh = 0,
        init = init_fun,
        sig_figs = 14
      ),
      user_sampling_params            # Overrides or adds choices from user
    )
    final_args <- final_args[!duplicated(names(final_args), fromLast = TRUE)]
    cmdstan_fit <- do.call(compiled_model$sample, final_args)

    draws_matrix <- cmdstan_fit$draws(format = "draws_matrix")
    attr(draws_matrix, "stan_code") <- stan_template

    if(test){
      lml <- bridgesampling::bridge_sampler(cmdstan_fit, silent = TRUE,
                                            maxiter = 10000)
    }else{
      lml <- NULL
    }
    return(list(draws = draws_matrix, lml = lml))
  }

  full_res <- fit_and_bridge(nmrdat$X2)

  bfs <- list()
  if(test){
    X2_full <- nmrdat$X2
    col_names <- colnames(X2_full)

    term_assign <- attr(X2_full, "assign")
    term_labels <- attr(terms(nmrdat$moderated_formula), "term.labels")

    has_intercept <- "(Intercept)" %in% col_names
    unique_terms <- unique(term_assign)
    if(has_intercept){
      unique_terms <- unique_terms[unique_terms != 0]
    }

    if(length(unique_terms) > 0){
      for(t_idx in unique_terms){
        if(t_idx <= length(term_labels)){
          term_name <- term_labels[t_idx]
        }else{
          term_name <- col_names[term_assign == t_idx]
        }

        cat(sprintf("\nTesting moderation of parameter group: %s...\n",
                    term_name))
        X2_reduced <- X2_full[, term_assign != t_idx, drop = FALSE]

        if(ncol(X2_reduced) == 0){
          cat("Cannot compute Bayes Factor:
              reduced model contains 0 random effects.\n")
          next
        }

        red_res <- fit_and_bridge(X2_reduced)
        bfs[[term_name]] <- honest_val <-
          as.numeric(bridgesampling::bf(full_res$lml, red_res$lml)$bf)
      }
    }
  }

  return(list(fit = full_res$draws, bfs = bfs))
}








#' Summary method for Network Moderated Regression models
#' @method summary nmr
#' @export
summary.nmr <- function(object, probs = c(0.025, 0.975), ...){
  p_dim <- ncol(object$data$X)
  m_dim <- ncol(object$data$X2)
  prob_names <- paste0(probs * 100, "%")

  if(object$method == "stan"){
    draws <- object$estimates$fit

    beta_cols   <- grep("^beta(\\[|$)", colnames(draws))
    lambda_cols <- grep("^lambda(\\[|$)", colnames(draws))

    # Compute descriptive columns directly from your matrix draws rows
    b_draws <- draws[, beta_cols, drop = FALSE]
    b_means <- colMeans(b_draws)
    b_quants <- apply(b_draws, 2, quantile, probs = probs, na.rm = TRUE)

    if(length(beta_cols) == 1){
      coef_table <- matrix(c(b_means, b_quants), nrow = 1)
    }else{
      coef_table <- cbind(b_means, t(b_quants))
    }
    colnames(coef_table) <- c("Estimate", prob_names)
    rownames(coef_table) <- colnames(object$data$X)

    # Process and transform random parameters matching your simulation scales
    l_draws <- draws[, lambda_cols, drop = FALSE]
    re_means <- colMeans(l_draws)
    re_conf  <- apply(l_draws, 2, quantile, probs = probs, na.rm = TRUE)

    if(object$family$family == "gaussian"){
      sigma_col <- which(colnames(draws) == "sigma")
      sig_draws <- draws[, sigma_col]
      sig_est   <- mean(sig_draws)
      sig_conf  <- quantile(sig_draws, probs = probs, na.rm = TRUE)

      # True scale standard deviations matrix element transform: lambda * sigma
      rand_means <- re_means * sig_est
      rand_intervals <- t(re_conf) * rep(sig_conf, each = m_dim)

      rand_table <- cbind(mean = rand_means, rand_intervals)
      rownames(rand_table) <- paste0("SD (Community): ",
                                     colnames(object$data$X2))

      sig_row <- matrix(c(sig_est, sig_conf), nrow = 1)
      colnames(sig_row) <- c("mean", prob_names)
      rownames(sig_row) <- "sigma (residual)"

      coef_table <- rbind(coef_table, rand_table, sig_row)
    }else{
      # Discrete models are scaled directly by lambda matrix vectors
      if(length(lambda_cols) == 1){
        rand_table <- matrix(c(re_means, re_conf), nrow = 1)
      }else{
        rand_table <- cbind(re_means, t(re_conf))
      }
      colnames(rand_table) <- c("mean", prob_names)
      rownames(rand_table) <- paste0("SD (Community): ",
                                     colnames(object$data$X2))
      coef_table <- rbind(coef_table, rand_table)
    }

    lp_col <- which(colnames(draws) == "lp__")
    logpost_val <- mean(draws[, lp_col])
    loglik_val  <- NA
    if(length(object$estimates$bfs) > 0){
      bayes_factors <- unlist(object$estimates$bfs)
    }else{
      bayes_factors <- NULL
    }

  }else{
    # Normal Approximation
    theta_star <- c(object$estimates$beta, object$estimates$psi)
    if(!is.null(object$estimates$sigma)) theta_star <- c(theta_star,
                                                         object$estimates$sigma)

    hinv <- tryCatch(solve(-object$estimates$raw_opt$hessian),
                     error = function(e) matrix(NA, length(theta_star),
                                                length(theta_star)))
    std_errs <- sqrt(diag(hinv))[1:length(theta_star)]

    z_crit <- qnorm(probs)
    intervals <- matrix(NA, nrow = length(theta_star), ncol = length(probs))
    for(i in 1:length(probs)) intervals[,i] <- theta_star + z_crit[i] * std_errs

    coef_table <- cbind(Estimate = theta_star, intervals)
    colnames(coef_table) <- c("Estimate", prob_names)

    param_names <- c(colnames(object$data$X), paste0("psi (log-scale): ",
                                                     colnames(object$data$X2)))
    if(!is.null(object$estimates$sigma)){
      param_names <- c(param_names, "sigma")
    }
    rownames(coef_table) <- param_names

    loglik_val  <- object$estimates$loglik
    logpost_val <- object$estimates$logpost
    bayes_factors <- NULL
  }

  res <- list(
    call          = object$call,
    coefficients  = coef_table,
    bayes_factors = bayes_factors,
    loglik        = loglik_val,
    logpost       = logpost_val,
    marg_lik      = NA,
    model         = "nmr",
    method        = object$method,
    converged = if(object$method == "normal") object$estimates$converged else T
  )

  class(res) <- "summary.netreg"
  return(res)
}


#' Print method for nmr objects
#' @method print nmr
#' @export
print.nmr <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  # Safely extract the estimates column from your unified summary object
  sm <- summary(x)
  coefs <- sm$coefficients[, 1]

  # Standardize the name of the column depending on what summary() labels it
  if ("Estimate" %in% colnames(sm$coefficients)) {
    coefs <- sm$coefficients[, "Estimate"]
  } else if ("mean" %in% colnames(sm$coefficients)) {
    coefs <- sm$coefficients[, "mean"]
  }

  cat("Coefficients & Parameters:\n")
  print(round(coefs, 5))
  cat("\n")

  invisible(x)
}


#' Plot method for Network Moderated Regression models
#' @method plot nmr
#' @export
plot.nmr <- function(x, pallete = "Plasma", ...) {
  # 1. Network Geometry Set Up
  g_network  <- convert_network(as.matrix(x$data$network), out = "igraph")$A
  net_layout <- igraph::layout_with_fr(g_network)

  sm       <- summary(x)
  coef_mat <- sm$coefficients

  beta_names   <- colnames(x$data$X)
  beta_vals    <- coef_mat[beta_names, "Estimate"]
  mod_vars     <- colnames(x$data$X2)
  comm_ids     <- x$data$communities
  unique_comms <- sort(unique(comm_ids))
  n_comm       <- length(unique_comms)

  # Allocate community effects matrix (Rows = Communities, Cols = Variables)
  comm_adjustments <- matrix(0, nrow = n_comm, ncol = length(mod_vars))
  colnames(comm_adjustments) <- mod_vars
  rownames(comm_adjustments) = as.character(unique_comms)

  # 2. Extract true random effects shifts from b_matrix[community_id, variable_id]
  if (x$method == "stan") {
    draws <- x$estimates$fit

    for (m_idx in seq_along(mod_vars)) {
      v_name <- mod_vars[m_idx]
      for (c_idx in seq_along(unique_comms)) {
        c_id <- unique_comms[c_idx]

        # Target the true scale parameter block
        pat <- paste0("^b_matrix\\[", c_id, ",", m_idx, "\\]")
        col_match <- grep(pat, colnames(draws))

        if (length(col_match) > 0) {
          comm_adjustments[c_idx, v_name] <- mean(draws[, col_match])
        }
      }
    }
  } else {
    # Normal Optimization Matrix Fallback
    if (!is.null(x$estimates$b_matrix)) {
      comm_adjustments           <- x$estimates$b_matrix
      colnames(comm_adjustments) <- mod_vars
      rownames(comm_adjustments) <- as.character(unique_comms)
    }
  }

  # Genuine, built-in  palette vector indices (Colorblind friendly)
  col_steps <- 100
  viridis_colors <- hcl.colors(col_steps, palette = pallete)

  # 3. Cycle Through Moderated Slope Variables
  for (i in seq_along(mod_vars)) {
    v_name <- mod_vars[i]
    base_val <- if (v_name %in% beta_names) beta_vals[v_name] else 0

    # Map nodes to their shared community group row index
    comm_row_indices <- match(as.character(comm_ids), rownames(comm_adjustments))
    node_slopes <- base_val + comm_adjustments[comm_row_indices, v_name]

    # Map slopes onto color indices cleanly
    slope_range <- range(node_slopes, na.rm = TRUE)
    if (diff(slope_range) == 0) {
      color_mapped_vec <- rep(viridis_colors[col_steps / 2], length(node_slopes))
    } else {
      color_mapped_vec <- viridis_colors[as.numeric(cut(node_slopes, breaks = col_steps))]
    }

    # A local helper function to cleanly draw both the network and its companion scale
    # This forces R to group them together as a single visual frame unit
    draw_frame <- function() {
      plot(g_network,
           layout = net_layout,
           vertex.color = color_mapped_vec,
           vertex.label = NA,
           vertex.size = 7,
           edge.color = "lightgray",
           edge.arrow.size = 0.3,
           ...)

      title(main = paste0("Community Slopes: ", v_name), cex.main = 0.8)

      legend_ticks <- seq(slope_range[1], slope_range[2], length.out = 5)
      legend("bottomright",
             legend = round(legend_ticks, 3),
             fill = hcl.colors(5, palette = pallete),
             title = "Group Slope",
             bty = "n",
             cex = 0.8)
    }

    # Execute the initial frame render
    draw_frame()

    cat("Showing plot for variable '", v_name, "'.\n", sep = "")
    readline(prompt = "Hit ENTER for next plot...")
  }


  # 4. Final Step: Residuals vs Fitted Diagnostic Plot
  X        <- x$data$X
  Y        <- x$data$Y
  beta     <- coef_mat[beta_names, "Estimate"]
  fam_name <- x$family$family
  link_fun <- x$family$linkinv

  # Calculate fitted values exactly matching your original structure logic
  if (fam_name == "binomial") {
    # Binomial branch: Uses deviance residuals to clean up the logistic tracks
    eta_vals <- as.vector(X %*% beta)
    comm_row_indices <- match(as.character(comm_ids), rownames(comm_adjustments))
    for (m_idx in seq_along(mod_vars)) {
      v_name <- mod_vars[m_idx]
      if (v_name %in% colnames(X)) {
        eta_vals <- eta_vals + comm_adjustments[comm_row_indices, v_name] * X[, v_name]
      } else if (v_name == "(Intercept)") {
        eta_vals <- eta_vals + comm_adjustments[comm_row_indices, v_name]
      }
    }

    y_success   <- Y[, 1]
    n_trials    <- Y[, 1] + Y[, 2]
    fitted_vals <- link_fun(eta_vals)
    mu_counts   <- n_trials * fitted_vals

    t1 <- ifelse(y_success > 0, y_success * log(y_success / mu_counts), 0)
    t2 <- ifelse((n_trials - y_success) > 0, (n_trials - y_success) * log((n_trials - y_success) / (n_trials - mu_counts)), 0)
    resids <- sign(y_success - mu_counts) * sqrt(2 * (t1 + t2))

    x_axis_label <- "Fitted Probabilities"
    y_axis_label <- "Deviance Residuals"

  } else {
    # RESTORED: Exactly your original code structure for Poisson and Gaussian
    y_vec <- if (is.matrix(Y)) Y[, 1] / rowSums(Y) else Y
    fitted_vals <- as.vector(X %*% beta)

    comm_row_indices <- match(as.character(comm_ids), rownames(comm_adjustments))
    for (m_idx in seq_along(mod_vars)) {
      v_name <- mod_vars[m_idx]
      if (v_name %in% colnames(X)) {
        fitted_vals <- fitted_vals + comm_adjustments[comm_row_indices, v_name] * X[, v_name]
      } else if (v_name == "(Intercept)") {
        fitted_vals <- fitted_vals + comm_adjustments[comm_row_indices, v_name]
      }
    }

    resids       <- as.vector(y_vec - fitted_vals)
    x_axis_label <- "Fitted values"
    y_axis_label <- "Residuals"
  }

  plot(fitted_vals, resids,
       xlab = x_axis_label,
       ylab = y_axis_label,
       main = "Residuals vs Fitted",
       panel.first = grid(),
       pch = 20, col = "darkgray", ...)

  abline(h = 0, lty = 3, col = "gray40", lwd = 1.5)
  if (length(resids) > 3) {
    lines(lowess(fitted_vals, resids), col = "black", lwd = 1.5)
  }

  cat("Showing Final Diagnostics Plot.\n")
  invisible(x)
}






