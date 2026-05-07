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
#' @param method A character indicating the method to use for estimation.
#' Options for NMR include \code{"stan"}, \code{"normal"}, and \code{"IS"}.
#' @param moderated_variables A character vector of column names from the
#' design matrix to be treated as random effects (moderated by community).
#' Defaults to NULL, which treats all covariates and the intercept as moderated.
#' @param communities An optional vector of community assignments. If NULL
#' (default), communities are detected using \code{igraph::cluster_label_prop}.
#' @param test Logical. If TRUE and \code{method="stan"}, the function calculates
#' Bayes Factors for each moderated variable to test the significance of the
#' community-level moderation using bridge sampling. Defaults to TRUE.
#' @param beta_prior The prior for the fixed effects (beta). Requires list
#' or function inputs. Default is a multivariate normal with mean=0 and sigma=5.
#' If a list, must contain 'mean' and 'sigma'. If \code{method="stan"}, or
#' \code{method="IS"}, only list input is accepted.
#' @param psi_prior The prior for the log-relative scale of random effects
#' (psi). Defaults to a normal distribution with mean=0 and sd=1. Requires a
#' list with 'mean' and 'sd'. Only used when \code{method="stan"} or
#' \code{method="IS"}.
#' @param lkj_prior The shape parameter for the LKJ prior on the correlation
#' matrix of random effects. Defaults to 1 (identity/uniform). Only used
#' when \code{method="stan"} or \code{method="IS"}.
#' @param ... Additional arguments to be passed to the specific estimation
#' function (e.g., \code{iter}, \code{chains} for Stan, or \code{control}
#' for \code{optim}).
#'
#' @return A list with two classes, \code{"nam"} and \code{"nmr"}, containing:
#' @return \code{estimates} The raw estimate output (Stan fit, optim list, or IS samples).
#' @return \code{formula} The symbolic formula used.
#' @return \code{model} String identifying the model type ("nmr").
#' @return \code{method} String identifying the estimation method.
#' @return \code{data} List containing the processed Y, X, Z matrices and community vector.
#' @return \code{call} The matched call.
#'
#' @import mvtnorm
#' @import rstan
#' @import igraph
#' @import lme4
#' @importFrom bridgesampling bridge_sampler bf
#' @importFrom stats model.matrix model.frame as.formula
#' @export



# Model-specific wrapper function for Network Moderated Regression
nmr <- function(formula, data, network, method = "stan",
                moderated_variables = NULL, communities = NULL,
                beta_prior = list(mean=0, sigma=5),
                psi_prior = list(mean=0, sd=1),
                lkj_prior = 1, test=TRUE, ...){

  X <- model.matrix(formula, data = data)
  Y <- model.frame(formula, data = data)[, 1]

  # Community logic
  if(is.null(communities)){
    A <- row_norm(row_norm(network), reverse = TRUE)
    g <- convert_adjacency(A, out="igraph")$A
    communities <- igraph::cluster_label_prop(g)$membership
  }

  # Moderator logic
  if(is.null(moderated_variables)) {
    Z_1 <- X
  } else {
    Z_1 <- X[, moderated_variables, drop = FALSE]
  }

  nmrdat <- list(Y=Y, X=X, Z_1=Z_1, communities=communities,
                 beta_prior=beta_prior, psi_prior=psi_prior,
                 lkj_prior=lkj_prior, p=ncol(X), m=ncol(Z_1),
                 test=test, ...)

  if(tolower(method) == "stan"){
    nmrout <- nmr_stan(nmrdat)
  }else if(tolower(method) == "normal"){
    nmrout <- nmr_normal(nmrdat)
  }else if(tolower(method) %in% c("importance", "is")){
    nmrout <- nmr_IS(nmrdat)
  }

  final_out <- list(
    estimates = nmrout,
    formula = formula,
    model = "nmr",
    method = tolower(method),
    data = list(Y = Y, X = X, Z = Z_1, communities = communities),
    call = match.call()
  )
  class(final_out) <- c("nam", "nmr")
  return(final_out)
}

nmr_likelihood <- function(theta, Y, X, Z_1, communities, log = TRUE) {
  # 1. Parameter Extraction
  # theta = [beta (K), psi (M_1), L_elements, random_effects (N_1*M_1), sigma]
  K <- ncol(X)
  M_1 <- ncol(Z_1)
  N_1 <- length(unique(communities))

  # Fixed effects
  beta <- theta[1:K]

  # Log-relative scales
  psi <- theta[(K + 1):(K + M_1)]
  lambda <- exp(psi)

  # Global scale (sigma)
  sigma <- theta[length(theta)]

  # Random effects (r_1)
  # We assume r_1 is passed as a flattened vector of the N_1 x M_1 matrix
  r_1_vec <- theta[(K + M_1 + 1):(K + M_1 + (N_1 * M_1))]
  r_1 <- matrix(r_1_vec, nrow = N_1, ncol = M_1)

  # 2. Compute Mean Structure
  mu <- as.vector(X %*% beta)
  J_1 <- as.integer(as.factor(communities))
  for (i in 1:length(Y)) {
    mu[i] <- mu[i] + sum(r_1[J_1[i], ] * Z_1[i, ])
  }

  # 3. Calculate Normal Density
  out <- dnorm(Y, mean = mu, sd = sigma, log = log)

  if (log) {
    return(sum(out))
  } else {
    return(prod(out))
  }
}

nmr_posterior <- function(theta, Y, X, Z_1, communities, log = TRUE,
                          bprior, psiprior, sigmaprior){
  K <- ncol(X)
  M_1 <- ncol(Z_1)
  groups <- as.factor(communities)
  N_1 <- length(levels(groups))

  # theta order: beta (K), psi (M_1), r_1 (N_1*M_1), sigma (1)
  beta <- theta[1:K]
  psi  <- theta[(K + 1):(K + M_1)]
  r_1_vec <- theta[(K + M_1 + 1):(K + M_1 + (N_1 * M_1))]
  r_1  <- matrix(r_1_vec, nrow = N_1, ncol = M_1)
  sigma <- theta[length(theta)]

  # 1. Likelihood: Y ~ N(mu, sigma)
  log_lik <- nmr_likelihood(theta, Y, X, Z_1, communities, log = TRUE)

  # 2. Hierarchical Penalty (Independent Random Effects)
  # Each column j in r_1 corresponds to moderated variable j with scale sigma * exp(psi[j])
  log_prior_re <- 0
  for(j in 1:M_1){
    col_sd <- sigma * exp(psi[j])
    log_prior_re <- log_prior_re + sum(dnorm(r_1[,j], 0, col_sd, log = TRUE))
  }

  # 3. Parameter Priors
  log_prior_beta  <- bprior(beta, log = TRUE)
  log_prior_psi   <- psiprior(psi, log = TRUE)
  log_prior_sigma <- sigmaprior(sigma, log = TRUE)

  total <- log_lik + log_prior_re + log_prior_beta + log_prior_psi + log_prior_sigma

  return(if(log) total else exp(total))
}


# Normal approximation estimation for NMR
nmr_normal <- function(nmrdat) {
  # 1. Setup Priors
  K <- nmrdat$p
  M_1 <- nmrdat$m
  groups <- as.factor(nmrdat$communities)
  N_1 <- length(levels(groups))

  bprior <- function(beta, log = TRUE) {
    mvtnorm::dmvnorm(beta, rep(nmrdat$beta_prior$mean, K), diag(K) * nmrdat$beta_prior$sigma, log = log)
  }
  psiprior <- function(psi, log = TRUE) {
    sum(dnorm(psi, nmrdat$psi_prior$mean, nmrdat$psi_prior$sd, log = log))
  }
  sigmaprior <- function(sigma, log = TRUE) {
    dt(sigma / 2.5, df = 3, log = log)
  }

  # 2. Initialize with lme4
  # Construct a formula: y ~ X + (Z | groups)
  lmer_df <- data.frame(Y = nmrdat$Y, G = groups)
  for(i in 1:K) lmer_df[[paste0("X", i)]] <- nmrdat$X[,i]
  for(i in 1:M_1) lmer_df[[paste0("Z", i)]] <- nmrdat$Z_1[,i]

  x_vars <- paste0("X", 1:K, collapse = " + ")
  z_vars <- paste0("(0 + Z", 1:M_1, " | G)", collapse = " + ")
  lmer_form <- as.formula(paste0("Y ~ 0 + ", x_vars, " + ", z_vars))

  lmer_fit <- lme4::lmer(lmer_form, data = lmer_df)

  # Extract initial values
  init_theta <- c(
    lme4::fixef(lmer_fit),
    log(as.data.frame(lme4::VarCorr(lmer_fit))$sdcor / sigma(lmer_fit)),
    as.vector(as.matrix(lme4::ranef(lmer_fit)$G)),
    sigma(lmer_fit)
  )

  # Optimize
  myopt <- optim(par = init_theta,
                 fn = nmr_posterior,
                 method = "BFGS",
                 control = list(fnscale = -1),
                 hessian = TRUE,
                 Y = nmrdat$Y,
                 X = nmrdat$X,
                 Z_1 = nmrdat$Z_1,
                 communities = nmrdat$communities,
                 bprior = bprior,
                 psiprior = psiprior,
                 sigmaprior = sigmaprior)

  # 5. Format Output
  beta_ests <- myopt$par[1:K]
  names(beta_ests) <- colnames(nmrdat$X)

  list(
    beta = beta_ests,
    psi = myopt$par[(K+1):(K+M_1)],
    r_1 = matrix(myopt$par[(K+M_1+1):(K+M_1+(N_1*M_1))], N_1, M_1),
    sigma = myopt$par[length(myopt$par)],
    loglik = nmr_likelihood(myopt$par, nmrdat$Y, nmrdat$X, nmrdat$Z_1, nmrdat$communities),
    converged = myopt$convergence == 0,
    raw_opt = myopt
  )
}

# Importance Sampler function for NMR
nmr_IS <- function(nmrdat){
  if(is.null(nmrdat$samples)){
    nmrdat$samples <- 5000
  }
  S <- nmrdat$samples

  # 1. Get the Proposal Distribution (using Normal method)
  fit_norm <- nmr_normal(nmrdat)
  mu <- fit_norm$raw_opt$par
  Sigma <- solve(-fit_norm$raw_opt$hessian)

  # 2. Sample from the Multivariate Normal Proposal
  cat("Generating", S, "samples from the proposal distribution...\n")
  theta_samples <- mvtnorm::rmvnorm(S, mean = mu, sigma = Sigma)

  # 3. Calculate Importance Weights
  log_weights <- numeric(S)

  # Re-setup priors for the weight calculation
  K <- nmrdat$p
  bprior <- function(beta){
    mvtnorm::dmvnorm(beta, rep(nmrdat$beta_prior$mean, K), diag(K) * nmrdat$beta_prior$sigma, log = TRUE)
  }
  psiprior <- function(psi){
    sum(dnorm(psi, nmrdat$psi_prior$mean, nmrdat$psi_prior$sd, log = TRUE))
  }
  sigmaprior <- function(sigma){
    dt(sigma / 2.5, df = 3, log = TRUE)
  }

  cat("Calculating weights...\n")
  for (i in 1:S) {
    # Log-Posterior (Target)
    log_post <- nmr_posterior(theta_samples[i,], nmrdat$Y, nmrdat$X, nmrdat$Z_1,
                              nmrdat$communities, log = TRUE,
                              bprior, psiprior, sigmaprior)

    # Log-Proposal Density
    log_prop <- mvtnorm::dmvnorm(theta_samples[i,], mean = mu,
                                 sigma = Sigma, log = TRUE)

    log_weights[i] <- log_post - log_prop
  }

  # 4. Normalize Weights
  max_log_w <- max(log_weights)
  weights <- exp(log_weights - max_log_w)
  weights <- weights / sum(weights)

  # 5. Resample to get the "Posterior"
  indices <- sample(1:S, size = S, replace = TRUE, prob = weights)
  resampled_theta <- theta_samples[indices, ]

  # 6. Format Output
  M_1 <- nmrdat$m
  N_1 <- length(unique(nmrdat$communities))

  # Extract means from resampled distribution
  beta_means <- colMeans(resampled_theta[, 1:K, drop = FALSE])
  names(beta_means) <- colnames(nmrdat$X)

  list(
    beta = beta_means,
    psi = colMeans(resampled_theta[, (K+1):(K+M_1), drop = FALSE]),
    sigma = mean(resampled_theta[, ncol(resampled_theta)]),
    samples = list(
      beta = resampled_theta[, 1:K],
      psi = resampled_theta[, (K+1):(K+M_1)],
      r_1 = resampled_theta[, (K+M_1+1):(K+M_1+(N_1*M_1))],
      sigma = resampled_theta[, ncol(resampled_theta)]
    ),
    log_marginal_lik = max_log_w + log(mean(exp(log_weights - max_log_w))),
    converged = TRUE
  )
}


# Stan Estimation function for NMR
nmr_stan <- function(nmrdat){
  test <- nmrdat$test
  if(!requireNamespace("rstan", quietly = TRUE)) stop("rstan required.")
  if(test && !requireNamespace("bridgesampling", quietly = TRUE)) stop("bridgesampling required.")

  # Standard setup
  K <- nmrdat$p
  if(length(nmrdat$beta_prior$mean) == 1){
    b_mean <- rep(nmrdat$beta_prior$mean, K)
  }else{
    b_mean <- nmrdat$beta_prior$mean
  }
  if(length(nmrdat$beta_prior$sigma) == 1){
    b_sigma <- diag(K) * nmrdat$beta_prior$sigma
  }else{
    nmrdat$beta_prior$sigma
  }

  # Helper to run sampling and bridge
  fit_and_bridge <- function(Z){
    s_data <- list(
      N = length(nmrdat$Y), Y = nmrdat$Y, K = K, X = nmrdat$X,
      N_1 = length(unique(nmrdat$communities)), M_1 = ncol(Z),
      J_1 = as.integer(as.factor(nmrdat$communities)), Z_1 = Z,
      b_mean = b_mean, b_sigma = b_sigma,
      psi_mean = nmrdat$psi_prior$mean, psi_sd = nmrdat$psi_prior$sd,
      lkj_shape = nmrdat$lkj_prior
    )
    # Higher iterations for bridge sampling stability
    if(test){
      iters <- 10000
    }else{
      iters <- (if(is.null(nmrdat$iterations)) 2000 else nmrdat$iterations)
    }

    fit <- rstan::sampling(stanmodels$nmr_model, data = s_data, iter = iters, ...)
    if(test){
      lml <- bridgesampling::bridge_sampler(fit, silent = TRUE)
    }else{
      lml <- NULL
    }
    return(list(fit = fit, lml = lml))
  }

  # 1. Run Full Model
  full_res <- fit_and_bridge(nmrdat$Z_1)

  bfs <- list()
  if(test){
    # 2. Run Reduced Models for each moderated variable
    Z_full <- nmrdat$Z_1
    M_1 <- ncol(Z_full)
    col_names <- colnames(Z_full)

    for(i in 1:M_1) {
      cat(sprintf("\nTesting moderation of: %s...\n", col_names[i]))
      Z_reduced <- Z_full[, -i, drop = FALSE]
      red_res <- fit_and_bridge(Z_reduced)

      bfs[[col_names[i]]] <- bridgesampling::bf(full_res$lml, red_res$lml)$bf
    }
  }

  return(list(fit = full_res$fit, bfs = bfs))
}




#' Summary method for nmr objects
#' @param object An object of secondary class "nmr"
#' @param probs A numeric vector of probabilities for the interval estimates.
#' @param ... Additional arguments
#' @export
summary.nmr <- function(object, probs = c(0.025, 0.975), ...){
  p <- ncol(object$data$X)
  m <- ncol(object$data$Z)
  prob_names <- paste0(probs * 100, "%")

  if(object$method %in% c("importance", "is")){
    samps <- object$estimates$samples
    # Combine beta, psi, and sigma for a unified table
    theta_samples <- cbind(samps$beta, samps$psi, samps$sigma)

    coef_table <- t(apply(theta_samples, 2, function(x){
      c(Estimate = mean(x), quantile(x, probs = probs))
    }))

    marg_lik <- object$estimates$log_marginal_lik
    loglik_val <- NA # Calculated at the mean if needed

  }else if(object$method == "normal"){
    # theta: beta (p), psi (m), r_1 (n_comm * m), sigma (1)
    theta_star <- object$estimates$raw_opt$par

    # Standard Errors from the Hessian
    hinv <- tryCatch(solve(-object$estimates$raw_opt$hessian),
                     error = function(e) matrix(NA, length(theta_star), length(theta_star)))
    std_errs <- sqrt(diag(hinv))

    # Calculate intervals
    z_crit <- qnorm(probs)
    intervals <- matrix(NA, nrow = length(theta_star), ncol = length(probs))
    for(i in 1:length(probs)) {
      intervals[,i] <- theta_star + z_crit[i] * std_errs
    }
    idx_to_show <- c(1:(p + m), length(theta_star))
    coef_table <- cbind(Estimate = theta_star[idx_to_show], intervals[idx_to_show, ])

    marg_lik <- NA
    loglik_val <- object$estimates$loglik
  }else if(object$method == "stan"){
    # Extract results from list structure
    fit <- object$estimates$fit
    bfs <- object$estimates$bfs

    s_mat <- rstan::summary(fit, probs = probs)$summary

    # Extract Beta, Psi, and Sigma
    beta_idx <- grep("^b\\[", rownames(s_mat))
    psi_idx <- grep("^psi\\[", rownames(s_mat))
    sigma_idx <- which(rownames(s_mat) == "sigma")

    coef_table <- s_mat[c(beta_idx, psi_idx, sigma_idx), c("mean", prob_names)]

    # Handle BF Column
    if(length(bfs) > 0){
      bf_col <- rep(NA, nrow(coef_table))
      psi_row_start <- p + 1
      for(i in 1:m) {
        var_name <- colnames(object$data$Z)[i]
        bf_col[psi_row_start + i - 1] <- bfs[[var_name]]
      }
      coef_table <- cbind(coef_table, BF = bf_col)
    }

    # Naming
    var_names <- c(colnames(object$data$X),
                   paste0("log_rel_sd_", colnames(object$data$Z)),
                   "sigma")
    rownames(coef_table) <- var_names

  }

  # Label the rows correctly
  var_names <- c(colnames(object$data$X),
                 paste0("log_rel_sd_", colnames(object$data$Z)),
                 "sigma")
  rownames(coef_table) <- var_names
  colnames(coef_table) <- c("Estimate", prob_names)

  res <- list(
    call = object$call,
    coefficients = coef_table,
    loglik = loglik_val,
    marg_lik = marg_lik,
    model = object$model,
    method = object$method,
    n_communities = length(unique(object$data$communities)),
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

  # --- New: NMR Bayes Factor Interpretation ---
  if(x$model == "nmr" && "BF" %in% colnames(x$coefficients)){
    cat("\nModeration Significance (Kass & Raftery, 1995):\n")

    # Identify the rows corresponding to log_rel_sd (moderation)
    bfs <- x$coefficients[grep("log_rel_sd_", rownames(x$coefficients)), "BF"]
    names(bfs) <- gsub("log_rel_sd_", "", names(bfs))

    for(var in names(bfs)) {
      val <- bfs[var]
      if(is.na(val)) next

      if(val > 100){
        msg <- "Decisive evidence for community moderation"
      }else if(val > 10){
        msg <- "Strong evidence for community moderation"
      }else if(val > 3){
        msg <- "Positive evidence for community moderation"
      }else if(val > 1/3){
        msg <- "Negligible/Inconclusive evidence"
      }else if(val > 1/10){
        msg <- "Positive evidence AGAINST moderation (treat as fixed)"
      }else{
        msg <- "Strong evidence AGAINST moderation (treat as fixed)"
      }
      cat(sprintf("  - %-15s : BF = %-8.2f [%s]\n", var, val, msg))
    }
  }

  cat("\nModel Fit Statistics:\n")
  cat("Log-Likelihood: ", round(x$loglik, 4))

  if(!is.null(x$logpost) && !is.na(x$logpost)){
    cat("\nLog-Posterior:  ", round(x$logpost, 4))
  }

  if(!is.null(x$marg_lik) && !is.na(x$marg_lik)){
    cat("\nLog-Marginal-Lik:", round(x$marg_lik, 4))
  }

  if(!is.null(x$converged) && isFALSE(x$converged)){
    cat("\n\n*** WARNING: Optimization did not converge! ***\n")
  }

  cat("\n")
  invisible(x)
}


