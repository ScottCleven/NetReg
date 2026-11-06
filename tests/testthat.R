# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(NetReg)
library(igraph)

test_check("NetReg")


sim_netdat <- function(n = 50, A = NULL, model = "effects", beta = c(1, 2, 3),
                       formula=NULL, moderated_formula = NULL, ...) {

  # 1. Base network generation (if not provided)
  if(is.null(A)){
    n <- 50
    g <- erdos.renyi.game(n, p = 0.15, directed = FALSE)
    V(g)$name <- paste0("Node_", 1:n)

    # 2. Add vertex-level covariates and simulated response outcomes
    V(g)$x1 <- runif(n, -2, 2)
    V(g)$x2 <- rnorm(n, mean = 0, sd = 1)
    V(g)$w1 <- rnorm(n, mean = 1, sd = 0.5)
    conv <- convert_adjacency(g)
    A <- conv$A
    mydat <- conv$df
  }

  # 2. Row-normalize and clean up original A for memory saving
  A_norm <- row_norm(A)
  I_n <- diag(n)

  p <- length(beta)
  X <- model.matrix(as.formula(formula), data = mydat)

  # 5. Extract core dot-dot-dot arguments with logical defaults
  args  <- list(...)
  if(!is.null(args$sigma)){
    sigma <- args$sigma
  }else{
    sigma <- 1
  }
  epsilon <- rnorm(n, mean = 0, sd = sigma)

  # --- Network Effects Model ---
  if(model == "effects"){
    if(!is.null(args$rho)){
      rho <- args$rho
    }else{
      rho <- 0.36
    }

    L <- solve(I_n - rho * A_norm)
    y <- as.vector(L %*% (X %*% beta + epsilon))

    return(list(
      data        = cbind(y = y, pool_df),
      A_norm      = A_norm,
      true_params = list(beta = setNames(beta, colnames(X)), rho = rho, sigma = sigma)
    ))
  }

  # --- Network Disturbances Model ---
  if(model == "disturbances"){
    if(!is.null(args$rho)){
      rho <- args$rho
    }else{
      rho <- 0.36
    }

    L <- solve(I_n - rho * A_norm)
    nu <- L %*% epsilon
    y <- as.vector((X %*% beta) + nu)

    return(list(
      data        = cbind(y = y, mydat),
      A_norm      = A_norm,
      true_params = list(beta = setNames(beta, colnames(X)), rho = rho, sigma = sigma)
    ))
  }

  # --- Heterogeneous Susceptibilities Matrix Setup (W) ---
  if(model %in% c("hs_effects", "hs_disturbances")){
    if(is.null(moderated_formula)){
      W <- X
    }else{
      W <- model.matrix(moderated_formula, data = mydat)
    }
    if(!is.null(args$gamma)){
      gamma_hs <- args$gamma
    }else{
      gamma_hs <- rep(0.36, ncol(W))
    }
    rho_vec <- as.vector(W %*% gamma_hs)
  }


  # --- HS Network Effects Model ---
  if(model == "hs_effects"){
    L <- solve(I_n - diag(rho_vec) %*% A_norm)
    y <- as.vector(L %*% (X %*% beta + epsilon))

    return(list(
      data        = cbind(y = y, mydat),
      A_norm      = A_norm,
      true_params = list(beta = setNames(beta, colnames(X)), gamma = setNames(gamma, colnames(W)), rho_vec = rho_vec, sigma = sigma)
    ))
  }

  # --- HS Network Disturbances Model ---
  if (model == "hs_disturbances") {
    L <- solve(I_n - diag(rho_vec) %*% A_norm)
    nu <- L %*% epsilon
    y <- as.vector((X %*% beta) + nu)

    return(list(
      data        = cbind(y = y, mydat),
      A_norm      = A_norm,
      true_params = list(beta = setNames(beta, colnames(X)), gamma = setNames(gamma, colnames(W)), rho_vec = rho_vec, sigma = sigma)
    ))
  }








  # --- Homophily-Adjusted Network Effects Model ("hane") ---
  if (model == "hane"){
    # Extract model-specific arguments with defaults
    args  <- list(...)
    if(!is.null(args$rho)){
      rho <- args$rho
    }else{
      rho <- 0.36
    }
    if(!is.null(args$gamma)){
      gamma <- args$gamma
    }else{
      gamma <- 0.5
    }
    if(!is.null(args$D)){
      D <- args$D
    }else{
      2
    }

    # Generate latent coordinate matrix U and its row-normalized distance matrix
    U <- matrix(rnorm(n * D), nrow = n, ncol = D)
    Dist_U <- as.matrix(dist(U, method = "euclidean"))

    # Solve structural equation:
    # y = rho * A_norm %*% y + Dist_U %*% gamma + X %*% beta + epsilon
    L <- solve(I_n - rho * A_norm)
    y <- as.vector(inv_env %*% (X %*% beta + epsilon))

    return(list(
      data   = cbind(y = y, pool_df),
      A_norm = A_norm,
      true_params = list(
        beta   = setNames(beta, colnames(X)),
        rho    = rho,
        gamma  = gamma,
        sigma  = sigma,
        latent_U = U
      )
    ))
  }


  # --- Homophily-Adjusted Network Disturbances Model ("hand") ---
  if (model == "hand") {
    # Extract model-specific arguments with defaults
    args  <- list(...)
    rho   <- if (!is.null(args$rho)) args$rho else 0.2
    gamma <- if (!is.null(args$gamma)) args$gamma else 0.15
    D     <- if (!is.null(args$D)) args$D else 2

    # Generate latent coordinate matrix U and its row-normalized distance matrix
    U <- matrix(rnorm(n * D), nrow = n, ncol = D)
    Dist_U <- as.matrix(dist(U, method = "euclidean"))
    Dist_U_norm <- row_norm(Dist_U)

    # Solve for homophily-structured disturbances:
    # nu = (I - rho * A_norm + gamma * Dist_U_norm)^(-1) %*% epsilon
    inv_env <- solve(I_n - rho * A_norm + gamma * Dist_U_norm)
    nu <- inv_env %*% epsilon

    # Compute final outcomes
    y <- as.vector((X %*% beta) + nu)

    return(list(
      data   = cbind(y = y, pool_df),
      A_norm = A_norm,
      true_params = list(
        beta   = setNames(beta, colnames(X)),
        rho    = rho,
        gamma  = gamma,
        sigma  = sigma,
        latent_U = U
      )
    ))
  }


  # --- Network Moderated Regression Model ("nmr") ---
  if (model == "nmr") {
    # 1. Extract model-specific arguments with defaults
    args     <- list(...)
    k_groups <- if (!is.null(args$k_groups)) args$k_groups else 4

    # 2. Recover the non-row-normalized version of A for community detection
    A_raw <- row_norm(A_norm, reverse = TRUE)
    km_res <- kmeans(A_raw, centers = min(k_groups, n - 1), nstart = 10)
    communities <- as.factor(km_res$cluster)
    num_communities <- length(levels(communities))
    rm(A_raw) # Clean up to save memory

    # 3. Build the Moderated Effects Design Matrix (Z_mod) via formula
    # If no formula is specified, default to moderating all variables present in X
    if (is.null(moderated_formula)) {
      Z_mod <- X
    } else {
      Z_mod <- model.matrix(moderated_formula, data = pool_df)
    }

    # Match the number of variance scale parameters (lambda) to chosen moderated variables
    q_mod <- ncol(Z_mod)
    lambda_diag <- if (!is.null(args$lambda)) args$lambda else rep(0.3, q_mod)

    # 4. Generate random effects coefficients matrix: (num_communities x q_mod)
    b_mat <- matrix(0, nrow = num_communities, ncol = q_mod)
    for (j in 1:q_mod) {
      b_mat[, j] <- rnorm(num_communities, mean = 0, sd = lambda_diag[j])
    }
    colnames(b_mat) <- colnames(Z_mod)

    # 5. Calculate the structural Khatri-Rao product contribution (Z %*% b)
    Z_b_contribution <- rep(0, n)
    for (i in 1:n) {
      comm_id <- as.numeric(communities[i])
      # Row-by-row dot product of a node's moderated traits and its cluster's random slopes
      Z_b_contribution[i] <- sum(Z_mod[i, ] * b_mat[comm_id, ])
    }

    # 6. Calculate outcome variable: y = X %*% beta + Z_mod %*% b + epsilon
    y <- as.vector(X %*% beta + Z_b_contribution + epsilon)

    # 7. Package outcomes, predictors, and community IDs into the dataframe
    df <- cbind(y = y, pool_df, community = communities)

    return(list(
      data   = df,
      A_norm = A_norm,
      true_params = list(
        beta        = setNames(beta, colnames(X)),
        lambda      = setNames(lambda_diag, colnames(Z_mod)),
        communities = communities,
        random_b    = b_mat
      )
    ))
  }





  stop("Invalid model specification.")
}














