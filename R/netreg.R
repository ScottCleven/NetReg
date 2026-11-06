#' netreg - General Wrapper for Network Analysis Models
#'
#' A unified wrapper function for executing various modeling frameworks within the
#' network autocorrelation and network regression family.
#'
#' @param formula An object of class \code{\link[stats]{formula}}: a symbolic
#'   description of the model to be fitted.
#' @param data An optional data frame, list, or environment containing model variables.
#' @param network An adjacency matrix or graph structure. Supported classes include:
#'   \code{matrix}, \code{dgCMatrix}, \code{igraph}, \code{network}, or an edgelist.
#' @param model A character string specifying the model type. Options include:
#'   \code{"effects"} (NEM), \code{"disturbances"} (NDM), \code{"ego_effects"},
#'   \code{"ego_disturbances"}, \code{"hs_effects"}, \code{"hs_disturbances"},
#'   \code{"hane"}, \code{"hand"}, and \code{"nmr"}.
#' @param method A character string specifying the estimation method (e.g., \code{"normal"}
#'   or \code{"stan"}). Options are model specific.
#' @param rownorm Logical. If \code{TRUE} (default), the adjacency matrix will be
#'   row-normalized before fitting.
#' @param ... Method-specific arguments passed to downstream subroutines.
#'
#' @return An object of class \code{c("[model]", "netreg")}.
#' @export
netreg <- function(formula, data = NULL, network, model = "effects",
                   method = "normal", rownorm = TRUE, ...){

  mod <- tolower(model)

  # Route directly to your polished, specific wrapper engines
  if(mod %in% c("effects", "nem")){
    out <- nem(formula=formula, data=data, network=network,
               method=method, rownorm=rownorm, ...)

  }else if(mod %in% c("disturbances", "ndm")){
    out <- ndm(formula=formula, data=data, network=network,
               method=method, rownorm=rownorm, ...)

  }else if(mod %in% c("hs_nem", "hs_effects",
                      "heterogeneous_effects")){
    out <- hs_effects(formula=formula, data=data, network=network,
                      method=method, rownorm=rownorm, ...)

  }else if(mod %in% c("hs_ndm", "hs_disturbances",
                      "heterogeneous_disturbances")){
    out <- hs_disturbances(formula=formula, data=data, network=network,
                           method=method, rownorm=rownorm, ...)

  }else if(mod %in% c("hane", "homophily_effects")){
    out <- hane(formula=formula, data=data, network=network,
                method=method, rownorm=rownorm, ...)

  }else if(mod %in% c("hand", "homophily_disturbances")){
    out <- hand(formula=formula, data=data, network=network,
                method=method, rownorm=rownorm, ...)

  }else if(mod %in% c("nmr", "network_moderated")){
    out <- nmr(formula=formula, data=data, network=network,
               method=method, rownorm=rownorm, ...)

  }else{
    stop(paste("Model token '", model, "' unrecognized. Options include: 'nem',
               'ndm', 'hs_effects', 'hs_disturbances', 'hane', 'hand', or
               'nmr'."))
  }

  return(out)
}

#' Base Print method for netreg objects
#' @method print netreg
#' @export
print.netreg <- function(x, ...) {
  cat("\nNetwork Regression Modeling Framework (NetReg)\n")
  cat("----------------------------------------------\n")
  cat("Model Type: ", x$model, "\n")
  cat("Method:     ", x$method, "\n")
  cat("\nCall:\n")
  print(x$call)

  cat("\nCoefficients (Beta):\n")
  if (!is.null(x$estimates$beta)) {
    print(round(x$estimates$beta, 4))
  } else {
    cat("Estimates not found.\n")
  }

  # Safely handle models without static global rho/s2 values
  if (!is.null(x$estimates$rho)) {
    if(length(x$estimates$rho) == 1) {
      cat("\nNetwork Dependence (Rho):", round(x$estimates$rho, 4), "\n")
    } else {
      cat("\nNetwork Dependence (Rho): [Heterogeneous Vector Estimated]\n")
    }
  }

  if (!is.null(x$estimates$s2)) {
    cat("Residual Variance (S2):   ", round(x$estimates$s2, 4), "\n")
  } else if (!is.null(x$estimates$sigma)) {
    cat("Residual Scale (Sigma):   ", round(x$estimates$sigma, 4), "\n")
  }

  cat("\n---\n")
  cat("Use summary() for intervals, model diagnostics, and p-values.\n")

  invisible(x)
}

#' Print method for netreg summaries
#' @method print summary.netreg
#' @export
print.summary.netreg <- function(x, ...) {
  cat("\nSummary of Network Analysis Model\n")
  cat("Model Type:", x$model, "| Estimation Method:", x$method, "\n")

  cat("\nCall:\n")
  print(x$call)

  # 1. Print heterogeneous individual rho summaries if present
  if(!is.null(x$rho_summary)) {
    cat("\nDistribution of Individual Susceptibilities (rho):\n")
    print(round(x$rho_summary, 4))
  }

  # 2. Print consolidated coefficient interval frames
  cat("\nParameter Estimates and Intervals:\n")
  print(round(x$coefficients, 4))

  # 3. Print bridge-sampled Bayes Factors for NMR model components
  if(!is.null(x$bayes_factors)) {
    cat("\nBayes Factors for Moderation Significance Tests (H1:
        Moderated vs H0: Constant):\n")
    print(round(x$bayes_factors, 4))
  }

  cat("\nModel Fit Statistics:\n")
  if(!is.null(x$loglik) && !is.na(x$loglik)) {
    cat("Log-Likelihood: ", round(x$loglik, 4))
  }

  if(!is.null(x$logpost) && !is.na(x$logpost)) {
    cat("\nLog-Posterior:  ", round(x$logpost, 4))
  }

  # Display Marginal Likelihood from Importance Samplers
  if(!is.null(x$marg_lik) && !is.na(x$marg_lik)) {
    cat("\nLog-Marginal-Lik:", round(x$marg_lik, 4))
  }

  if(!is.null(x$converged) && isFALSE(x$converged)) {
    cat("\n\n*** WARNING: Sampler/Optimization did not reach full convergence!
        ***\n")
  }

  cat("\n")
  invisible(x)
}




