#' nam - General Wrapper for Network Analysis Models
#'
#' A wrapper function for performing various models within the
#' network autocorrelation family.
#'
#' @param formula An object of class \code{\link[stats]{formula}}: a symbolic
#'   description of the model to be fitted.
#' @param data An optional data frame, list, or environment. If \code{data} is an
#'   \code{igraph} object, the variables and adjacency matrix are extracted
#'   automatically (ignoring the \code{network} argument).
#' @param network An adjacency matrix or network object. Supported classes include:
#'   \code{matrix}, \code{dgCMatrix}, \code{igraph}, \code{network}, or an edgelist.
#' @param model A character string specifying the model type. Options include:
#'   \code{"effects"} (NEM), \code{"disturbances"} (NDM), \code{"nmr"}, \code{"iv"},
#'   \code{"hane"}, and \code{"hand"}. See \code{?models.nam} for details and citations.
#' @param method A character string specifying the estimation method. Options
#'   are model specific but can potentially include \code{"normal"}
#'   (Maximum Likelihood/Optimization), \code{"stan"}
#'   (MCMC via Hamiltonian Monte Carlo), and \code{"lnam"}.
#' @param rownorm Logical. If \code{TRUE} (default), the adjacency matrix will be
#'   row-normalized. This may be ignored by models requiring specific normalization.
#' @param ... Method-specific arguments passed to underlying estimation functions
#'   (e.g., \code{rstan::sampling} or \code{stats::optim}).
#'
#' @return An object of class \code{c("nam", "subclass")}, containing:
#' \describe{
#'   \item{estimates}{A list of parameter estimates (beta, s2, rho) and estimation metadata.}
#'   \item{formula}{The formula used in the model.}
#'   \item{model}{The string identifier for the model type.}
#'   \item{method}{The string identifier for the estimation method.}
#'   \item{data}{A list containing the response, predictors, and adjacency matrix used.}
#'   \item{priors}{The prior settings used (for Bayesian methods).}
#'   \item{call}{The matched call.}
#' }
#'
#' @seealso \code{\link{summary.nam}}, \code{\link{print.nam}}
#' @export


nam <- function(formula, data=NULL, network, model="effects", method="normal", rownorm=TRUE, ...){

  # 1. Uniform Data Handling
  # Convert network to matrix early so we don't have to do it inside every if-block
  if(inherits(data, "igraph")){
    df <- convert_adjacency(data, "matrix")
    data <- df$df
    network <- df$A
  }else{
    network <- convert_adjacency(network, "matrix")$A
  }

  # 2. Dispatch to specific model functions
  # We use tolower() once to simplify the logic
  mod <- tolower(model)

  if(mod %in% c("effects", "nem")){
    out <- nem(formula=formula, data=data, network=network, method=method, rownorm=rownorm, ...)

  }else if(mod %in% c("disturbances", "ndm")){
    out <- ndm(formula=formula, data=data, network=network, method=method, rownorm=rownorm, ...)

  }else if(mod %in% c("iv")){
    out <- iv(formula=formula, data=data, network=network, method=method, rownorm=rownorm, ...)

  }else if(mod %in% c("nmr")){
    out <- nmr(formula=formula, data=data, network=network, method=method, rownorm=rownorm, ...)

  }else if(mod %in% c("hane", "homopholy-adjusted network effects")){
    out <- hane(formula=formula, data=data, network=network, method=method, rownorm=rownorm, ...)

  }else if(mod %in% c("hand", "homopholy-adjusted network disturbances")){
    out <- hand(formula=formula, data=data, network=network, method=method, rownorm=rownorm, ...)

  }else{
    stop("You did not provide a valid model name. Options are: 'nem', 'ndm', 'iv', 'nmr', or 'ha'.")
  }

  return(out)
}


#' Print method for nam objects
#' @param x An object of class "nam"
#' @param ... Additional arguments (ignored)
#' @export
print.nam <- function(x, ...) {
  cat("\nNetwork Analysis Model (NAM)\n")
  cat("----------------------------\n")
  cat("Model Type: ", x$model, "\n")
  cat("Method:     ", x$method, "\n")
  cat("\nCall:\n")
  print(x$call)

  cat("\nCoefficients (Beta):\n")
  # Pulling from your nested 'estimates' list
  if (!is.null(x$estimates$beta)) {
    print(round(x$estimates$beta, 4))
  } else {
    cat("Estimates not found in output.\n")
  }

  cat("\nNetwork Correlation (Rho):", round(x$estimates$rho, 4), "\n")
  cat("Residual Variance (S2):   ", round(x$estimates$s2, 4), "\n")

  cat("\n---\n")
  cat("Use summary() for detailed statistics and p-values.\n")

  invisible(x)
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



