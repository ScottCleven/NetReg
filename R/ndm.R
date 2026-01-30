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
#' Options vary by model but all potential options are c("normal", 
#' "greta", "stan", "lnam"). Run ?methods.nam() for model-specific designations.
#' @param beta_prior The prior for the beta vector of the model. Requires list 
#' or function inputs. Default is a multivariate normal with mean=0 and sigma=1 
#' If you supply a list, you must have one element named 'mean' which is the 
#' mean of the multivariate normal and one named 'sigma' which is the 
#' variance-covariance matrix. If mean is a single integer, 
#' that mean will be used for every element in the mean vector and if sigma is 
#' a single integer, it will multiply that integer by the identity matrix. 
#' Otherwise, both mean and sigma will be what you specify them to be. Or the 
#' user can specify their own prior function with which to input. If 
#' method="greta", the user-specified function must be a function from the 
#' greta package.
#' @param s2_prior Prior for the variance parameter. Requires list or function 
#' inputs. Defaults to a gamma(2,0.5). If Input is a list, the list must have 
#' two elements named 'shape' and 'scale' to represent the shape and scale of 
#' the gamma distribution respectively. Or the user can specify their own prior 
#' function with which to input. If method="greta", the user-specified function 
#' must be a function from the greta package. 
#' @param rho_prior The prior for the correlation of the adjacency matrix of 
#' the model. Requires list or function inputs. Default is a N(0.36,0.49) 
#' (Dittrich et al. 2017). If input is a list, the list must have two elements 
#' named 'mean' and 'sd' which are the mean and standard deviation of the normal 
#' distribution respectively. Or the user can specify their own prior 
#' function with which to input. If method="greta", the user-specified function 
#' must be a function from the greta package. If method="stan" only list input 
#' will be accepted.
#' @param ... Additional arguments to be passed into the optim() function for
#' when method="normal".
#' @return A list with two classes, the first being "nam" and the second is the 
#' class of the method that contains values:
#' @return \code{estimates} A list containing the estimates of the parameters 
#' of the chosen method.
#' @return \code{A} The adjacency matrix. Will be row-normalized if rownorm=TRUE.
#' @return \code{loglik} The log of the likelihood of the model.
#' @return \code{AIC} The Akaike Information Criterion of the model.
#' @return \code{BIC} The Bayesian Information Criterion of the model. 
#' 
#' @import greta
#' @import mvtnorm
#' @import rstan
#' @import truncnorm 
#' @export


# Model-specific wrapper function to call different methods 
ndm <- function(formula, data=NULL, network, method, 
                beta_prior = list(mean=0, sigma = 1), 
                s2_prior = list(shape=2, scale=0.5), 
                rho_prior = list(mean=0.36, sd=0.7^2), 
                ...){
  # Cleaning inputs
  X = model.matrix(formula, data=data)
  y = model.frame(formula, data)[,1]
  n = NROW(y)
  p = NCOL(X)

  # Check stan prior inputs are list only. 
  if(tolower(method)=="stan" && typeof(beta_prior) != "list" |
     tolower(method)=="stan" && typeof(s2_prior) != "list" |
     tolower(method)=="stan" && typeof(rho_prior) != "list"){
    stop("stan method requires list input priors if you want to change the 
         priors using MCMC use method=\"greta\"")
  }

  # Creating dataset input for individual methods
  ndmdat <- list(y=y, X=X, network=network, method=method, p=p, n=n, I=diag(n),
                 beta_prior=beta_prior, 
                 s2_prior = s2_prior, 
                 rho_prior = rho_prior, 
                 ...)

  if(tolower(method)=="normal"){
    ndmout <- ndm_normal(ndmdat)
  }else if(tolower(method)=="greta"){
    ndmout <- ndm_greta(ndmdat)
  }else if(tolower(method)=="stan"){
    ndmout <- ndm_stan(ndmdat)
  }else if(tolower(method)=="lnam"){
    ndmout <- ndm_lnam(ndmdat)
  }
  ndmout  
}

# Likelihood of the disturbances model. 
ndm_likelihood <- function(theta, A, y, X, log=TRUE){
  beta <- theta[1:(NROW(theta)-2)]
  s2 <- theta[NROW(theta)-1]
  rho <- theta[NROW(theta)]
  
  n <- NROW(y)
  p <- length(beta)
  I <- diag(n)
  
  Launtif <- qr.solve(I - rho*A)
  
  mvtnorm::dmvnorm(y, X %*% beta, s2 * tcrossprod(Launtif), log = log)
  
}

# Posterior of the effects model.
ndm_posterior <- function(theta, A, y, X, log=TRUE){
  beta <- theta[1:(NROW(theta)-2)]
  s2 <- theta[NROW(theta)-1]
  rho <- theta[NROW(theta)]
  
  n <- NROW(y)
  p <- length(beta)
  I <- diag(n)
  
  Launtif <- qr.solve(I - rho*A)
  
  out <- ndm_likelihood(theta, A, y, X, log=FALSE) *
    bprior(beta) * s2prior(s2) * rhoprior(rho)
  if(log){
    out <- log(out)
  }
  out
}

# Code for when method="normal"
ndm_normal <- function(ndmdat){
  
  # Checking Prior Inputs
  if(typeof(ndmdat$beta_prior) == "list"){
    if(length(ndmdat$beta_prior$mean) == 1){
      bmean <- rep(ndmdat$beta_prior$mean,ndmdat$p)
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
      dgamma(ndmdat$s2_prior$shape, ndmdat$s2_prior$scale)
    }
  }else{
    s2prior <- ndmdat$s2_prior
  }
  
  if(typeof(ndmdat$rho_prior) == "list"){
    rhoprior <- function(rho){
      truncnorm::dtruncnorm(rho, a=-1, b=1, mean=ndmdat$rho_prior$mean, 
                            sd=ndmdat$rho_prior$sd)
    }
  }else{
    rhoprior <- ndmdat$rho_prior
  }
  
  
  # Checking if additional Arguments were input
  if(!is.null(ndmdat$gr)){
    gr=ndmdat$gr
  }else{
    gr=NULL
  }
  
  if(!is.null(ndmdat$lower)){
    lower = ndmdat$lower
  }else{
    lower=-Inf
  }
  
  if(!is.null(ndmdat$upper)){
    upper = ndmdat$upper
  }else{
    upper=Inf
  }
  
  if(!is.null(ndmdat$control)){
    control=ndmdat$control
  }else{
    control=list()
  }
  
  # Finding good initial beta values
  init.betas <- lm(ndmdat$y ~ ndmdat$X[,-1])$coefficients
  
  # Finding posterior estimates using the normal approximation of the posterior
  myopt <- optim(par=c(init.betas, 1, 0.36), fn=ndm_posterior, 
                 method = "BFGS",
                 lower=lower,
                 upper=upper,
                 control=control,
                 hessian=TRUE,
                 y=nemdat$y,
                 X=nemdat$X,
                 A=nemdat$network)
  
  # Outputting the results
  list(beta = myopt$par[1:ndmdat$p], s2 = myopt$par[ndmdat$p + 1],
       rho = myopt$par[ndmdat$p + 2], network = ndmdat$network, 
       H = myopt$hessian, Hinv = qr.solve(myopt$hessian), logpost = myopt$value, 
       loglik = ndm_likelihood(myopt$par, ndmdat$network, ndmdat$y, ndmdat$X),
       y = ndmdat$y, X = ndmdat$X, method="normal", model="effects")
  
}

# Code for when method="greta"
ndm_greta <- function(ndmdat){
  # Checking/Creating Priors 
  if(typeof(ndmdat$beta_prior) == "list"){
    if(length(ndmdat$beta_prior$mean) == 1){
      bmean <- rep(ndmdat$beta_prior$mean,ndmdat$p)
    }else{
      bmean <- ndmdat$beta_prior$mean
    }
    
    if(length(ndmdat$beta_prior$sigma) == 1){
      bsigma <- diag(ndmdat$p) * ndmdat$beta_prior$sigma
    }else{
      bsigma <- ndmdat$beta_prior$sigma
    }
    
    beta = greta::multivariate_normal(matrix(bmean, nrow = 1), bsigma)

  }else{
    beta = ndmdat$beta_prior
  }
  
  if(typeof(ndmdat$s2_prior) == "list"){
    s2 = greta::gamma(2, 2)
  }else{
    s2 =  ndmdat$s2_prior
  }
  
  if(typeof(ndmdat$rho_prior) == "list"){
    rho = greta::normal(0.36, 0.7^2, truncation = c(-1,1))
  }else{
    rho = ndmdat$rho_prior
  }
  
  
  # Checking if additional Arguments were input
  if(!is.null(ndmdat$n_samples)){
    n_samples=ndmdat$n_samples
  }else{
    n_samples=5000
  }
  
  if(!is.null(ndmdat$chains)){
    chains = ndmdat$chains
  }else{
    chains=4
  }
  
  if(!is.null(ndmdat$one_by_one)){
    one_by_one=ndmdat$one_by_one
  }else{
    one_by_one=TRUE
  }
  
  if(!is.null(ndmdat$initial_values)){
    initial_values=ndmdat$initial_values
  }else{
    initial_values=greta::initials(s2=1, rho=0.36)
  }
  
  # set likelihood
  Launtif <- solve(ndmdat$I - rho * ndmdat$network)
  mu <- ndmdat$X %*% t(beta) 
  Sigma <- s2*tcrossprod(Launtif)
  y <- t(matrix(ndmdat$y, ncol = 1))
  greta::distribution(y) <- greta::multivariate_normal(
    mean = t(mu),
    Sigma = Sigma,
    dimension = NROW(ndmdat$n))
  
  # Run model and get draws
  m <- greta::model(beta, s2, rho)
  
  draws <- greta::mcmc(m, n_samples=n_samples, chains=chains, 
                       one_by_one = one_by_one, initial_values = initial_values)
  draws
  
}

# Code for when method="stan"
#ndm_stan <- function(ndmdat){
  
#}

# Code for when method="lnam"
ndm_lnam <- function(ndmdat){
  
  # Checking if additional arguments were input
  if(!is.null(ndmdat$theta.seed)){
    theta.seed=ndmdat$theta.seed
  }else{
    theta.seed=NULL
  }
  
  if(!is.null(ndmdat$null.model)){
    null.model = ndmdat$null.model
  }else{
    null.model="meanstd"
  }
  
  if(!is.null(ndmdat$control)){
    control=ndmdat$control
  }else{
    control=list()
  }
  
  if(!is.null(ndmdat$tol)){
    tol=ndmdat$tol
  }else{
    tol=1e-10
  }
  
  # Running lnam
  sna::lnam(y=ndmdat$y, x=ndmdat$X, W2=ndmdat$network, theta.seed=theta.seed,
            null.model=null.model, method="BFGS", control=control, tol=tol)
}
