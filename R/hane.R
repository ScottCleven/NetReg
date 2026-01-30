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
#' "greta", "stan", "lnam"). Run ?methods.nam() for model-specific designations.
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
#' method="greta", the user-specified function must be a function from the
#' greta package.
#' @param gamma_prior The prior for the gamma vector of the model. Requires list
#' or function inputs. Default is a multivariate normal with mean=0 and
#' sigma=2.25^2. If you supply a list, you must have one element named
#' mean' which is the mean of the multivariate normal and one named 'sigma'
#' which is the variance-covariance matrix. If mean is a single integer,
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
#' @import greta
#' @import mvtnorm
#' @import rstan
#' @import truncnorm
#' @export


# Model-specific wrapper function to call different methods
hane <- function(formula, data, network, method, Usample,
                beta_prior = list(mean=0, sigma = 1),
                gamma_prior=list(mean=0, sigma=2.25^2),
                s2_prior = list(shape=2, scale=0.5),
                rho_prior = list(mean=0.36, sd=0.7^2),
                Usample.maxiter=100, Usample.eps=1e-8, init_vec=NULL,
                ...){

  X = model.matrix(formula, data=data)
  y = model.frame(formula, data)[,1]
  n = NROW(y)
  p = NCOL(X)
  nUsamp = dim(Usample)[1]
  D = dim(Usample)[3]

  # Check stan prior inputs are list only.
  if(tolower(method)=="stan" && typeof(beta_prior) != "list" |
     tolower(method)=="stan" && typeof(s2_prior) != "list" |
     tolower(method)=="stan" && typeof(rho_prior) != "list"){
    stop("stan method requires list input priors if you want to change the
         priors using MCMC use method=\"greta\"")
  }

  # Check if U is an array
  if(length(dim(Usample)) != 3 | dim(Usample)[2] != n){
    stop("Error: Input of Usample not in the right array format (size, n, D)")
  }

  # Calculate matrix normal approximation to prior of U
  Uapprox <- Uprior_appox(Usample, Usample.maxiter, Usample.eps)
  if(!Uapprox$conv) stop('Approximation to Usample failed! Consider increasing
                         Usample.maxiter or decreasing Usample.eps')
  Lambda  <- Uapprox$Lambda
  Omega   <- Uapprox$Omega
  Psi     <- Uapprox$Psi


  hanedat <- list(y=y, X=X, network=network, method=method, p=p, n=n, I=diag(n),
                 Lambda=Lambda, Omega=Omega, Psi=Psi, D=D,
                 beta_prior=beta_prior,
                 gamma_prior=gamma_prior,
                 s2_prior = s2_prior,
                 rho_prior = rho_prior,
                 ...)

  if(tolower(method)=="normal"){
    haneout <- hane_normal(hanedat)
  }else if(tolower(method)=="greta"){
    haneout <- hane_greta(hanedat)
  }else if(tolower(method)=="stan"){
    haneout <- hane_stan(hanedat)
  }
  haneout
}

# Likelihood of the HANE model.
hane_likelihood <- function(theta, A, y, X, Lambda, Omega, log=TRUE){

  n <- NROW(y)
  p <- length(beta)
  I <- diag(n)
  D <- ncol(Lambda)

  beta <- theta[1:p]
  gamma <- theta[(p+1):(p+D)]
  s2 <- theta[p+D+1]
  rho <- theta[p+D+2]



  Launtif <- qr.solve(I - rho*A)

  mvtnorn::dmvnorm(y,
                   (Launtif %*% (X %*% beta + Lambda%*%gamma))[,1],
                   as.matrix(Launtif%*%tcrossprod(crossprod(gamma)[1]*
                                                    Omega + diag(s2, n),
                                                  Launtif)),
                   log = log)

}

# Posterior of the HANE model.
hane_posterior <- function(theta, A, y, X, Lambda, Omega, log=TRUE){

  n <- NROW(y)
  p <- length(beta)
  I <- diag(n)
  D <- ncol(Lambda)

  beta    <- theta[1:p]
  gamma   <- theta[(p+1):(p+D)]
  sigmasq <- theta[p+D+1]
  rho     <- theta[p+D+2]


  out <- hane_likelihood(theta, A, y, X, Lambda, Omega, log=FALSE) *
    bprior(beta) * gprior(gamma) * s2prior(s2) * rhoprior(rho)
  if(log){
    out <- log(out)
  }
  out
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

  # Checking Prior Inputs
  if(typeof(hanedat$beta_prior) == "list"){
    if(length(hanedat$beta_prior$mean) == 1){
      bmean <- rep(hanedat$beta_prior$mean,p)
    }else{
      bmean <- hanedat$beta_prior$mean
    }

    if(length(hanedat$beta_prior$sigma) == 1){
      bsigma <- diag(p) * hanedat$beta_prior$sigma
    }else{
      bsigma <- hanedat$beta_prior$sigma
    }

    bprior <- function(beta){
      mvtnorm::dmvnorm(beta, bmean, bsigma)
    }
  }else{
    bprior <- hanedat$beta_prior
  }

  if(typeof(hanedat$gamma_prior) == "list"){
    if(length(hanedat$gamma_prior$mean) == 1){
      gmean <- rep(hanedat$gamma_prior$mean,p)
    }else{
      gmean <- hanedat$gamma_prior$mean
    }

    if(length(hanedat$gamma_prior$sigma) == 1){
      gsigma <- diag(p) * hanedat$gamma_prior$sigma
    }else{
      gsigma <- hanedat$gamma_prior$sigma
    }

    gprior <- function(gamma){
      mvtnorm::dmvnorm(gamma, gmean, gsigma)
    }
  }else{
    gprior <- hanedat$gamma_prior
  }

  if(typeof(hanedat$s2_prior) == "list"){
    s2prior <- function(s2){
      dgamma(hanedat$s2_prior$shape, hanedat$s2_prior$scale)
    }
  }else{
    s2prior <- hanedat$s2_prior
  }

  if(typeof(hanedat$rho_prior) == "list"){
    rhoprior <- function(rho){
      truncnorm::dtruncnorm(rho, a=-1, b=1, mean=hanedat$rho_prior$mean,
                            sd=hanedat$rho_prior$sd)
    }
  }else{
    rhoprior <- hanedat$rho_prior
  }


  # Checking if additional Arguments were input
  if(!is.null(hanedat$gr)){
    gr=hanedat$gr
  }else{
    gr=NULL
  }

  if(!is.null(hanedat$lower)){
    lower = hanedat$lower
  }else{
    lower=c(rep(-Inf, hanedat$p+hanedat$D), 1e-5,  -0.999)
  }

  if(!is.null(hanedat$upper)){
    upper = hanedat$upper
  }else{
    upper=c(rep(Inf, hanedat$p+hanedat$D+1), 0.999)
  }

  if(!is.null(hanedat$control)){
    control=hanedat$control
  }else{
    control=c(list(fnscale=-1), optim.control)
  }

  # Finding good initial values
  if(is.null(init_vec)){
    fit2sls  <- lnam2sls_effect(hanedat$y, cbind(hanedat$X, hanedat$Lambda),
                                hanedat$network)
    init_vec <- c(fit2sls$coefs[1:(p+D)], fit2sls$s2, fit2sls$coefs[p+D+1])
  }

  # Finding posterior estimates using the normal approximation of the posterior
  myopt <- optim(init_vec,
                 fn=hane_posterior,
                 method = "L-BFGS-B",
                 control=control,
                 hessian = T,
                 lower=lower, upper=upper,
                 y=hanedat$y, X=hanedat$X, A=hanedat$A,
                 Lambda=hanedat$Lambda, Omega=hanedat$Omega, Psi=hanedat$Psi,
                 ...)

  # Outputting the results
  list(beta = myopt$par[1:hanedat$p], s2 = myopt$par[hanedat$p + 1],
       rho = myopt$par[hanedat$p + 2], network = hanedat$network,
       H = myopt$hessian, Hinv = qr.solve(myopt$hessian), logpost = myopt$value,
       loglik = hane_likelihood(myopt$par, hanedat$network, hanedat$y, hanedat$X),
       y = hanedat$y, X = hanedat$X, method="normal", model="effects")

}

# Code for when method="greta"
hane_greta <- function(hanedat){
  # Checking/Creating Priors
  if(typeof(hanedat$beta_prior) == "list"){
    if(length(hanedat$beta_prior$mean) == 1){
      bmean <- rep(hanedat$beta_prior$mean,hanedat$p)
    }else{
      bmean <- hanedat$beta_prior$mean
    }

    if(length(hanedat$beta_prior$sigma) == 1){
      bsigma <- diag(hanedat$p) * hanedat$beta_prior$sigma
    }else{
      bsigma <- hanedat$beta_prior$sigma
    }

    bprior <- function(beta){
      beta = greta::multivariate_normal(matrix(bmean, nrow = 1), bsigma)
    }
  }else{
    beta = hanedat$beta_prior
  }

  if(typeof(hanedat$gamma_prior) == "list"){
    if(length(hanedat$gamma_prior$mean) == 1){
      gmean <- rep(hanedat$gamma_prior$mean,hanedat$D)
    }else{
      gmean <- hanedat$gamma_prior$mean
    }

    if(length(hanedat$gamma_prior$sigma) == 1){
      gsigma <- diag(hanedat$D) * hanedat$gamma_prior$sigma
    }else{
      gsigma <- hanedat$gamma_prior$sigma
    }

    gprior <- function(gamma){
      gamma = greta::multivariate_normal(matrix(gmean, nrow = 1), gsigma)
    }
  }else{
    gamma = hanedat$gamma_prior
  }


  if(typeof(hanedat$s2_prior) == "list"){
    s2 = greta::gamma(2, 2)
  }else{
    s2 =  hanedat$s2_prior
  }

  if(typeof(hanedat$rho_prior) == "list"){
    rho = greta::normal(0.36, 0.7^2, truncation = c(-1,1))
  }else{
    rho = hanedat$rho_prior
  }


  # Checking if additional Arguments were input
  if(!is.null(hanedat$n_samples)){
    n_samples=hanedat$n_samples
  }else{
    n_samples=5000
  }

  if(!is.null(hanedat$chains)){
    chains = hanedat$chains
  }else{
    chains=4
  }

  if(!is.null(hanedat$one_by_one)){
    one_by_one=hanedat$one_by_one
  }else{
    one_by_one=TRUE
  }

  if(!is.null(hanedat$initial_values)){
    initial_values=hanedat$initial_values
  }else{
    initial_values=greta::initials(s2=1, rho=0.36)
  }

  # set likelihood
  Launtif <- solve(hanedat$I - rho * hanedat$network)
  mu <- Launtif %*% hanedat$X %*% beta
  Sigma <- s2*Launtif %*% t(Launtif)
  y <- t(matrix(hanedat$y, ncol = 1))
  greta::distribution(y) <- greta::multivariate_normal(
    mean = t(mu),
    Sigma = Sigma,
    dimension = NROW(hanedat$n))

  # Run model and get draws
  m <- greta::model(beta, s2, rho)

  draws <- greta::mcmc(m, n_samples=n_samples, chains=chains,
                       one_by_one = one_by_one, initial_values = initial_values)
  draws

}

# Code for when method="stan"
hane_stan <- function(hanedat){

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
