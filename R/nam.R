#' nam - A wrapper function for performing the different network autocorrelation
#'  family of models
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
#' @param model A character indicating the model to be run. options are 
#' c("effects", "disturbances", "NMR", "IV", "HA"). Run ?models.nam()
#' for a more detailed explanation and citations on each method. 
#' @param method A character indicating the method to use for simulation. 
#' Options vary by model but all potential options are c("normal", 
#' "greta", "stan", "lnam"). Run ?methods.nam() for model-specific designations. 
#' @param rownorm A logical indicating whether the inputted adjacency matrix 
#' should be row-normalized (TRUE). Ignored by certain models that require 
#' normalization.
#' @param ... Additional arguments that are method-specific. Run 
#' ?args.nam() for a more detailed explanation on the different arguments for
#' each method.   
#' @return A list with two classes, the first being "nam" and the second is the 
#' class of the method that contains values:
#' @return \code{estimates} A list containing the estimates of the parameters 
#' of the chosen method.
#' @return \code{A} The adjacency matrix. Will be row-normalized if rownorm=TRUE.
#' @return \code{loglik} The log of the likelihood of the model.
#' @return \code{AIC} The Akaike Information Criterion of the model.
#' @return \code{BIC} The Bayesian Information Criterion of the model.  
#' @export

nam <- function(formula, data=NULL, network, model="effects", method="normal", 
                rownorm=TRUE, ...){
  
  # Formatting the data.
  if(class(data)=="igraph"){
    df <- convert.adjacency(data)
    data <- df$df
    network <- df$A
  }else{
    network <- convert.adjacency(network)$A
  }
  
  
  # Wrapper to distinguish the model-specific cleaning and function call.
  if(tolower(model)=="effects" | tolower(model)=="nem"){
    if(rownorm){
      network <- row.norm(network)
    }
    out <- nem(formula=formula, data=data, network=network, model=model, 
                       method=method, rownorm=rownorm, ...)
  }else if(tolower(model) == "disturbances" | tolower(model)=="ndm"){
    if(rownorm){
      network <- row.norm(network)
    }
    out <- ndm(formula=formula, data=data, network=network, model=model, 
                       method=method, rownorm=rownorm, ...)
  }else if(tolower(model) == "iv"){
    if(rownorm){
      network <- row.norm(network)
    }
    out <- iv(formula=formula, data=data, network=network, model=model, 
                       method=method, rownorm=rownorm, ...)
  }else if(tolower(model) == "nmr"){
    network <- row.norm(network)
    out <- nmr(formula=formula, data=data, A=network, model=model, 
                       method=method, rownorm=rownorm, ...)
  }else if(tolower(model) == "ha" | tolower(model) == "homopholy-adjusted"){
    network <- row.norm(network)
    out <- ha(formula=formula, data=data, network=network, model=model, 
               method=method, rownorm=rownorm, ...)
  }else{
    stop("You did not give a valid model name.")
  }

  # Class assignment and output.
  class(out) <- "nam"
  return(out)
}





