
#' convert.adjacency -  A function where you input some kind of adjacency object 
#' and you out the object class of your choice. 
#' 
#' Will fail if the object is too large to store in memory as class "matrix" 
#' and the original object is a network or out="matrix". (no method
#' for coercing a "network" object directly into a dgCMatrix)
#' 
#' Additionally, you cannot input a row-normalized matrix and output an igraph 
#' or edgelist, it will not be correct: try running row.norm(A, reverse=TRUE). 
#' @param A The adjacency object to convert of class type 
#' c("matrix", "dgCMatrix", "edgelist", "igraph", "network").
#' @param out The form of the adjacency object to output. Options are 
#' c("matrix", "Matrix", "edgelist", "igraph", "network"). 
#' @param df An optional data frame of node attributes. Only potentially 
#' important if out="igraph" so the outputted igraph object has node attributes.  
#' @return An object of class=out:
#' @return \code{A} The adjacency matrix.
#' @return \code{df} A data frame of node attributes. Mostly applicable if input
#' is an "igraph" object with node attributes. 
#' @import Matrix
#' @import igraph
#' @import intergraph
#' @export

convert.adjacency <- function(A, out="matrix", df=NULL){
  library(Matrix)
  # It's easiest to convert from a matrix object into everything else.
  # So I first convert to a matrix object then convert to the class of choice.
  convert.to.igraph <- function(A2){
  
    if(class(A2)[1]=="matrix" & NCOL(A2) == NROW(A2)){
      if("names" %in% colnames(df) | "name" %in% colnames(df)){
        vert <- df
      }else{
        vert <- cbind(names=1:NROW(A2),df)
      }
      
      if(all(abs(rowSums(A2) - 1) < 2*.Machine$double.eps | rowSums(A2) == 0)){
        A2 <- row.norm(A2, TRUE)
      }
      
      step1 <- igraph::as_edgelist(igraph::graph_from_adjacency_matrix(A2))
      outA2 <- igraph::graph_from_data_frame(step1, directed = TRUE, 
                                             vertices = vert)
      
    }else if(class(A2)[1]=="network"){
      outA2 <- intergraph::asIgraph(A2)
      
    }else if(class(A2)[1]=="dgCMatrix" | class(A2)[1]=="dgeMatrix"){
      if("names" %in% colnames(df) | "name" %in% colnames(df)){
        vert <- df
      }else{
        vert <- cbind(names=1:NROW(A2),df)
      }
      if(all(abs(rowSums(A2) - 1) < 2*.Machine$double.eps | rowSums(A2) == 0)){
        A2 <- row.norm(A2, TRUE)
      }
      step1 <- igraph::as_edgelist(
        igraph::graph_from_adjacency_matrix(as.matrix(A2)))
      outA2 <- igraph::graph_from_data_frame(step1, directed = TRUE, 
                                             vertices = vert)
      
    }else if(class(A2)[1] == "igraph"){
      outA2 <- A2
      
    }else if(class(A2)[1] == "matrix" | class(A2)[1] == "data.frame"){
      if("names" %in% colnames(df) | "name" %in% colnames(df)){
        vert <- df
      }else if(!(is.null(df))){
        vert <- cbind(names=1:NROW(df),df)
      }else{
        vert <- df
      }
      outA2 <- igraph::graph_from_data_frame(A2, directed = TRUE, 
                                             vertices = vert)
    }
    
    outA2
  }
  
  
  if(out=="matrix"){
    step1 <- convert.to.igraph(A)
    outdf <- igraph::vertex_attr(step1)
    step2 <- igraph::as_adjacency_matrix(step1)
    outA <- list(A=as.matrix(step2), df=outdf)
  }else if(out == "Matrix"){
    step1 <- convert.to.igraph(A)
    outdf <- igraph::vertex_attr(step1)
    outA <- list(A=igraph::as_adjacency_matrix(step1), df=outdf)
  }else if(out == "edgelist"){
    step1 <- convert.to.igraph(A)
    outdf <- igraph::vertex_attr(step1)
    outA <- list(A=igraph::as_edgelist(step1), df=outdf)
  }else if(out == "igraph"){
    step1 <- convert.to.igraph(A)
    outdf <- igraph::vertex_attr(step1)
    outA <- list(A=step1, df=outdf)
  }else if(out == "network"){
    step1 <- convert.to.igraph(A)
    outdf <- igraph::vertex_attr(step1)
    outA <- list(A=intergraph::asNetwork(step1), df=outdf)
  }
  outA
}
