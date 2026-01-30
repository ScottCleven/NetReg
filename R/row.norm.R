#' row.norm -  A function where you input a matrix or dgCMatrix and you receive
#' the row-normalized version of it (or the reverse). 
#' @param Netowrk The adjacency object to row normalize (or un-row normalize).
#' @param reverse If TRUE, converts a row-normalized matrix into one that isn't
#' row normalized. 
#' @param sparse Only applies if reverse=TRUE. Conditional if you want the 
#' output to be a sparse matrix. If reverse=FALSE, the output will maintain the 
#' class of the input. 
#' @return A matrix or dgCMatrix that is either row normalized or not depending 
#' on reverse:
#' @export


row.norm <- function(Network, reverse=FALSE, sparse=FALSE){
  if(reverse){
    if(sparse){
      x <- rowSums(Network > 0)
      sparseMatrix(i=seq(x), j=seq(x), x=x)
      
    }else{
      diag(rowSums(Network > 0)) %*% Network
    }
    
  }else{
    Network/ifelse(rowSums(Network) == 0,1,rowSums(Network))
  }
}
