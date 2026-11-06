#' row_norm -  A function where you input a matrix or dgCMatrix and you receive
#' the row-normalized version of it (or the reverse).
#' @param network The adjacency object to row normalize (or un-row normalize).
#' @param reverse If TRUE, converts a row-normalized matrix into one that isn't
#' row normalized.
#' @param sparse Conditional if you want the
#' output to be a sparse matrix. If sparse=FALSE, the output will maintain the
#' class of the input.
#' @return A matrix or dgCMatrix that is either row normalized or not depending
#' on reverse:
#' @import Matrix
#' @export


row_norm <- function(network, reverse=FALSE, sparse=FALSE){
  if(reverse){
    x <- Matrix::rowSums(network > 0)
    if(inherits(network, "dgCMatrix") || sparse){
      return(Matrix::Matrix(Matrix::Diagonal(x = x) %*% network, sparse=TRUE))

    }else{
      return(diag(x) %*% as.matrix(network))
    }

  }else{
    row_sums <- Matrix::rowSums(network)
    row_sums[row_sums == 0] <- 1
    if(inherits(network, "dgCMatrix")){
      return(Matrix::Diagonal(x = 1 / row_sums) %*% network)
    }else{
      if(sparse){
        return(Matrix::Matrix(network / row_sums, sparse = TRUE))
      }else{
        return(network / row_sums)
      }
    }
  }
}
