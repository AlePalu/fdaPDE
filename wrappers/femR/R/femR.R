#' @importFrom Rcpp evalCpp
#' @useDynLib femR, .registration = TRUE
#' @name femR
NULL

#' Prova
#' @return A vector
#' @description Prova Prova Prova
#' @export
#' @examples
#' library(femR)
#' z <- hello_world()
hello_world <- function(points, edges, elements, neigh, boundary, diffusion= 1., transport= matrix(0,nrow=2,ncol=1),  reaction= 0.) {

  
  return(rcpp_hello_world(points, edges, elements, neigh, boundary, diffusion, transport, reaction))
}
