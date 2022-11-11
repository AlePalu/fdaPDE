
#' @useDynLib femR, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @name femR

#' Prova
#' @return A vector
#' @description Prova Prova Prova
#' @export
#' @examples
#' library(femR)
#' z <- hello_world()
hello_world <- function() {
  return(rcpp_hello_world())
}
