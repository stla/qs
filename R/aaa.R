#' @useDynLib qs, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @import onion
#' @noRd
NULL

#' Title
#'
#' @param n x
#'
#' @return x
#' @export
#'
#' @examples
#' 3
rversor <- function(n){
  as.quaternion(rversor_cpp(n))
}