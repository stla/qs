getQMatrix <- function(quaternions){
  stopifnot(is.quaternion(quaternions))
  as.matrix(quaternions)
}

isNumericVector <- function(x){
  is.numeric(x) && !anyNA(x)
}

#' Title
#'
#' @param segments x
#' @param keyTimes x
#' @param times x
#'
#' @return x
#' @export
#'
#' @examples
#' 2+2
DeCasteljau <- function(
  segments, keyTimes = NULL, times
){
  stopifnot(is.list(segments))
  stopifnot(is.null(keyTimes) || isNumericVector(keyTimes))
  stopifnot(isNumericVector(times))
  segments <- lapply(segments, getQMatrix)
  Q <- DeCasteljau_cpp(segments, keyTimes, times)
  as.quaternion(Q)
}
  