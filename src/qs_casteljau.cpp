#include "qs.h"

std::vector<qtrn> _select_segment_and_normalize_t(
    std::vector<std::vector<qtrn>> segments,
    Rcpp::NumericVector keyTimes,
    double t,
    double* time,
    double* difftime) {
  const std::size_t idx = _check_time(t, keyTimes, false);
  const double t0 = keyTimes[idx];
  const double t1 = keyTimes[idx+1];
  const double delta_t = t1 - t0;
  *difftime = delta_t;
  *time = (t - t0) / delta_t;
  return segments[idx];
}

qtrn getQuaternion(Rcpp::NumericVector qR){
  qtrn q(qR(0), qR(1), qR(2), qR(3));
  return q;
}

std::vector<qtrn> _reduce_de_casteljau(std::vector<qtrn> segment, double t){
  size_t l = segment.size();
  if(l < 2){
    Rcpp::stop("Segment must have at least two quaternions.");
  }
  while(l > 2){
    std::vector<qtrn> newsegment(l - 1);
    for(std::size_t i = 0; i < l-1; i++){
      qtrn one = segment[i];
      qtrn two = segment[i+1];
      newsegment[i] = one.slerp(t, two);
    }
    segment = newsegment;
    l--;
  }
  return segment;
}

//   while(length(segment) > 2L){
//     newsegment <- quaternion(length.out = length(segment) - 1L)
//     for(i in seq_len(length(segment)-1L)){
//       one <- segment[i]
//       two <- segment[i+1L]
//       newsegment[i] <- .slerp(one, two, t)
//     }
//     segment <- newsegment
//   }
//   segment
// }

// {}
// // [[Rcpp::export]]