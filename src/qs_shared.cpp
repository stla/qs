#include "qs.h"

std::vector<qtrn> _canonicalized(std::vector<qtrn> quaternions) {
  const std::size_t n = quaternions.size();
  std::vector<qtrn> out(n);
  qtrn p(1.0, 0.0, 0.0, 0.0);
  for(std::size_t i = 0; i < n; i++) {
    qtrn q = quaternions[i];
    if(p.dot(q) < 0.0) {
      q.w() = -q.w();
      q.vec() = -q.vec();
    }
    out[i] = q;
    p = q;
  }
  return out;
}

std::vector<qtrn> _check_keyRotors(std::vector<qtrn> keyRotors, bool closed) {
  if(keyRotors.size() < 2) {
    Rcpp::stop("At least two key rotors are required.");
  }
  if(closed) {
    keyRotors.push_back(keyRotors[0]);
  }
  return keyRotors;
}

Rcpp::NumericVector _seq_len(std::size_t n) {
  Rcpp::NumericVector seq(n);
  for(std::size_t i = 1; i <= n; i++) {
    seq(i) = (double)(i);
  }
  return seq;
}

Rcpp::NumericVector _check_keyTimes(Rcpp::NumericVector keyTimes,
                                    std::size_t n_quaternions) {
  if(keyTimes == R_NilValue) {
    return (_seq_len(n_quaternions));
  }
  std::size_t n = keyTimes.size();
  for(std::size_t i = 1; i < n; i++) {
    if(keyTimes(i) - keyTimes(i - 1) <= 0) {
      Rcpp::stop("`keyTimes` must be an increasing vector of numbers.");
    }
  }
  return keyTimes;
}

std::size_t _findInterval(double x, Rcpp::NumericVector vec) {
  size_t n = vec.size();
  if(x > vec(n - 1)) {
    return n;
  }
  std::size_t idx = 0;
  for(std::size_t i = 0; i < n - 1; i++) {
    if(x >= vec(i)) {
      idx++;
    } else {
      break;
    }
  }
  return idx;
}

std::size_t _check_time(double t, Rcpp::NumericVector keyTimes, bool special) {
  std::size_t n_keyTimes = keyTimes.size();
  double lastKeyTime = keyTimes(n_keyTimes - 1);
  if(t < keyTimes(0) || t > lastKeyTime) {
    Rcpp::stop(
        "The interpolating times must be within the range of `keyTimes`.");
  }
  std::size_t idx;
  if(t < lastKeyTime) {
    idx = _findInterval(t, keyTimes);
  } else {  // t = lastKeyTime
    if(special) {
      idx = n_keyTimes - 3;  // 2?
    } else {
      idx = n_keyTimes - 2;  // 1?
    }
  }
  return idx;
}

Rcpp::NumericVector _head(Rcpp::NumericVector v) {
  const std::size_t n = v.size();
  return Rcpp::head(v, n - 1);
}

Rcpp::NumericVector _seq(double a, double b, std::size_t l) {
  Rcpp::NumericVector out(l);
  const double delta = (b - a) / ((double)(l - 1));
  double current = a;
  for(size_t i = 0; i < l; i++) {
    out(i) = current;
    current += delta;
  }
  return out;
}

Rcpp::NumericVector _interpolateTimes(Rcpp::NumericVector times,
                                      std::size_t n,
                                      bool last) {
  const std::size_t n_times = times.size();
  const std::size_t len_out =
      last ? (n * (n_times - 1) + 1) : (n * (n_times - 1));
  Rcpp::NumericVector newtimes(len_out);
  std::size_t k = 0;
  for(std::size_t i = 0; i < n_times - 1; i++) {
    const Rcpp::NumericVector vi = _head(_seq(times(i), times(i + 1), n + 1));
    for(std::size_t j = 0; j < n; i++) {
      newtimes(k) = vi(j);
      k++;
    }
  }
  if(last) {
    newtimes(k) = times(n_times - 1);
  }
  return newtimes;
}

//   newtimes
// }

// {}