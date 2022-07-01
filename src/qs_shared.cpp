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
    seq(i-1) = (double)(i);
  }
  return seq;
}

Rcpp::NumericVector _check_keyTimes(Rcpp::NumericVector keyTimes,
                                    std::size_t n_quaternions) {
  if(keyTimes == R_NilValue) {
    return _seq_len(n_quaternions);
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
    idx = _findInterval(t, keyTimes) - 1;
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

// qtrn slerp(qtrn q0, qtrn q1, double t) {
//   qtrn q2 = q1 * q0.inverse();
//   qtrn H1(1.0, 0.0, 0.0, 0.0);
//   qtrn q2powt = q2.slerp(t, H1);
//   return q2powt * q0;
// }

qtrn qpower(qtrn q, double t) {
  //qtrn H1(1.0, 0.0, 0.0, 0.0);
  double w = q.w();
  if(w == 1){
    return q; // ok because versor assumption => q = H1
  }
  double alpha = t*acos(w);
  double a = sin(alpha)/sqrt(1-w*w);
  double b = cos(alpha);
  qtrn qpowt(b, a * q.x(), a * q.y(), a * q.z());
  return qpowt;
}

qtrn slerp(qtrn q0, qtrn q1, double t) {
  qtrn q2 = q1 * q0.inverse();
  qtrn q2powt = qpower(q2, t);
  return q2powt * q0;
}

// // [[Rcpp::export]]  c(v*sin(t*acos(w))/sqrt(1-w*w), cos(t*acos(w))) 
// Rcpp::NumericMatrix	slerp_(const Rcpp::NumericVector & q1, 
//                            const Rcpp::NumericVector & q2, 
//                            const Rcpp::NumericVector & t)
// {
//   if(q1.size() != 4 || q2.size() != 4){
//     throw Rcpp::exception("q1 and q2 must be quaternions");
//   }
//   Eigen::Quaterniond qa(q1[0], q1[1], q1[2], q1[3]);
//   Eigen::Quaterniond qb(q2[0], q2[1], q2[2], q2[3]);
//   Rcpp::NumericMatrix out(t.size(), 4);
//   for(unsigned i=0; i<t.size(); i++){
//     Eigen::Quaterniond q = qa.slerp(t[i], qb);
//     out(i,0) = q.w();
//     out(i,1) = q.x();
//     out(i,2) = q.y();
//     out(i,3) = q.z();
//   }
//   return out;
// }

// [[Rcpp::export]]
Rcpp::NumericMatrix rversor_cpp(std::size_t n){
  Rcpp::NumericMatrix out(4, n);
  for(std::size_t i = 0; i < n; i++){
    qtrn q = qtrn::UnitRandom();
    out(0, i) = q.w();
    out(1, i) = q.x();
    out(2, i) = q.y();
    out(3, i) = q.z();
  }
  return out;
}  

//   newtimes
// }

// {}