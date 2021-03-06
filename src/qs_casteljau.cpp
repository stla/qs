#include "qs.h"

std::vector<qtrn> _select_segment_and_normalize_t(
    std::vector<std::vector<qtrn>> segments,
    Rcpp::NumericVector keyTimes,
    double t,
    double* time,
    double* difftime) {
  const std::size_t idx = _check_time(t, keyTimes, false);
  const double t0 = keyTimes[idx];
  const double t1 = keyTimes[idx + 1];
  const double delta_t = t1 - t0;
  *difftime = delta_t;
  *time = (t - t0) / delta_t;
  return segments[idx];
}

qtrn _getRQuaternion(Rcpp::NumericVector qR) {
  qtrn quat(qR(0), qR(1), qR(2), qR(3));
  return quat;
}

std::vector<qtrn> _getRQuaternions(Rcpp::NumericMatrix Q) {
  std::size_t n = Q.ncol();
  std::vector<qtrn> quats(n);
  for(std::size_t j = 0; j < n; j++) {
    quats[j] = _getRQuaternion(Q(Rcpp::_, j));
  }
  return quats;
}

std::vector<std::vector<qtrn>> _getRSegments(Rcpp::List rsegments) {
  std::size_t nsegments = rsegments.size();
  std::vector<std::vector<qtrn>> segments(nsegments);
  for(std::size_t i = 0; i < nsegments; i++) {
    Rcpp::NumericMatrix segment = Rcpp::as<Rcpp::NumericMatrix>(rsegments(i));
    segments[i] = _getRQuaternions(segment);
  }
  return segments;
}

Rcpp::NumericVector _getCQuaternion(qtrn quat) { 
  return Rcpp::NumericVector::create(quat.w(), quat.x(), quat.y(), quat.z());
} 

Rcpp::NumericMatrix _getCQuaternions(std::vector<qtrn> quats) { 
  std::size_t n = quats.size();
  Rcpp::NumericMatrix Q(4, n);
  for(std::size_t j = 0; j < n; j++) {
    Rcpp::NumericVector qR = _getCQuaternion(quats[j]);
    Q(Rcpp::_, j) = qR;
  }
  return Q;
} 

std::vector<qtrn> _reduce_de_casteljau(std::vector<qtrn> segment, double t) {
  size_t l = segment.size();
  if(l < 2) {
    Rcpp::stop("Segment must have at least two quaternions.");
  }
  while(l > 2) {
    std::vector<qtrn> newsegment(l - 1);
    for(std::size_t i = 0; i < l - 1; i++) {
      qtrn one = segment[i];
      qtrn two = segment[i + 1];
      newsegment[i] = slerp(one, two, t);
    }
    segment = newsegment;
    l--;
  }
  return segment;
}

qtrn _eval_casteljau_single(double t,
                            std::vector<std::vector<qtrn>> segments,
                            Rcpp::NumericVector keyTimes) {
  double time, difftime;  // difftime not used here
  std::vector<qtrn> segment =
      _select_segment_and_normalize_t(segments, keyTimes, t, &time, &difftime);
  std::vector<qtrn> quats = _reduce_de_casteljau(segment, time);
  return slerp(quats[0], quats[1], time);
}

std::vector<qtrn> _eval_casteljau_vector(
    Rcpp::NumericVector times,
    std::vector<std::vector<qtrn>> segments,
    Rcpp::NumericVector keyTimes) {
  std::size_t n = times.size();
  std::vector<qtrn> quats(n);
  for(std::size_t i = 0; i < n; i++) {
    quats[i] = _eval_casteljau_single(times(i), segments, keyTimes);
  }
  return quats;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix DeCasteljau_cpp(
  Rcpp::List rsegments, Rcpp::NumericVector keyTimes, Rcpp::NumericVector times
){
  std::size_t nsegments = rsegments.size();
  std::size_t nkeyTimes = keyTimes.size();
  if(nkeyTimes == 0){
    keyTimes = _seq_len(nsegments + 1);
  }else if(nkeyTimes != nsegments + 1){
    Rcpp::stop("Number of key times must be one more than number of segments.");
  }
  std::vector<std::vector<qtrn>> segments = _getRSegments(rsegments);
  std::vector<qtrn> quats = _eval_casteljau_vector(times, segments, keyTimes);
  return _getCQuaternions(quats);
}

std::array<qtrn, 2> _calculate_control_quaternions(
  std::array<qtrn, 3> quaternions, std::array<double, 3> times, 
  double t, double c, double b
){
  qtrn q_1 = quaternions[0];
  qtrn q0  = quaternions[1];
  qtrn q1  = quaternions[2];
  double t_1 = times[0];
  double t0  = times[1];
  double t1  = times[2];
  if(t_1 == t0 || t0 == t1){
    return {q0, q0};
  }
  double A = (1 - t) * (1 + c) * (1 + b);
  double B = (1 - t) * (1 - c) * (1 - b);
  double C = (1 - t) * (1 - c) * (1 + b);
  double D = (1 - t) * (1 + c) * (1 - b);
  qtrn lq_in = qlog(q0 * q_1.inverse());
  qtrn lq_out = qlog(q1 * q0.inverse()); 
  Rcpp::NumericVector v_in = 
    Rcpp::NumericVector::create(0.0, lq_in.x(), lq_in.y(), lq_in.z());
  Rcpp::NumericVector v_out = 
    Rcpp::NumericVector::create(0.0, lq_out.x(), lq_out.y(), lq_out.z());
  v_in  = v_in / (t0 - t_1);
  v_out = v_out / (t1 - t0);
  Rcpp::NumericVector v0CD = 
    (C * (t1 - t0) * v_in + D * (t0 - t_1) * v_out) / (t1 - t_1);
  Rcpp::NumericVector v0AB = 
    (A * (t1 - t0) * v_in + B * (t0 - t_1) * v_out) / (t1 - t_1);
  std::array<qtrn, 2> out = {
    qexp(_getRQuaternion((t_1 - t0) * v0CD / 3.0)) * q0,
    qexp(_getRQuaternion((t1 - t0) * v0AB / 3.0)) * q0
  };
  return out;
}


template<typename T>
std::vector<std::array<T, 3>> makeTriplets(std::vector<T> vec) {
  const std::size_t n = vec.size();
  std::vector<std::array<T, 3>> triplets(n-2);
  for(size_t i = 0; i < n - 2; i ++) {
    triplets[i] = {vec[i], vec[i+1], vec[i+2]};
  }
  return triplets;
}

std::vector<std::array<double, 3>> makeTriplets_times(std::vector<double> times, bool closed) {
  if(closed){
    const std::size_t ntimes = times.size();
    times.insert(times.begin(), times[0] - (times[ntimes-1] - times[ntimes-2]));
    times.push_back(times[ntimes-1] + (times[1] - times[0]));
  }
  return makeTriplets<double>(times);
}

std::vector<std::array<qtrn, 3>> makeTriplets_rotors(std::vector<qtrn> rotors, bool closed) {
  if(closed){
    const std::size_t nrotors = rotors.size();
    qtrn prefix = rotors[nrotors - 2];
    if(prefix.dot(rotors[0]) < 0.0){
      prefix.w() = -prefix.w();
      prefix.vec() = -prefix.vec();
    }
    qtrn suffix = rotors[1];
    if(suffix.dot(rotors[nrotors-1]) < 0.0){
      suffix.w() = -suffix.w();
      suffix.vec() = -suffix.vec();
    }
    rotors.insert(rotors.begin(), prefix);
    rotors.push_back(suffix);
  }
  return makeTriplets<qtrn>(rotors);
}

qtrn _natural_control_quaternion(qtrn outer, qtrn inner_control, qtrn inner){
  return qpower((inner_control * inner.inverse()) * (inner * outer.inverse()), 0.5) * outer;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix KochanekBartels_cpp(
  Rcpp::NumericMatrix keyRotorsR, Rcpp::NumericVector keyTimes, double t, double c, 
  double b, Rcpp::NumericVector timesR, std::size_t nintertimes, bool closed
){
  std::vector<qtrn> keyRotors = _getRQuaternions(keyRotorsR);
  std::vector<double> times(timesR.begin(), timesR.end());
  keyRotors = _check_keyRotors(keyRotors, closed);
  const std::size_t nkeyRotors = keyRotors.size();
  const std::size_t nkeyTimes = keyTimes.size();
  if(nkeyTimes == 0 && nintertimes > 0){
    times = _seqvec(1.0, nkeyRotors, nintertimes * (nkeyRotors - 1) + 1);
    if(closed){
      times.pop_back();
    }
  }
  keyTimes = _check_keyTimes(keyTimes, nkeyRotors);
  std::vector<double> kTimes(keyTimes.begin(), keyTimes.end());
  std::vector<std::array<double, 3>> triplets_times = 
    makeTriplets_times(kTimes, closed);
  std::vector<std::array<qtrn, 3>> triplets_rotors = 
    makeTriplets_rotors(keyRotors, closed);
  
  Rcpp::Rcout << "TTTTTTTTTT  " << triplets_times[7][0] << " --- ";
  Rcpp::Rcout << "TTTTTTTTTT  " << triplets_times[7][1] << " --- ";
  Rcpp::Rcout << "TTTTTTTTTT  " << triplets_times[7][2] << " --- ";

  std::vector<qtrn> control_points(0);
  const std::size_t MMM = triplets_rotors.size(); // nkeyRotors-2 ?
  Rcpp::Rcout << "MMMMMMMMMMMMMMM " << MMM << "  " << triplets_times.size();
  for(std::size_t i = 0; i < MMM; i++){ 
    std::array<qtrn, 3> qs = triplets_rotors[i];
    std::array<qtrn, 2> qb_qa = _calculate_control_quaternions(
      qs, triplets_times[i], t, c, b
    );
    qtrn q_before = qb_qa[0];
    qtrn q_after  = qb_qa[1];
    Rcpp::Rcout << "UUUUUUUUUUU  " << q_before.norm() << " --- ";
    Rcpp::Rcout << "UUUUUUUUUUU  " << qs[1].norm() << " --- ";
    Rcpp::Rcout << "UUUUUUUUUUU  " << q_after.norm() << " --- ";
    control_points.push_back(q_before);
    control_points.push_back(qs[1]);
    control_points.push_back(qs[1]);
    control_points.push_back(q_after);
  }
  const std::size_t ncontrol_points = control_points.size();
  
  if(closed){
    // stopifnot(4*length(keyTimes) == n_control_points)
    Rcpp::Rcout << "XXXXXXXXXXXX" << ncontrol_points << " --- ";
    control_points.pop_back();
    control_points.pop_back();
    control_points.erase(control_points.begin());
    control_points.erase(control_points.begin());
  }else if(ncontrol_points == 0){
    // two quaternions -> slerp
    //stopifnot(n_keyRotors == 2L)
    //stopifnot(length(keyTimes) == 2L)
    qtrn q0 = keyRotors[0];
    qtrn q1 = keyRotors[1];
    qtrn offset = qpower(q1 * q0.inverse(), 1.0/3.0);
    control_points.push_back(q0);
    control_points.push_back(offset * q0);
    control_points.push_back(offset.inverse() * q1);
    control_points.push_back(q1);
  }else{ // natural
    qtrn ncq1 = _natural_control_quaternion(
      keyRotors[0], control_points[0], control_points[1]
    );
    qtrn ncq2 = _natural_control_quaternion(
      keyRotors[nkeyRotors-1], control_points[ncontrol_points-1], 
      control_points[ncontrol_points-2]
    );
    control_points.push_back(ncq2);
    control_points.push_back(keyRotors[nkeyRotors-1]);
    control_points.insert(control_points.begin(), ncq1);
    control_points.insert(control_points.begin(), keyRotors[0]);
  }

  const Rcpp::IntegerVector indices = 
    4 * Rcpp::seq_len(control_points.size() / 4) - 4;
  const std::size_t nsegments = indices.size();
  Rcpp::List Segments(nsegments);
  for(size_t i = 0; i < nsegments; i++) {
    const int j = indices(i);
        Rcpp::Rcout << "YYYYYXXXXXXXXXXXX" << control_points[j+3].norm() << " --- ";
    const std::vector<qtrn> segment = {
      control_points[j], control_points[j+1], 
      control_points[j+2], control_points[j+3]
    };
    Segments(i) = _getCQuaternions(segment);
  }
  Rcpp::NumericVector Rtimes(times.begin(), times.end());
  
  //keyTimes = Rcpp::NumericVector(0);
  return DeCasteljau_cpp(Segments, keyTimes, Rtimes);
}
// {} []
// // [[Rcpp::export]]


// // [[Rcpp::export]]
// void test(){
//   qtrn q1(1.0/sqrt(2),1.0/sqrt(2), 0.0, 0.0);
//   qtrn q2(0.0, 0.0, 0.0, 1.0);
//   qtrn q3(0.5, 0.5, 0.5, 0.5);
//   std::vector<qtrn> quats = {q1, q2, q3};
//   std::array<qtrn, 2> arr = _calculate_control_quaternions(quats, 0.2, 0.4, 0.6, 0.1, 0.1, 0.1);
//   qtrn qq1 = arr[0];
//   qtrn qq2 = arr[1];
//   Rcpp::Rcout << qq1.w() << " " << qq1.x() << " " << qq1.y() << " " << qq1.z() << " --- ";
//   Rcpp::Rcout << qq2.w() << " " << qq2.x() << " " << qq2.y() << " " << qq2.z() << " --- ";
// }

// {} []
// // [[Rcpp::export]]