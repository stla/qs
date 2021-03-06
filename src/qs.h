#ifndef _QSHEADER_
#define _QSHEADER_

#include <RcppEigen.h>
#include <Eigen/Geometry> 
// #include "quaternion.h"

typedef Eigen::Quaterniond qtrn;



std::size_t _check_time(double, Rcpp::NumericVector, bool);

Rcpp::NumericVector _seq_len(std::size_t);

qtrn slerp(qtrn, qtrn, double);

qtrn qexp(qtrn);
qtrn qlog(qtrn);
qtrn qpower(qtrn, double);

std::vector<double> _seqvec(double, double, std::size_t);

Rcpp::NumericVector _check_keyTimes(Rcpp::NumericVector, std::size_t);
std::vector<qtrn> _check_keyRotors(std::vector<qtrn>, bool);

#endif
