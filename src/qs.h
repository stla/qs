#ifndef _QSHEADER_
#define _QSHEADER_

#include <RcppEigen.h>
#include <Eigen/Geometry> 

typedef Eigen::Quaterniond qtrn;



std::size_t _check_time(double, Rcpp::NumericVector, bool);

#endif
