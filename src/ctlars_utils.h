// header file (.h) for ctlars_utils
// File stores C++ functions for direct use in R.
// ----------------------------------------------------------------------------
#ifndef ctlars_utils_H
#define ctlars_utils_H

// ----------------------------------------------------------------------------
#include <string.h>
#include "Rcpp.h"

//' @title
//' Hello ctlars from your C++ friends :).
//'
//' @name
//' hello_ctlars_from_Cpp
//'
//' @description
//' This function returns a greeting message.
//'
//' @return A string with "Hello, ctlars!".
//'
//' @export
// -------------------------------------------------------------------------
// [[Rcpp::export]]
std::string hello_ctlars_from_Cpp();


// ----------------------------------------------------------------------------
#endif
