#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector count_bins(NumericVector counts, NumericVector bins){
  int i;
  for (i=0; i<bins.size(); i++){ counts[bins[i] - 1]++; }
  return counts;
}


