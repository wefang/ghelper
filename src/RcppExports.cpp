// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// count_bins
NumericVector count_bins(NumericVector counts, NumericVector bins);
RcppExport SEXP _ghelper_count_bins(SEXP countsSEXP, SEXP binsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bins(binsSEXP);
    rcpp_result_gen = Rcpp::wrap(count_bins(counts, bins));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ghelper_count_bins", (DL_FUNC) &_ghelper_count_bins, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ghelper(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
