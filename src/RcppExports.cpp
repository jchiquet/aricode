// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// countPairs
List countPairs(IntegerVector classi1, IntegerVector classi2, IntegerVector order);
RcppExport SEXP _aricode_countPairs(SEXP classi1SEXP, SEXP classi2SEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type classi1(classi1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type classi2(classi2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(countPairs(classi1, classi2, order));
    return rcpp_result_gen;
END_RCPP
}
// expected_MI
double expected_MI(IntegerVector ni_, IntegerVector n_j);
RcppExport SEXP _aricode_expected_MI(SEXP ni_SEXP, SEXP n_jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ni_(ni_SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_j(n_jSEXP);
    rcpp_result_gen = Rcpp::wrap(expected_MI(ni_, n_j));
    return rcpp_result_gen;
END_RCPP
}
// getRank
List getRank(IntegerVector classi);
RcppExport SEXP _aricode_getRank(SEXP classiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type classi(classiSEXP);
    rcpp_result_gen = Rcpp::wrap(getRank(classi));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_aricode_countPairs", (DL_FUNC) &_aricode_countPairs, 3},
    {"_aricode_expected_MI", (DL_FUNC) &_aricode_expected_MI, 2},
    {"_aricode_getRank", (DL_FUNC) &_aricode_getRank, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_aricode(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
