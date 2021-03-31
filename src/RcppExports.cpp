// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// bresenham
IntegerVector bresenham(int x1, int y1, int x2, int y2);
RcppExport SEXP _GVI_bresenham(SEXP x1SEXP, SEXP y1SEXP, SEXP x2SEXP, SEXP y2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< int >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< int >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< int >::type y2(y2SEXP);
    rcpp_result_gen = Rcpp::wrap(bresenham(x1, y1, x2, y2));
    return rcpp_result_gen;
END_RCPP
}
// isVisible
LogicalVector isVisible(NumericVector x);
RcppExport SEXP _GVI_isVisible(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(isVisible(x));
    return rcpp_result_gen;
END_RCPP
}
// tangents
NumericVector tangents(int x1, int y1, double height0, IntegerVector xy2, NumericVector dsm_profile);
RcppExport SEXP _GVI_tangents(SEXP x1SEXP, SEXP y1SEXP, SEXP height0SEXP, SEXP xy2SEXP, SEXP dsm_profileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< int >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< double >::type height0(height0SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type xy2(xy2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dsm_profile(dsm_profileSEXP);
    rcpp_result_gen = Rcpp::wrap(tangents(x1, y1, height0, xy2, dsm_profile));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GVI_bresenham", (DL_FUNC) &_GVI_bresenham, 4},
    {"_GVI_isVisible", (DL_FUNC) &_GVI_isVisible, 1},
    {"_GVI_tangents", (DL_FUNC) &_GVI_tangents, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_GVI(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
