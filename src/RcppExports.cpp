// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pcoriaccel_NW_basic
[[nodiscard]] NumericMatrix pcoriaccel_NW_basic(NumericVector Xb, NumericVector Y, NumericVector xb, NumericVector y_seq, double h);
RcppExport SEXP _SensIAT_pcoriaccel_NW_basic(SEXP XbSEXP, SEXP YSEXP, SEXP xbSEXP, SEXP y_seqSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Xb(XbSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xb(xbSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y_seq(y_seqSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(pcoriaccel_NW_basic(Xb, Y, xb, y_seq, h));
    return rcpp_result_gen;
END_RCPP
}
// pcoriaccel_NW
[[nodiscard]] NumericMatrix pcoriaccel_NW(NumericVector Xb, NumericVector Y, NumericVector xb, NumericVector y_seq, double h, String kernel);
RcppExport SEXP _SensIAT_pcoriaccel_NW(SEXP XbSEXP, SEXP YSEXP, SEXP xbSEXP, SEXP y_seqSEXP, SEXP hSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Xb(XbSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xb(xbSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y_seq(y_seqSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< String >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(pcoriaccel_NW(Xb, Y, xb, y_seq, h, kernel));
    return rcpp_result_gen;
END_RCPP
}
// pcoriaccel_mmul
[[nodiscard]] SEXP pcoriaccel_mmul(SEXP matrA, SEXP matrB);
RcppExport SEXP _SensIAT_pcoriaccel_mmul(SEXP matrASEXP, SEXP matrBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type matrA(matrASEXP);
    Rcpp::traits::input_parameter< SEXP >::type matrB(matrBSEXP);
    rcpp_result_gen = Rcpp::wrap(pcoriaccel_mmul(matrA, matrB));
    return rcpp_result_gen;
END_RCPP
}
// pcoriaccel_inner_prod
[[nodiscard]] double pcoriaccel_inner_prod(NumericVector vecA, NumericVector vecB);
RcppExport SEXP _SensIAT_pcoriaccel_inner_prod(SEXP vecASEXP, SEXP vecBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vecA(vecASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vecB(vecBSEXP);
    rcpp_result_gen = Rcpp::wrap(pcoriaccel_inner_prod(vecA, vecB));
    return rcpp_result_gen;
END_RCPP
}
// pcoriaccel_outer_sum
[[nodiscard]] NumericMatrix pcoriaccel_outer_sum(NumericVector vecA, NumericVector vecB);
RcppExport SEXP _SensIAT_pcoriaccel_outer_sum(SEXP vecASEXP, SEXP vecBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vecA(vecASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vecB(vecBSEXP);
    rcpp_result_gen = Rcpp::wrap(pcoriaccel_outer_sum(vecA, vecB));
    return rcpp_result_gen;
END_RCPP
}
// pcoriaccel_outer_prod
[[nodiscard]] NumericMatrix pcoriaccel_outer_prod(NumericVector vecA, NumericVector vecB);
RcppExport SEXP _SensIAT_pcoriaccel_outer_prod(SEXP vecASEXP, SEXP vecBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vecA(vecASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vecB(vecBSEXP);
    rcpp_result_gen = Rcpp::wrap(pcoriaccel_outer_prod(vecA, vecB));
    return rcpp_result_gen;
END_RCPP
}
// pcoriaccel_sorted_unique
[[nodiscard]] NumericVector pcoriaccel_sorted_unique(NumericVector vec);
RcppExport SEXP _SensIAT_pcoriaccel_sorted_unique(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(pcoriaccel_sorted_unique(vec));
    return rcpp_result_gen;
END_RCPP
}
// pcoriaccel_compute_influence_term_2_quadv_sim_via_matrix
[[nodiscard]] NumericMatrix pcoriaccel_compute_influence_term_2_quadv_sim_via_matrix(NumericMatrix X, NumericVector Y, NumericVector times, NumericMatrix individual_X, NumericVector x_slope, NumericVector alpha, NumericVector beta, S4 spline_basis, double bandwidth, double tol, String kernel);
RcppExport SEXP _SensIAT_pcoriaccel_compute_influence_term_2_quadv_sim_via_matrix(SEXP XSEXP, SEXP YSEXP, SEXP timesSEXP, SEXP individual_XSEXP, SEXP x_slopeSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP spline_basisSEXP, SEXP bandwidthSEXP, SEXP tolSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type times(timesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type individual_X(individual_XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_slope(x_slopeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< S4 >::type spline_basis(spline_basisSEXP);
    Rcpp::traits::input_parameter< double >::type bandwidth(bandwidthSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< String >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(pcoriaccel_compute_influence_term_2_quadv_sim_via_matrix(X, Y, times, individual_X, x_slope, alpha, beta, spline_basis, bandwidth, tol, kernel));
    return rcpp_result_gen;
END_RCPP
}
// pcoriaccel_estimate_pmf
[[nodiscard]] NumericVector pcoriaccel_estimate_pmf(NumericVector Xb, NumericVector Y, double xi, NumericVector y_seq, double h, String kernel);
RcppExport SEXP _SensIAT_pcoriaccel_estimate_pmf(SEXP XbSEXP, SEXP YSEXP, SEXP xiSEXP, SEXP y_seqSEXP, SEXP hSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Xb(XbSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y_seq(y_seqSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< String >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(pcoriaccel_estimate_pmf(Xb, Y, xi, y_seq, h, kernel));
    return rcpp_result_gen;
END_RCPP
}
// pcoriaccel_integrate_simp
[[nodiscard]] List pcoriaccel_integrate_simp(Function integrand, double lo, double hi, double tol);
RcppExport SEXP _SensIAT_pcoriaccel_integrate_simp(SEXP integrandSEXP, SEXP loSEXP, SEXP hiSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Function >::type integrand(integrandSEXP);
    Rcpp::traits::input_parameter< double >::type lo(loSEXP);
    Rcpp::traits::input_parameter< double >::type hi(hiSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(pcoriaccel_integrate_simp(integrand, lo, hi, tol));
    return rcpp_result_gen;
END_RCPP
}
// pcoriaccel_evaluate_basis
[[nodiscard]] NumericVector pcoriaccel_evaluate_basis(S4 spline_basis, double x);
RcppExport SEXP _SensIAT_pcoriaccel_evaluate_basis(SEXP spline_basisSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type spline_basis(spline_basisSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(pcoriaccel_evaluate_basis(spline_basis, x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SensIAT_pcoriaccel_NW_basic", (DL_FUNC) &_SensIAT_pcoriaccel_NW_basic, 5},
    {"_SensIAT_pcoriaccel_NW", (DL_FUNC) &_SensIAT_pcoriaccel_NW, 6},
    {"_SensIAT_pcoriaccel_mmul", (DL_FUNC) &_SensIAT_pcoriaccel_mmul, 2},
    {"_SensIAT_pcoriaccel_inner_prod", (DL_FUNC) &_SensIAT_pcoriaccel_inner_prod, 2},
    {"_SensIAT_pcoriaccel_outer_sum", (DL_FUNC) &_SensIAT_pcoriaccel_outer_sum, 2},
    {"_SensIAT_pcoriaccel_outer_prod", (DL_FUNC) &_SensIAT_pcoriaccel_outer_prod, 2},
    {"_SensIAT_pcoriaccel_sorted_unique", (DL_FUNC) &_SensIAT_pcoriaccel_sorted_unique, 1},
    {"_SensIAT_pcoriaccel_compute_influence_term_2_quadv_sim_via_matrix", (DL_FUNC) &_SensIAT_pcoriaccel_compute_influence_term_2_quadv_sim_via_matrix, 11},
    {"_SensIAT_pcoriaccel_estimate_pmf", (DL_FUNC) &_SensIAT_pcoriaccel_estimate_pmf, 6},
    {"_SensIAT_pcoriaccel_integrate_simp", (DL_FUNC) &_SensIAT_pcoriaccel_integrate_simp, 4},
    {"_SensIAT_pcoriaccel_evaluate_basis", (DL_FUNC) &_SensIAT_pcoriaccel_evaluate_basis, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_SensIAT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
