// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pcoriaccel_hello
[[nodiscard]] String pcoriaccel_hello();
RcppExport SEXP _pcoriRPackage_pcoriaccel_hello() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(pcoriaccel_hello());
    return rcpp_result_gen;
END_RCPP
}
// pcoriaccel_NW_basic
[[nodiscard]] NumericMatrix pcoriaccel_NW_basic(NumericVector Xb, NumericVector Y, NumericVector xb, NumericVector y_seq, double h);
RcppExport SEXP _pcoriRPackage_pcoriaccel_NW_basic(SEXP XbSEXP, SEXP YSEXP, SEXP xbSEXP, SEXP y_seqSEXP, SEXP hSEXP) {
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
// pcoriaccel_estimate_pmf
[[nodiscard]] NumericVector pcoriaccel_estimate_pmf(NumericVector X, NumericVector Y, double xi, NumericVector y_seq, double h);
RcppExport SEXP _pcoriRPackage_pcoriaccel_estimate_pmf(SEXP XSEXP, SEXP YSEXP, SEXP xiSEXP, SEXP y_seqSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y_seq(y_seqSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(pcoriaccel_estimate_pmf(X, Y, xi, y_seq, h));
    return rcpp_result_gen;
END_RCPP
}
// pcoriaccel_NW
[[nodiscard]] NumericMatrix pcoriaccel_NW(NumericVector Xb, NumericVector Y, NumericVector xb, NumericVector y_seq, double h, String kernel);
RcppExport SEXP _pcoriRPackage_pcoriaccel_NW(SEXP XbSEXP, SEXP YSEXP, SEXP xbSEXP, SEXP y_seqSEXP, SEXP hSEXP, SEXP kernelSEXP) {
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

static const R_CallMethodDef CallEntries[] = {
    {"_pcoriRPackage_pcoriaccel_hello", (DL_FUNC) &_pcoriRPackage_pcoriaccel_hello, 0},
    {"_pcoriRPackage_pcoriaccel_NW_basic", (DL_FUNC) &_pcoriRPackage_pcoriaccel_NW_basic, 5},
    {"_pcoriRPackage_pcoriaccel_estimate_pmf", (DL_FUNC) &_pcoriRPackage_pcoriaccel_estimate_pmf, 5},
    {"_pcoriRPackage_pcoriaccel_NW", (DL_FUNC) &_pcoriRPackage_pcoriaccel_NW, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_pcoriRPackage(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
