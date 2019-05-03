// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// overwrite_as_adjacency_cpp
NumericMatrix overwrite_as_adjacency_cpp(NumericMatrix m, double tol);
RcppExport SEXP _SeqNet_overwrite_as_adjacency_cpp(SEXP mSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(overwrite_as_adjacency_cpp(m, tol));
    return rcpp_result_gen;
END_RCPP
}
// check_adjacency_cpp
bool check_adjacency_cpp(NumericMatrix m);
RcppExport SEXP _SeqNet_check_adjacency_cpp(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(check_adjacency_cpp(m));
    return rcpp_result_gen;
END_RCPP
}
// components_in_adjacency
IntegerVector components_in_adjacency(NumericMatrix adj);
RcppExport SEXP _SeqNet_components_in_adjacency(SEXP adjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type adj(adjSEXP);
    rcpp_result_gen = Rcpp::wrap(components_in_adjacency(adj));
    return rcpp_result_gen;
END_RCPP
}
// ecdf_cpp
NumericVector ecdf_cpp(NumericVector x);
RcppExport SEXP _SeqNet_ecdf_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(ecdf_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// edges_from_adjacency_cpp
NumericMatrix edges_from_adjacency_cpp(NumericMatrix adj);
RcppExport SEXP _SeqNet_edges_from_adjacency_cpp(SEXP adjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type adj(adjSEXP);
    rcpp_result_gen = Rcpp::wrap(edges_from_adjacency_cpp(adj));
    return rcpp_result_gen;
END_RCPP
}
// ring_lattice_cpp
NumericMatrix ring_lattice_cpp(int p, int neig_size);
RcppExport SEXP _SeqNet_ring_lattice_cpp(SEXP pSEXP, SEXP neig_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type neig_size(neig_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(ring_lattice_cpp(p, neig_size));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SeqNet_overwrite_as_adjacency_cpp", (DL_FUNC) &_SeqNet_overwrite_as_adjacency_cpp, 2},
    {"_SeqNet_check_adjacency_cpp", (DL_FUNC) &_SeqNet_check_adjacency_cpp, 1},
    {"_SeqNet_components_in_adjacency", (DL_FUNC) &_SeqNet_components_in_adjacency, 1},
    {"_SeqNet_ecdf_cpp", (DL_FUNC) &_SeqNet_ecdf_cpp, 1},
    {"_SeqNet_edges_from_adjacency_cpp", (DL_FUNC) &_SeqNet_edges_from_adjacency_cpp, 1},
    {"_SeqNet_ring_lattice_cpp", (DL_FUNC) &_SeqNet_ring_lattice_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_SeqNet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
