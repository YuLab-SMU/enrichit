#ifndef ENRICHIT_H
#define ENRICHIT_H

#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>

namespace enrichit {

// ORA function declaration
Rcpp::DataFrame ora(const Rcpp::CharacterVector& gene,
                    const Rcpp::CharacterVector& universe,
                    const Rcpp::List& gene_sets);

// GSEA function declaration
Rcpp::DataFrame gsea(const Rcpp::NumericVector& stats,
                     const Rcpp::List& gene_sets,
                     const Rcpp::CharacterVector& gene_set_names,
                     int nPerm,
                     double exponent,
                     std::string method);

// Adaptive GSEA function declaration
Rcpp::DataFrame gsea_adaptive(const Rcpp::NumericVector& stats,
                              const Rcpp::List& gene_sets,
                              const Rcpp::CharacterVector& gene_set_names,
                              int minPerm,
                              int maxPerm,
                              double pvalThreshold,
                              double exponent,
                              std::string method);

// Helper function for hypergeometric distribution (log scale)
double dhyper(int k, int m, int n, int k_plus_m_minus_n, bool log_p = false);

} // namespace enrichit

#endif // ENRICHIT_H
