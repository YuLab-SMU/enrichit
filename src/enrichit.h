#ifndef ENRICHIT_H
#define ENRICHIT_H

#include <Rcpp.h>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <algorithm>

namespace enrichit {

// Hypergeometric probability mass function
inline double dhyper(int x, int m, int n, int k) {
    if (x < std::max(0, k - n) || x > std::min(m, k)) {
        return 0.0;
    }
    // log combination: lchoose(n, k) = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
    double log_prob = (lgamma(m + 1) - lgamma(x + 1) - lgamma(m - x + 1)) +
                      (lgamma(n + 1) - lgamma(k - x + 1) - lgamma(n - k + x + 1)) -
                      (lgamma(m + n + 1) - lgamma(k + 1) - lgamma(m + n - k + 1));
    return std::exp(log_prob);
}



// Main ORA function
Rcpp::DataFrame ora(const Rcpp::CharacterVector& gene_set,
                    const Rcpp::CharacterVector& background,
                    const Rcpp::List& gene_sets,
                    const Rcpp::CharacterVector& gene_set_names);

} // namespace enrichit

#endif // ENRICHIT_H
