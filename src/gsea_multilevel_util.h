#ifndef GSEA_MULTILEVEL_UTIL_H
#define GSEA_MULTILEVEL_UTIL_H

#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <cstdint>

namespace enrichit {

// Random engine type
using random_engine_t = std::mt19937;

// Score class for precise ES calculation using integer arithmetic
class score_t {
private:
    int64_t NS;        // Normalization sum
    int64_t numerator; // Numerator
    int64_t nMiss;     // Number of misses
    int64_t nHit;      // Number of hits before current position
    
public:
    score_t() : NS(0), numerator(0), nMiss(0), nHit(0) {}

    score_t(int64_t NS_, int64_t num_, int64_t nMiss_, int64_t nHit_)
        : NS(NS_), numerator(num_), nMiss(nMiss_), nHit(nHit_) {}
    
    // Get maximum NS value (use 2^30 to avoid overflow)
    static int64_t getMaxNS() {
        return (1LL << 30);
    }
    
    // Convert to double
    double getDouble() const {
        if (NS == 0) return 0.0;
        double p_hit = static_cast<double>(numerator) / NS;
        double p_miss = (nMiss == 0) ? 0.0 : static_cast<double>(nHit) / nMiss;
        return p_hit - p_miss;
    }
    
    // Get numerator (used for sign check)
    int64_t getNumerator() const {
        return numerator;
    }
    
    // Comparison operators (for sorting and threshold checking)
    bool operator<(const score_t& other) const {
        long double lhs_phit = (NS == 0) ? 0.0L : (static_cast<long double>(numerator) / static_cast<long double>(NS));
        long double lhs_pmiss = (nMiss == 0) ? 0.0L : (static_cast<long double>(nHit) / static_cast<long double>(nMiss));
        long double rhs_phit = (other.NS == 0) ? 0.0L : (static_cast<long double>(other.numerator) / static_cast<long double>(other.NS));
        long double rhs_pmiss = (other.nMiss == 0) ? 0.0L : (static_cast<long double>(other.nHit) / static_cast<long double>(other.nMiss));
        return (lhs_phit - lhs_pmiss) < (rhs_phit - rhs_pmiss);
    }
    
    bool operator<=(const score_t& other) const {
        long double lhs_phit = (NS == 0) ? 0.0L : (static_cast<long double>(numerator) / static_cast<long double>(NS));
        long double lhs_pmiss = (nMiss == 0) ? 0.0L : (static_cast<long double>(nHit) / static_cast<long double>(nMiss));
        long double rhs_phit = (other.NS == 0) ? 0.0L : (static_cast<long double>(other.numerator) / static_cast<long double>(other.NS));
        long double rhs_pmiss = (other.nMiss == 0) ? 0.0L : (static_cast<long double>(other.nHit) / static_cast<long double>(other.nMiss));
        return (lhs_phit - lhs_pmiss) <= (rhs_phit - rhs_pmiss);
    }
    
    bool operator>(const score_t& other) const {
        long double lhs_phit = (NS == 0) ? 0.0L : (static_cast<long double>(numerator) / static_cast<long double>(NS));
        long double lhs_pmiss = (nMiss == 0) ? 0.0L : (static_cast<long double>(nHit) / static_cast<long double>(nMiss));
        long double rhs_phit = (other.NS == 0) ? 0.0L : (static_cast<long double>(other.numerator) / static_cast<long double>(other.NS));
        long double rhs_pmiss = (other.nMiss == 0) ? 0.0L : (static_cast<long double>(other.nHit) / static_cast<long double>(other.nMiss));
        return (lhs_phit - lhs_pmiss) > (rhs_phit - rhs_pmiss);
    }
    
    bool operator>=(const score_t& other) const {
        long double lhs_phit = (NS == 0) ? 0.0L : (static_cast<long double>(numerator) / static_cast<long double>(NS));
        long double lhs_pmiss = (nMiss == 0) ? 0.0L : (static_cast<long double>(nHit) / static_cast<long double>(nMiss));
        long double rhs_phit = (other.NS == 0) ? 0.0L : (static_cast<long double>(other.numerator) / static_cast<long double>(other.NS));
        long double rhs_pmiss = (other.nMiss == 0) ? 0.0L : (static_cast<long double>(other.nHit) / static_cast<long double>(other.nMiss));
        return (lhs_phit - lhs_pmiss) >= (rhs_phit - rhs_pmiss);
    }
    
    bool operator==(const score_t& other) const {
        long double lhs_phit = (NS == 0) ? 0.0L : (static_cast<long double>(numerator) / static_cast<long double>(NS));
        long double lhs_pmiss = (nMiss == 0) ? 0.0L : (static_cast<long double>(nHit) / static_cast<long double>(nMiss));
        long double rhs_phit = (other.NS == 0) ? 0.0L : (static_cast<long double>(other.numerator) / static_cast<long double>(other.NS));
        long double rhs_pmiss = (other.nMiss == 0) ? 0.0L : (static_cast<long double>(other.nHit) / static_cast<long double>(other.nMiss));
        return (lhs_phit - lhs_pmiss) == (rhs_phit - rhs_pmiss);
    }
};

// Calculate ES for a given sample (returns score_t for precise comparison)
score_t calcES(const std::vector<int64_t>& ranks, const std::vector<int>& sample);

// Calculate positive ES (absolute value, for sorting)
score_t calcPositiveES(const std::vector<int64_t>& ranks, const std::vector<int>& sample);

// Beta mean log: E[log(Beta(a, b))] = digamma(a) - digamma(a+b)
double betaMeanLog(unsigned long a, unsigned long b);

// Generate random combination of k indices from [a, b]
std::vector<int> combination(int a, int b, int k, random_engine_t& rng);

// Uniform integer distribution wrapper for efficiency
struct uid_wrapper {
    int from, len;
    random_engine_t& rng;
    
    uid_wrapper(int from_, int to_, random_engine_t& rng_)
        : from(from_), len(to_ - from_ + 1), rng(rng_) {}
    
    int operator()();
};

// Variance per level for error estimation
double getVarPerLevel(unsigned long k, unsigned long n);

} // namespace enrichit

#endif // GSEA_MULTILEVEL_UTIL_H
