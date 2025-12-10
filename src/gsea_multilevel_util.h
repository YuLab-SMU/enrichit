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
// Matches fgsea's score_t structure
class score_t {
public:
    int64_t NS;          // Normalization sum
    int64_t coef_NS;     // Coefficient for NS (numerator)
    int64_t diff;        // n - k (number of misses)
    int64_t coef_const;  // Coefficient for constant (number of misses before current position)
    
    score_t() : NS(0), coef_NS(0), diff(0), coef_const(0) {}

    score_t(int64_t NS_, int64_t coef_NS_, int64_t diff_, int64_t coef_const_)
        : NS(NS_), coef_NS(coef_NS_), diff(diff_), coef_const(coef_const_) {}
    
    // Get maximum NS value (use 2^30 to avoid overflow)
    static int64_t getMaxNS() {
        return (1LL << 30);
    }
    
    // Convert to double
    // score = coef_NS / NS - coef_const / diff
    double getDouble() const {
        if (NS == 0) return 0.0;
        double p_hit = static_cast<double>(coef_NS) / NS;
        double p_miss = (diff == 0) ? 0.0 : static_cast<double>(coef_const) / diff;
        return p_hit - p_miss;
    }
    
    // Get numerator (used for sign check)
    // numerator = coef_NS * diff - coef_const * NS
    int64_t getNumerator() const {
        return coef_NS * diff - coef_const * NS;
    }
    
    // Comparison using exact integer arithmetic (like fgsea)
    // Avoids floating point errors
    bool operator<(const score_t& other) const {
        // Compare: coef_NS/NS - coef_const/diff < other.coef_NS/other.NS - other.coef_const/other.diff
        // Rearrange to avoid division:
        // (coef_NS * diff - coef_const * NS) / (NS * diff) < (other.coef_NS * other.diff - other.coef_const * other.NS) / (other.NS * other.diff)
        // Cross multiply (assuming positive denominators):
        // (coef_NS * diff - coef_const * NS) * other.NS * other.diff < (other.coef_NS * other.diff - other.coef_const * other.NS) * NS * diff
        
        int64_t lhs_num = coef_NS * diff - coef_const * NS;
        int64_t rhs_num = other.coef_NS * other.diff - other.coef_const * other.NS;
        
        // Use long double for comparison to avoid overflow
        long double lhs = static_cast<long double>(lhs_num) / (static_cast<long double>(NS) * diff);
        long double rhs = static_cast<long double>(rhs_num) / (static_cast<long double>(other.NS) * other.diff);
        return lhs < rhs;
    }
    
    bool operator<=(const score_t& other) const {
        int64_t lhs_num = coef_NS * diff - coef_const * NS;
        int64_t rhs_num = other.coef_NS * other.diff - other.coef_const * other.NS;
        long double lhs = static_cast<long double>(lhs_num) / (static_cast<long double>(NS) * diff);
        long double rhs = static_cast<long double>(rhs_num) / (static_cast<long double>(other.NS) * other.diff);
        return lhs <= rhs;
    }
    
    bool operator>(const score_t& other) const {
        int64_t lhs_num = coef_NS * diff - coef_const * NS;
        int64_t rhs_num = other.coef_NS * other.diff - other.coef_const * other.NS;
        long double lhs = static_cast<long double>(lhs_num) / (static_cast<long double>(NS) * diff);
        long double rhs = static_cast<long double>(rhs_num) / (static_cast<long double>(other.NS) * other.diff);
        return lhs > rhs;
    }
    
    bool operator>=(const score_t& other) const {
        int64_t lhs_num = coef_NS * diff - coef_const * NS;
        int64_t rhs_num = other.coef_NS * other.diff - other.coef_const * other.NS;
        long double lhs = static_cast<long double>(lhs_num) / (static_cast<long double>(NS) * diff);
        long double rhs = static_cast<long double>(rhs_num) / (static_cast<long double>(other.NS) * other.diff);
        return lhs >= rhs;
    }
    
    bool operator==(const score_t& other) const {
        int64_t lhs_num = coef_NS * diff - coef_const * NS;
        int64_t rhs_num = other.coef_NS * other.diff - other.coef_const * other.NS;
        long double lhs = static_cast<long double>(lhs_num) / (static_cast<long double>(NS) * diff);
        long double rhs = static_cast<long double>(rhs_num) / (static_cast<long double>(other.NS) * other.diff);
        return std::abs(lhs - rhs) < 1e-15;
    }
    
    // Unary minus operator
    score_t operator-() const {
        return score_t(NS, -coef_NS, diff, -coef_const);
    }
    
    // Absolute value
    score_t abs() const {
        return std::max(*this, -(*this));
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
