/*
 * This file is based on the fgsea package (https://github.com/ctlab/fgsea)
 * Copyright (c) 2016-2024 Alexey Sergushichev
 *
 * It has been adapted for the enrichit package.
 */

#include "gsea_multilevel_util.h"
#include <Rcpp.h>

namespace enrichit {

uid_wrapper::uid_wrapper(int min, int max, random_engine_t& rng)
    : dist(min, max), rng(rng) {}

int uid_wrapper::operator()() {
    return dist(rng);
}

std::vector<int> combination(const int &a, const int &b, const int &k, random_engine_t& rng) {
    uid_wrapper uni(a, b, rng);
    std::vector<int> v;
    v.reserve(k);

    int n = b - a + 1;
    std::vector<char> used(n, 0); // using char as bool
    
    // If k is small, use rejection sampling
    if (k < n * 1.0 / 2) {
        for (int i = 0; i < k; i++) {
            while (true) {
                int x = uni();
                if (!used[x - a]) {
                    v.push_back(x);
                    used[x - a] = 1;
                    break;
                }
            }
        }
    } else {
        // Floyd's algorithm variant
        for (int r = n - k; r < n; ++r){
            int x = std::uniform_int_distribution<int>(0, r)(rng);
            if (!used[x]){
                v.push_back(a + x);
                used[x] = 1;
            } else{
                v.push_back(a + r);
                used[r] = 1;
            }
        }
        
        // Shuffle result
        for (int i = v.size() - 1; i > 0; --i) {
            int j = std::uniform_int_distribution<int>(0, i)(rng);
            std::swap(v[i], v[j]);
        }
    }

    return v;
}

double betaMeanLog(unsigned long k, unsigned long n) {
    return R::digamma((double)k) - R::digamma((double)(n + 1));
}

double getVarPerLevel(unsigned long k, unsigned long n) {
    return R::trigamma((double)k) - R::trigamma((double)(n + 1));
}

} // namespace enrichit
