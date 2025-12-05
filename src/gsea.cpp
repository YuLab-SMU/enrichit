#include "enrichit.h"
#include <random>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace enrichit {

// Structure to hold ES and leading edge info
struct GSEA_Result {
    double ES;
    int rank;
    double tags;
    double list;
    double signal;
    std::string core_enrichment;
};

// Helper to calculate Enrichment Score (ES) and leading edge info
GSEA_Result calculate_es_details(const std::vector<double>& stats,
                                 const std::vector<bool>& in_set,
                                 const Rcpp::CharacterVector& gene_names,
                                 double exponent) {
    
    double N_R = 0.0; // Sum of weights for genes in set
    int N_H = 0;      // Number of hits (genes in set)
    int N = stats.size();
    
    for (int i = 0; i < N; ++i) {
        if (in_set[i]) {
            N_R += std::pow(std::abs(stats[i]), exponent);
            N_H++;
        }
    }
    
    GSEA_Result res;
    res.ES = 0.0;
    res.rank = 0;
    res.tags = 0.0;
    res.list = 0.0;
    res.signal = 0.0;
    res.core_enrichment = "";
    
    if (N_H == 0) return res;
    
    double P_hit = 0.0;
    double P_miss = 0.0;
    double max_dev = 0.0;
    int peak_idx = 0;
    
    double N_miss = (double)(N - N_H);
    
    for (int i = 0; i < N; ++i) {
        if (in_set[i]) {
            P_hit += std::pow(std::abs(stats[i]), exponent) / N_R;
        } else {
            P_miss += 1.0 / N_miss;
        }
        
        double dev = P_hit - P_miss;
        if (std::abs(dev) > std::abs(max_dev)) {
            max_dev = dev;
            peak_idx = i;
        }
    }
    
    res.ES = max_dev;
    
    if (res.ES >= 0) {
        res.rank = peak_idx + 1;
    } else {
        res.rank = N - peak_idx;
    }
    
    int hits_in_ledge = 0;
    std::vector<std::string> core_genes;
    
    if (res.ES >= 0) {
        for (int i = 0; i <= peak_idx; ++i) {
            if (in_set[i]) {
                hits_in_ledge++;
                core_genes.push_back(Rcpp::as<std::string>(gene_names[i]));
            }
        }
        res.tags = (double)hits_in_ledge / N_H;
        res.list = (double)(peak_idx + 1) / N;
    } else {
        for (int i = peak_idx; i < N; ++i) {
            if (in_set[i]) {
                hits_in_ledge++;
                core_genes.push_back(Rcpp::as<std::string>(gene_names[i]));
            }
        }
        res.tags = (double)hits_in_ledge / N_H;
        res.list = (double)(N - peak_idx) / N;
    }
    
    res.signal = res.tags * (1.0 - res.list) * ((double)N / (N - N_H));
    
    std::stringstream ss;
    for (size_t i = 0; i < core_genes.size(); ++i) {
        if (i > 0) ss << "/";
        ss << core_genes[i];
    }
    res.core_enrichment = ss.str();
    
    return res;
}

// Simplified ES calculation for "permute" method (label permutation)
double calculate_es_permute(const std::vector<double>& stats,
                            const std::vector<bool>& in_set,
                            const std::vector<int>& perm_idx,
                            double exponent) {
    
    double N_R = 0.0;
    int N_H = 0;
    int N = stats.size();
    
    for (int i = 0; i < N; ++i) {
        if (in_set[perm_idx[i]]) {
            N_R += std::pow(std::abs(stats[i]), exponent);
            N_H++;
        }
    }
    
    if (N_H == 0) return 0.0;
    
    double P_hit = 0.0;
    double P_miss = 0.0;
    double max_dev = 0.0;
    double N_miss = (double)(N - N_H);
    
    for (int i = 0; i < N; ++i) {
        if (in_set[perm_idx[i]]) {
            P_hit += std::pow(std::abs(stats[i]), exponent) / N_R;
        } else {
            P_miss += 1.0 / N_miss;
        }
        double dev = P_hit - P_miss;
        if (std::abs(dev) > std::abs(max_dev)) {
            max_dev = dev;
        }
    }
    return max_dev;
}

// Optimized ES calculation for "sample" method (random sampling)
// Takes sorted indices of hits
double calculate_es_sparse(const std::vector<double>& stats,
                           const std::vector<int>& hit_indices,
                           double exponent) {
    
    int N = stats.size();
    int N_H = hit_indices.size();
    
    if (N_H == 0) return 0.0;
    
    double N_R = 0.0;
    for (int idx : hit_indices) {
        N_R += std::pow(std::abs(stats[idx]), exponent);
    }
    
    double P_hit = 0.0;
    double P_miss = 0.0;
    double max_dev = 0.0;
    double N_miss = (double)(N - N_H);
    double dec_per_miss = 1.0 / N_miss;
    
    int last_idx = -1;
    
    for (int idx : hit_indices) {
        // Process misses before this hit
        int misses = idx - last_idx - 1;
        if (misses > 0) {
            P_miss += misses * dec_per_miss;
            double dev = P_hit - P_miss;
            if (std::abs(dev) > std::abs(max_dev)) max_dev = dev;
        }
        
        // Process hit
        P_hit += std::pow(std::abs(stats[idx]), exponent) / N_R;
        double dev = P_hit - P_miss;
        if (std::abs(dev) > std::abs(max_dev)) max_dev = dev;
        
        last_idx = idx;
    }
    
    // Process remaining misses
    int remaining_misses = N - 1 - last_idx;
    if (remaining_misses > 0) {
        P_miss += remaining_misses * dec_per_miss;
        double dev = P_hit - P_miss;
        if (std::abs(dev) > std::abs(max_dev)) max_dev = dev;
    }
    
    return max_dev;
}

Rcpp::DataFrame gsea(const Rcpp::NumericVector& stats,
                     const Rcpp::List& gene_sets,
                     const Rcpp::CharacterVector& gene_set_names,
                     int nPerm,
                     double exponent,
                     std::string method) {
    
    int n_sets = gene_sets.size();
    int n_genes = stats.size();
    
    std::vector<double> gene_stats = Rcpp::as<std::vector<double>>(stats);
    Rcpp::CharacterVector gene_names = stats.names();
    
    std::unordered_map<std::string, int> gene_map;
    for (int i = 0; i < n_genes; ++i) {
        gene_map[Rcpp::as<std::string>(gene_names[i])] = i;
    }
    
    std::vector<std::vector<bool>> gs_bools(n_sets, std::vector<bool>(n_genes, false));
    std::vector<int> gs_sizes(n_sets, 0);
    
    for (int i = 0; i < n_sets; ++i) {
        Rcpp::CharacterVector gs = gene_sets[i];
        int count = 0;
        for (int j = 0; j < gs.size(); ++j) {
            std::string g = Rcpp::as<std::string>(gs[j]);
            if (gene_map.find(g) != gene_map.end()) {
                gs_bools[i][gene_map[g]] = true;
                count++;
            }
        }
        gs_sizes[i] = count;
    }
    
    Rcpp::NumericVector es(n_sets);
    Rcpp::IntegerVector rank(n_sets);
    Rcpp::NumericVector tags(n_sets);
    Rcpp::NumericVector list(n_sets);
    Rcpp::NumericVector signal(n_sets);
    Rcpp::CharacterVector core_enrichment(n_sets);
    
    for (int i = 0; i < n_sets; ++i) {
        if (gs_sizes[i] > 0) {
            GSEA_Result res = calculate_es_details(gene_stats, gs_bools[i], gene_names, exponent);
            es[i] = res.ES;
            rank[i] = res.rank;
            tags[i] = res.tags;
            list[i] = res.list;
            signal[i] = res.signal;
            core_enrichment[i] = res.core_enrichment;
        } else {
            es[i] = 0.0;
            rank[i] = 0;
            tags[i] = 0.0;
            list[i] = 0.0;
            signal[i] = 0.0;
            core_enrichment[i] = "";
        }
    }
    
    std::vector<std::vector<double>> perm_es(n_sets, std::vector<double>(nPerm));
    std::mt19937 rng(12345);
    
    if (method == "permute") {
        // Label permutation: shuffle indices once per permutation
        std::vector<int> perm_idx(n_genes);
        std::iota(perm_idx.begin(), perm_idx.end(), 0);
        
        for (int p = 0; p < nPerm; ++p) {
            std::shuffle(perm_idx.begin(), perm_idx.end(), rng);
            for (int i = 0; i < n_sets; ++i) {
                if (gs_sizes[i] > 0) {
                    perm_es[i][p] = calculate_es_permute(gene_stats, gs_bools[i], perm_idx, exponent);
                } else {
                    perm_es[i][p] = 0.0;
                }
            }
        }
    } else {
        // Random sampling: sample K positions for each set
        std::vector<int> universe(n_genes);
        std::iota(universe.begin(), universe.end(), 0);
        
        for (int i = 0; i < n_sets; ++i) {
            int k = gs_sizes[i];
            if (k == 0) {
                for (int p = 0; p < nPerm; ++p) perm_es[i][p] = 0.0;
                continue;
            }
            
            std::vector<int> sample(k);
            for (int p = 0; p < nPerm; ++p) {
                // Partial shuffle to get k random indices
                for (int j = 0; j < k; ++j) {
                    std::uniform_int_distribution<> dis(j, n_genes - 1);
                    int swap_idx = dis(rng);
                    std::swap(universe[j], universe[swap_idx]);
                    sample[j] = universe[j];
                }
                std::sort(sample.begin(), sample.end());
                perm_es[i][p] = calculate_es_sparse(gene_stats, sample, exponent);
            }
        }
    }
    
    Rcpp::NumericVector pvalues(n_sets);
    Rcpp::NumericVector nes(n_sets);
    
    for (int i = 0; i < n_sets; ++i) {
        double obs_es = es[i];
        
        if (gs_sizes[i] == 0) {
            pvalues[i] = 1.0;
            nes[i] = 0.0;
            continue;
        }
        
        int count_better = 0;
        double sum_pos_es = 0.0;
        double sum_neg_es = 0.0;
        int count_pos = 0;
        int count_neg = 0;
        
        for (int p = 0; p < nPerm; ++p) {
            double pes = perm_es[i][p];
            if (obs_es > 0) {
                if (pes >= obs_es) count_better++;
                if (pes >= 0) { sum_pos_es += pes; count_pos++; }
            } else {
                if (pes <= obs_es) count_better++;
                if (pes < 0) { sum_neg_es += pes; count_neg++; }
            }
        }
        
        pvalues[i] = (double)(count_better + 1) / (double)(nPerm + 1);
        
        if (obs_es > 0) {
            double mean_pos = (count_pos > 0) ? (sum_pos_es / count_pos) : 1.0;
            nes[i] = obs_es / mean_pos;
        } else if (obs_es < 0) {
            double mean_neg = (count_neg > 0) ? (sum_neg_es / count_neg) : -1.0;
            nes[i] = obs_es / std::abs(mean_neg);
        } else {
            nes[i] = 0.0;
        }
    }
    
    Rcpp::CharacterVector leading_edge_str(n_sets);
    for (int i = 0; i < n_sets; ++i) {
        std::stringstream ss;
        ss << "tags=" << std::round(tags[i] * 100) << "%"
           << ", list=" << std::round(list[i] * 100) << "%"
           << ", signal=" << std::round(signal[i] * 100) << "%";
        leading_edge_str[i] = ss.str();
    }
    
    return Rcpp::DataFrame::create(
        Rcpp::Named("GeneSet") = gene_set_names,
        Rcpp::Named("ES") = es,
        Rcpp::Named("NES") = nes,
        Rcpp::Named("PValue") = pvalues,
        Rcpp::Named("Size") = Rcpp::wrap(gs_sizes),
        Rcpp::Named("rank") = rank,
        Rcpp::Named("leading_edge") = leading_edge_str,
        Rcpp::Named("core_enrichment") = core_enrichment
    );
}

} // namespace enrichit

// Adaptive GSEA with early-stopping and geometric batch scaling
namespace enrichit {

Rcpp::DataFrame gsea_adaptive(const Rcpp::NumericVector& stats,
                              const Rcpp::List& gene_sets,
                              const Rcpp::CharacterVector& gene_set_names,
                              int minPerm,
                              int maxPerm,
                              double pvalThreshold,
                              double exponent,
                              std::string method) {
    
    int n_sets = gene_sets.size();
    int n_genes = stats.size();
    
    std::vector<double> gene_stats = Rcpp::as<std::vector<double>>(stats);
    Rcpp::CharacterVector gene_names = stats.names();
    
    std::unordered_map<std::string, int> gene_map;
    for (int i = 0; i < n_genes; ++i) {
        gene_map[Rcpp::as<std::string>(gene_names[i])] = i;
    }
    
    std::vector<std::vector<bool>> gs_bools(n_sets, std::vector<bool>(n_genes, false));
    std::vector<int> gs_sizes(n_sets, 0);
    
    for (int i = 0; i < n_sets; ++i) {
        Rcpp::CharacterVector gs = gene_sets[i];
        int count = 0;
        for (int j = 0; j < gs.size(); ++j) {
            std::string g = Rcpp::as<std::string>(gs[j]);
            if (gene_map.find(g) != gene_map.end()) {
                gs_bools[i][gene_map[g]] = true;
                count++;
            }
        }
        gs_sizes[i] = count;
    }
    
    // Calculate observed ES and details for all gene sets
    Rcpp::NumericVector es(n_sets);
    Rcpp::IntegerVector rank(n_sets);
    Rcpp::NumericVector tags(n_sets);
    Rcpp::NumericVector list(n_sets);
    Rcpp::NumericVector signal(n_sets);
    Rcpp::CharacterVector core_enrichment(n_sets);
    
    for (int i = 0; i < n_sets; ++i) {
        if (gs_sizes[i] > 0) {
            GSEA_Result res = calculate_es_details(gene_stats, gs_bools[i], gene_names, exponent);
            es[i] = res.ES;
            rank[i] = res.rank;
            tags[i] = res.tags;
            list[i] = res.list;
            signal[i] = res.signal;
            core_enrichment[i] = res.core_enrichment;
        } else {
            es[i] = 0.0;
            rank[i] = 0;
            tags[i] = 0.0;
            list[i] = 0.0;
            signal[i] = 0.0;
            core_enrichment[i] = "";
        }
    }
    
    // Adaptive permutation results
    Rcpp::NumericVector pvalues(n_sets);
    Rcpp::NumericVector nes(n_sets);
    Rcpp::IntegerVector actual_perms(n_sets);
    
    // Process each gene set with adaptive permutation
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int i = 0; i < n_sets; ++i) {
        double obs_es = es[i];
        
        if (gs_sizes[i] == 0) {
            pvalues[i] = 1.0;
            nes[i] = 0.0;
            actual_perms[i] = 0;
            continue;
        }
        
        // Thread-local RNG for thread safety
        std::mt19937 rng(12345 + i * 1000);
        
        int total_perms = 0;
        int count_better = 0;
        double sum_pos_es = 0.0;
        double sum_neg_es = 0.0;
        int count_pos = 0;
        int count_neg = 0;
        
        int batch_size = minPerm;
        bool converged = false;
        
        // Prepare for sampling method
        std::vector<int> universe(n_genes);
        std::iota(universe.begin(), universe.end(), 0);
        int k = gs_sizes[i];
        std::vector<int> sample(k);
        
        // Prepare for permute method
        std::vector<int> perm_idx(n_genes);
        std::iota(perm_idx.begin(), perm_idx.end(), 0);
        
        while (!converged && total_perms < maxPerm) {
            // Run batch of permutations
            for (int p = 0; p < batch_size && (total_perms + p) < maxPerm; ++p) {
                double perm_es;
                
                if (method == "permute") {
                    std::shuffle(perm_idx.begin(), perm_idx.end(), rng);
                    perm_es = calculate_es_permute(gene_stats, gs_bools[i], perm_idx, exponent);
                } else {
                    // Sample method - faster
                    for (int j = 0; j < k; ++j) {
                        std::uniform_int_distribution<> dis(j, n_genes - 1);
                        int swap_idx = dis(rng);
                        std::swap(universe[j], universe[swap_idx]);
                        sample[j] = universe[j];
                    }
                    std::sort(sample.begin(), sample.end());
                    perm_es = calculate_es_sparse(gene_stats, sample, exponent);
                }
                
                // Update counts
                if (obs_es > 0) {
                    if (perm_es >= obs_es) count_better++;
                    if (perm_es >= 0) { sum_pos_es += perm_es; count_pos++; }
                } else {
                    if (perm_es <= obs_es) count_better++;
                    if (perm_es < 0) { sum_neg_es += perm_es; count_neg++; }
                }
            }
            
            total_perms += batch_size;
            
            // Calculate current p-value
            double current_pval = (double)(count_better + 1) / (double)(total_perms + 1);
            
            // Early stopping conditions
            if (total_perms >= minPerm) {
                if (current_pval > pvalThreshold) {
                    // Not significant, stop early
                    converged = true;
                } else if (total_perms >= maxPerm) {
                    // Reached max, stop
                    converged = true;
                } else {
                    // Significant, increase batch size geometrically
                    batch_size = std::min(batch_size * 2, maxPerm - total_perms);
                    if (batch_size <= 0) converged = true;
                }
            }
        }
        
        // Final p-value and NES calculation
        pvalues[i] = (double)(count_better + 1) / (double)(total_perms + 1);
        actual_perms[i] = total_perms;
        
        if (obs_es > 0) {
            double mean_pos = (count_pos > 0) ? (sum_pos_es / count_pos) : 1.0;
            nes[i] = obs_es / mean_pos;
        } else if (obs_es < 0) {
            double mean_neg = (count_neg > 0) ? (sum_neg_es / count_neg) : -1.0;
            nes[i] = obs_es / std::abs(mean_neg);
        } else {
            nes[i] = 0.0;
        }
    }
    
    // Create leading edge strings
    Rcpp::CharacterVector leading_edge_str(n_sets);
    for (int i = 0; i < n_sets; ++i) {
        std::stringstream ss;
        ss << "tags=" << std::round(tags[i] * 100) << "%"
           << ", list=" << std::round(list[i] * 100) << "%"
           << ", signal=" << std::round(signal[i] * 100) << "%";
        leading_edge_str[i] = ss.str();
    }
    
    return Rcpp::DataFrame::create(
        Rcpp::Named("GeneSet") = gene_set_names,
        Rcpp::Named("ES") = es,
        Rcpp::Named("NES") = nes,
        Rcpp::Named("PValue") = pvalues,
        Rcpp::Named("Size") = Rcpp::wrap(gs_sizes),
        Rcpp::Named("nPerm") = actual_perms,
        Rcpp::Named("rank") = rank,
        Rcpp::Named("leading_edge") = leading_edge_str,
        Rcpp::Named("core_enrichment") = core_enrichment
    );
}

} // namespace enrichit

// [[Rcpp::export]]
Rcpp::DataFrame gsea_cpp(const Rcpp::NumericVector& stats,
                         const Rcpp::List& gene_sets,
                         const Rcpp::CharacterVector& gene_set_names,
                         int nPerm = 1000,
                         double exponent = 1.0,
                         std::string method = "sample") {
    return enrichit::gsea(stats, gene_sets, gene_set_names, nPerm, exponent, method);
}

// [[Rcpp::export]]
Rcpp::DataFrame gsea_adaptive_cpp(const Rcpp::NumericVector& stats,
                                  const Rcpp::List& gene_sets,
                                  const Rcpp::CharacterVector& gene_set_names,
                                  int minPerm = 1000,
                                  int maxPerm = 100000,
                                  double pvalThreshold = 0.1,
                                  double exponent = 1.0,
                                  std::string method = "sample") {
    return enrichit::gsea_adaptive(stats, gene_sets, gene_set_names, minPerm, maxPerm, pvalThreshold, exponent, method);
}
