#include "enrichit.h"

namespace enrichit {

Rcpp::DataFrame ora(const Rcpp::CharacterVector& gene,
                    const Rcpp::CharacterVector& universe,
                    const Rcpp::List& gene_sets,
                    const Rcpp::CharacterVector& gene_set_names) {
    
    int n_sets = gene_sets.size();
    Rcpp::IntegerVector overlap(n_sets);
    Rcpp::IntegerVector set_size(n_sets);
    Rcpp::IntegerVector de_size(n_sets);
    Rcpp::NumericVector p_value(n_sets);
    
    // Convert gene to hash set for fast lookup
    Rcpp::CharacterVector de_genes = Rcpp::unique(gene);
    std::unordered_set<std::string> de_set;
    for (int i = 0; i < de_genes.size(); ++i) {
        de_set.insert(Rcpp::as<std::string>(de_genes[i]));
    }
    
    // Convert universe to hash set
    Rcpp::CharacterVector bg_genes = Rcpp::unique(universe);
    std::unordered_set<std::string> bg_set;
    for (int i = 0; i < bg_genes.size(); ++i) {
        bg_set.insert(Rcpp::as<std::string>(bg_genes[i]));
    }
    
    // Total genes in universe
    int N = bg_set.size();
    // DE genes in universe
    int K = 0;
    for (const auto& gene : de_set) {
        if (bg_set.find(gene) != bg_set.end()) {
            ++K;
        }
    }
    
    // Process each gene set
    for (int i = 0; i < n_sets; ++i) {
        Rcpp::CharacterVector gs = gene_sets[i];
        std::unordered_set<std::string> gs_set;
        for (int j = 0; j < gs.size(); ++j) {
            gs_set.insert(Rcpp::as<std::string>(gs[j]));
        }
        
        // Intersection with universe
        std::unordered_set<std::string> gs_bg;
        for (const auto& gene : gs_set) {
            if (bg_set.find(gene) != bg_set.end()) {
                gs_bg.insert(gene);
            }
        }
        
        int M = gs_bg.size(); // gene set size in universe
        set_size[i] = M;
        
        // Count overlap (DE genes in this gene set)
        int count = 0;
        for (const auto& gene : de_set) {
            if (gs_bg.find(gene) != gs_bg.end()) {
                ++count;
            }
        }
        overlap[i] = count;
        de_size[i] = K;
        
        // Calculate Fisher's exact test p-value
        // Hypergeometric parameters:
        // m = total successes in population (gene set size in universe)
        // n = total failures in population (genes NOT in gene set)
        // k = number of draws (total DE genes in universe)
        // x = observed successes (DE genes in gene set)
        int m = M;      // total in gene set (in universe)
        int n = N - M;  // total NOT in gene set (in universe)
        int k = K;      // total DE genes (in universe)
        int min_ak = std::min(m, k);
        
        double p_val = 0.0;
        for (int x = count; x <= min_ak; ++x) {
            p_val += dhyper(x, m, n, k);
        }
        p_value[i] = p_val;
    }
    
    // Create result DataFrame (only p-values, no correction or significance filtering)
    Rcpp::DataFrame result = Rcpp::DataFrame::create(
        Rcpp::Named("GeneSet") = gene_set_names,
        Rcpp::Named("SetSize") = set_size,
        Rcpp::Named("DEInSet") = overlap,
        Rcpp::Named("DESize") = de_size,
        Rcpp::Named("PValue") = p_value
    );
    
    return result;
}

} // namespace enrichit

// [[Rcpp::export]]
Rcpp::DataFrame ora_cpp(const Rcpp::CharacterVector& gene,
                        const Rcpp::CharacterVector& universe,
                        const Rcpp::List& gene_sets,
                        const Rcpp::CharacterVector& gene_set_names) {
    return enrichit::ora(gene, universe, gene_sets, gene_set_names);
}
