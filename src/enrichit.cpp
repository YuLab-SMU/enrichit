#include "enrichit.h"
#include <sstream>

namespace enrichit {

// Helper function for hypergeometric distribution
double dhyper(int k, int m, int n, int k_plus_m_minus_n, bool log_p) {
    return R::dhyper((double)k, (double)m, (double)n, (double)k_plus_m_minus_n, (int)log_p);
}

Rcpp::DataFrame ora(const Rcpp::CharacterVector& gene,
                    const Rcpp::CharacterVector& universe,
                    const Rcpp::List& gene_sets,
                    const Rcpp::CharacterVector& gene_set_names) {
    
    int n_sets = gene_sets.size();
    Rcpp::IntegerVector overlap(n_sets);
    Rcpp::IntegerVector set_size(n_sets);
    Rcpp::IntegerVector de_size(n_sets);
    Rcpp::NumericVector p_value(n_sets);
    Rcpp::CharacterVector gene_id(n_sets);
    Rcpp::IntegerVector universe_size(n_sets);
    
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
        universe_size[i] = N;
        
        // Count overlap (DE genes in this gene set) and collect IDs
        int count = 0;
        std::stringstream ss;
        bool first = true;
        
        for (const auto& gene : de_set) {
            if (gs_bg.find(gene) != gs_bg.end()) {
                ++count;
                if (!first) ss << "/";
                ss << gene;
                first = false;
            }
        }
        overlap[i] = count;
        de_size[i] = K;
        gene_id[i] = ss.str();
        
        // Calculate Fisher's exact test p-value
        int m = M;      // total in gene set (in universe)
        int n = N - M;  // total NOT in gene set (in universe)
        int k = K;      // total DE genes (in universe)
        double p_val = R::phyper((double)(count - 1), (double)m, (double)n, (double)k, 0, 0);
        if (p_val < 0.0) p_val = 0.0;
        if (p_val > 1.0) p_val = 1.0;
        p_value[i] = p_val;
    }
    
    // Create result DataFrame
    Rcpp::DataFrame result = Rcpp::DataFrame::create(
        Rcpp::Named("GeneSet") = gene_set_names,
        Rcpp::Named("SetSize") = set_size,
        Rcpp::Named("DEInSet") = overlap,
        Rcpp::Named("DESize") = de_size,
        Rcpp::Named("UniverseSize") = universe_size,
        Rcpp::Named("PValue") = p_value,
        Rcpp::Named("geneID") = gene_id
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
