#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <unordered_map>

// [[Rcpp::export]]
Rcpp::DataFrame gsea_scores_cpp(const Rcpp::NumericVector& stats,
                                const Rcpp::LogicalVector& in_set,
                                double exponent) {
    int N = stats.size();
    Rcpp::NumericVector running_es(N);
    Rcpp::IntegerVector position(N);
    Rcpp::NumericVector x(N);

    double N_R = 0.0;
    int N_H = 0;

    for (int i = 0; i < N; ++i) {
        if (in_set[i]) {
            N_R += std::pow(std::abs(stats[i]), exponent);
            N_H++;
        }
        x[i] = i + 1;
    }

    double P_hit = 0.0;
    double P_miss = 0.0;
    double N_miss = (double)(N - N_H);

    for (int i = 0; i < N; ++i) {
        if (in_set[i]) {
            if (N_R != 0) {
                P_hit += std::pow(std::abs(stats[i]), exponent) / N_R;
            }
            position[i] = 1;
        } else {
            if (N_miss != 0) {
                P_miss += 1.0 / N_miss;
            }
            position[i] = 0;
        }
        running_es[i] = P_hit - P_miss;
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("x") = x,
        Rcpp::Named("runningScore") = running_es,
        Rcpp::Named("position") = position
    );
}

// Compute the classic weighted enrichment score (ES) for many gene sets at once.
// `stats` is a named numeric vector sorted in decreasing order; `gene_sets` is a
// list of character vectors. `scoreType` is one of "std", "pos", "neg":
//   - "std": ES is the running-sum deviation with maximum absolute value (sign kept);
//   - "pos": ES is the maximum (positive) deviation;
//   - "neg": ES is the minimum (negative) deviation.
// The ES definition matches the classic GSEA weighted running-sum statistic used
// by enrichit's gsea_cpp/calculate_es_details (hit weight |stat|^exponent / N_R,
// miss penalty 1/(N - N_H)). Computing ES from the hit positions only makes it
// cheap enough to evaluate on every permutation of the whole-pipeline NSEA null.
// [[Rcpp::export]]
Rcpp::NumericVector gsea_es_all_cpp(const Rcpp::NumericVector& stats,
                                    const Rcpp::List& gene_sets,
                                    double exponent,
                                    std::string scoreType) {
    int n_sets = gene_sets.size();
    int N = stats.size();
    Rcpp::NumericVector es_out(n_sets);

    Rcpp::CharacterVector gene_names = stats.names();
    std::unordered_map<std::string, int> gene_map;
    gene_map.reserve(N * 2);
    for (int i = 0; i < N; ++i) {
        gene_map[Rcpp::as<std::string>(gene_names[i])] = i;
    }

    std::vector<double> gene_stats = Rcpp::as<std::vector<double>>(stats);

    for (int s = 0; s < n_sets; ++s) {
        Rcpp::CharacterVector gs = gene_sets[s];
        std::vector<int> hits;
        hits.reserve(gs.size());
        for (int j = 0; j < gs.size(); ++j) {
            auto it = gene_map.find(Rcpp::as<std::string>(gs[j]));
            if (it != gene_map.end()) {
                hits.push_back(it->second);
            }
        }
        if (hits.empty()) {
            es_out[s] = 0.0;
            continue;
        }
        std::sort(hits.begin(), hits.end());

        int N_H = static_cast<int>(hits.size());
        double N_R = 0.0;
        for (int idx : hits) {
            N_R += std::pow(std::abs(gene_stats[idx]), exponent);
        }
        if (N_R == 0) {
            es_out[s] = 0.0;
            continue;
        }

        double N_miss = static_cast<double>(N - N_H);
        double dec_per_miss = (N_miss > 0) ? 1.0 / N_miss : 0.0;
        double P_hit = 0.0, P_miss = 0.0;
        double max_dev = 0.0, min_dev = 0.0;
        int last_idx = -1;

        auto update = [&](double dev) {
            if (dev > max_dev) max_dev = dev;
            if (dev < min_dev) min_dev = dev;
        };

        for (int idx : hits) {
            int misses = idx - last_idx - 1;
            if (misses > 0 && N_miss > 0) {
                P_miss += misses * dec_per_miss;
                update(P_hit - P_miss);
            }
            P_hit += std::pow(std::abs(gene_stats[idx]), exponent) / N_R;
            update(P_hit - P_miss);
            last_idx = idx;
        }
        int remaining = N - 1 - last_idx;
        if (remaining > 0 && N_miss > 0) {
            P_miss += remaining * dec_per_miss;
            update(P_hit - P_miss);
        }

        if (scoreType == "pos") {
            es_out[s] = max_dev;
        } else if (scoreType == "neg") {
            es_out[s] = min_dev;
        } else {
            es_out[s] = (std::abs(max_dev) >= std::abs(min_dev)) ? max_dev : min_dev;
        }
    }
    return es_out;
}
