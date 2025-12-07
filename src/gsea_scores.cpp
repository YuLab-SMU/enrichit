#include <Rcpp.h>
#include <cmath>

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
