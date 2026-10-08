// Native set-test engine: camera / roast / fry, mirroring limma's
// algorithms (camera.default, .lmEffects, .roastEffects, .fryEffects,
// squeezeVar/fitFDist) so results agree with the limma backend.
// The Householder QR reproduces LINPACK dqrdc2's sign convention
// (beta = -dsign(rnormk, x_kk)) so that coefficient signs -- and hence
// Down/Up directions -- match limma exactly.
#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

// ------------------------------------------------------------------
// special functions (limma::logmdigamma, limma::trigammaInverse)
// ------------------------------------------------------------------

static double logmdigamma(double x) {
    // log(x) - digamma(x); series split at 5 as in limma
    if (x < 5.0) {
        double acc = 1.0 / x + 1.0 / (x + 1.0) + 1.0 / (x + 2.0) +
            1.0 / (x + 3.0) + 1.0 / (x + 4.0);
        return std::log(x / (x + 5.0)) + logmdigamma(x + 5.0) + acc;
    }
    double z2 = 1.0 / (x * x);
    double tail = z2 * (-1.0 / 12.0 + z2 * (1.0 / 120.0 +
        z2 * (-1.0 / 252.0 + z2 * (1.0 / 240.0 + z2 * (-1.0 / 132.0 +
        z2 * (691.0 / 32760.0 + z2 * (-1.0 / 12.0 + 3617.0 * z2 / 8160.0)))))));
    return 1.0 / (2.0 * x) - tail;
}

static double trigamma_inv(double x) {
    // limma::trigammaInverse
    if (x > 1e7) return 1.0 / std::sqrt(x);
    if (x < 1e-6) return 1.0 / x;
    double y = 0.5 + 1.0 / x;
    for (int it = 1; it <= 50; ++it) {
        double tri = Rf_trigamma(y);
        double dif = tri * (1.0 - tri / x) / Rf_psigamma(y, 2.0);
        y += dif;
        if (std::fabs(-dif / y) < 1e-8) break;
        if (it == 50) Rf_warning("Iteration limit exceeded");
    }
    return y;
}

static double sgn(double x) { return (x > 0.0) - (x < 0.0); }

// limma::zscoreT approximations
static double zscore_t_hill(double x, double df) {
    double A = df - 0.5;
    double B = 48.0 * A * A;
    double z = A * std::log1p(x / df * x);
    double w = (((((-0.4 * z - 3.3) * z - 24.0) * z - 85.5) /
        (0.8 * z * z + 100.0 + B) + z + 3.0) / B + 1.0) * std::sqrt(z);
    return w * sgn(x);
}

static double zscore_t_bailey(double x, double df) {
    double s = (df + 0.125) / (df + 1.125) *
        std::sqrt((df + 19.0 / 12.0) * std::log1p(x / (df + 1.0 / 12.0) * x));
    return s * sgn(x);
}

// ------------------------------------------------------------------
// eBayes prior fitting: limma::fitFDist (scalar df1, no covariate)
// ------------------------------------------------------------------

static void fit_fdist(std::vector<double> x, double df1,
                      double& scale, double& df2) {
    int n = (int)x.size();
    if (n == 0) { scale = NA_REAL; df2 = NA_REAL; return; }
    if (n == 1) { scale = x[0]; df2 = 0.0; return; }
    for (int i = 0; i < n; ++i) x[i] = std::max(x[i], 0.0);
    std::vector<double> xs(x);
    int mid = n / 2;
    std::nth_element(xs.begin(), xs.begin() + mid, xs.end());
    double m = xs[mid];
    if (n % 2 == 0) {
        std::nth_element(xs.begin(), xs.begin() + mid - 1, xs.end());
        m = 0.5 * (m + xs[mid - 1]);
    }
    if (m == 0.0) {
        Rf_warning("More than half of residual variances are exactly zero: eBayes unreliable");
        m = 1.0;
    } else {
        for (int i = 0; i < n; ++i)
            if (x[i] == 0.0) {
                Rf_warning("Zero sample variances detected, have been offset away from zero");
                break;
            }
    }
    for (int i = 0; i < n; ++i) x[i] = std::max(x[i], 1e-5 * m);
    double lgm = logmdigamma(df1 / 2.0);
    double emean = 0.0;
    for (int i = 0; i < n; ++i) emean += std::log(x[i]) + lgm;
    emean /= n;
    double evar = 0.0;
    for (int i = 0; i < n; ++i) {
        double e = std::log(x[i]) + lgm - emean;
        evar += e * e;
    }
    evar /= (n - 1);
    evar -= Rf_trigamma(df1 / 2.0);
    if (evar > 0.0) {
        df2 = 2.0 * trigamma_inv(evar);
        scale = std::exp(emean - logmdigamma(df2 / 2.0));
    } else {
        df2 = R_PosInf;
        double s = 0.0;
        for (int i = 0; i < n; ++i) s += x[i];
        scale = s / n;
    }
}

struct SVResult {
    double var_prior;
    double df_prior;
    std::vector<double> var_post;
};

// limma::squeezeVar, legacy path (scalar df, robust = FALSE, no covariate)
static SVResult squeeze_var(const std::vector<double>& var, double df) {
    SVResult out;
    int n = (int)var.size();
    if (n < 3) {
        out.var_post = var;
        // limma reports var.prior = var here; any finite value is
        // equivalent downstream (df.prior = 0 zeroes the term), but NA
        // would poison 0 * var_prior in the roast moderation
        out.var_prior = 0.0;
        out.df_prior = 0.0;
        return out;
    }
    fit_fdist(var, df, out.var_prior, out.df_prior);
    out.var_post.resize(n);
    if (!std::isfinite(out.df_prior) || out.df_prior > 1e100) {
        for (int i = 0; i < n; ++i) out.var_post[i] = out.var_prior;
    } else {
        for (int i = 0; i < n; ++i)
            out.var_post[i] = (df * var[i] + out.df_prior * out.var_prior) /
                (df + out.df_prior);
    }
    return out;
}

// ------------------------------------------------------------------
// Householder QR with LINPACK dqrdc2's sign convention
// ------------------------------------------------------------------

// A: n*p col-major, destroyed; Y: n*m col-major, becomes Q^T Y.
// Returns R[p-1, p-1] (the beta at the final reflector).
// rank_tol: also reject a column whose remaining norm falls below 1e-7
// of its original norm, LINPACK dqrdc2's collapse criterion; limma stops
// whenever any column collapses, so the verdict matches without
// pivoting. Weighted per-gene calls pass false, mirroring .lmEffects,
// which only rank-checks the unweighted design.
static double qty_householder(int n, int p, double* A, int m, double* Y,
                              bool rank_tol = false) {
    std::vector<double> u(n);
    std::vector<double> onorm;
    if (rank_tol) {
        onorm.resize(p);
        for (int j = 0; j < p; ++j) {
            double s = 0.0;
            for (int i = 0; i < n; ++i)
                s += A[i + (size_t)j * n] * A[i + (size_t)j * n];
            onorm[j] = std::sqrt(s);
        }
    }
    double rpp = 0.0;
    for (int k = 0; k < p; ++k) {
        double rnorm2 = 0.0;
        for (int i = k; i < n; ++i) rnorm2 += A[i + k * n] * A[i + k * n];
        double rnorm = std::sqrt(rnorm2);
        if (rnorm == 0.0 || (rank_tol && rnorm <= 1e-7 * onorm[k]))
            Rcpp::stop("design matrix is not of full rank");
        double xkk = A[k + k * n];
        double beta = (xkk >= 0.0) ? -rnorm : rnorm;
        for (int i = k; i < n; ++i) u[i - k] = A[i + k * n];
        u[0] = xkk - beta;
        double uu = 0.0;
        for (int i = k; i < n; ++i) uu += u[i - k] * u[i - k];
        for (int j = k + 1; j < p; ++j) {
            double t = 0.0;
            for (int i = k; i < n; ++i) t += u[i - k] * A[i + j * n];
            t *= 2.0 / uu;
            for (int i = k; i < n; ++i) A[i + j * n] -= t * u[i - k];
        }
        for (int i = k + 1; i < n; ++i) A[i + k * n] = 0.0;
        A[k + k * n] = beta;
        for (int j = 0; j < m; ++j) {
            double* Yj = Y + (size_t)j * n;
            double t = 0.0;
            for (int i = k; i < n; ++i) t += u[i - k] * Yj[i];
            t *= 2.0 / uu;
            for (int i = k; i < n; ++i) Yj[i] -= t * u[i - k];
        }
        rpp = beta;
    }
    return rpp;
}

// ------------------------------------------------------------------
// .lmEffects equivalent: effects matrix G x (n - p + 1); column 1 is
// the contrast coordinate (signed as limma does), columns 2.. the
// n - p residual coordinates.
// ------------------------------------------------------------------

static void check_weights(const Rcpp::NumericMatrix& W, int G, int n) {
    if (W.nrow() != G || W.ncol() != n)
        Rcpp::stop("weights must have the same dimensions as E");
    for (int i = 0; i < W.size(); ++i)
        if (!(W[i] > 0.0))
            Rcpp::stop("weights must be positive");
}

static void lm_effects(const Rcpp::NumericMatrix& E,
                       const Rcpp::NumericMatrix& W, bool wtd,
                       const Rcpp::NumericMatrix& design, int contrast_idx,
                       Rcpp::NumericMatrix& Eff) {
    int G = E.rows(), n = E.cols(), p = design.cols();
    if (n <= p) Rcpp::stop("No residual degrees of freedom");
    for (int i = 0; i < E.size(); ++i)
        if (!R_finite(E[i]))
            Rcpp::stop("All y values must be finite and non-NA");
    // contrast column last, preserving the order of the others
    Rcpp::NumericMatrix X(n, p);
    int c = 0;
    for (int j = 0; j < p; ++j)
        if (j != contrast_idx - 1) X(Rcpp::_, c++) = design(Rcpp::_, j);
    X(Rcpp::_, p - 1) = design(Rcpp::_, contrast_idx - 1);

    Eff = Rcpp::NumericMatrix(G, n - p + 1);
    std::vector<double> A((size_t)n * p);
    if (!wtd) {
        for (int j = 0; j < p; ++j)
            for (int i = 0; i < n; ++i) A[i + (size_t)j * n] = X(i, j);
        std::vector<double> Y((size_t)n * G);
        for (int g = 0; g < G; ++g)
            for (int i = 0; i < n; ++i) Y[i + (size_t)g * n] = E(g, i);
        double rpp = qty_householder(n, p, A.data(), G, Y.data(), true);
        double s = (rpp < 0.0) ? -1.0 : 1.0;
        for (int g = 0; g < G; ++g) {
            Eff(g, 0) = s * Y[(p - 1) + (size_t)g * n];
            for (int k = 1; k <= n - p; ++k)
                Eff(g, k) = Y[(p - 1 + k) + (size_t)g * n];
        }
    } else {
        std::vector<double> y(n);
        for (int g = 0; g < G; ++g) {
            for (int j = 0; j < p; ++j)
                for (int i = 0; i < n; ++i)
                    A[i + (size_t)j * n] = X(i, j) * std::sqrt(W(g, i));
            for (int i = 0; i < n; ++i) y[i] = E(g, i) * std::sqrt(W(g, i));
            double rpp = qty_householder(n, p, A.data(), 1, y.data());
            double s = (rpp < 0.0) ? -1.0 : 1.0;
            Eff(g, 0) = s * y[p - 1];
            for (int k = 1; k <= n - p; ++k) Eff(g, k) = y[p - 1 + k];
        }
    }
}

static std::vector<std::vector<int>> as_index(const Rcpp::List& index) {
    std::vector<std::vector<int>> out(index.size());
    for (int i = 0; i < index.size(); ++i)
        out[i] = Rcpp::as<std::vector<int>>(index[i]);
    return out;
}

// ------------------------------------------------------------------
// Hestenes one-sided Jacobi SVD, singular values only (descending)
// ------------------------------------------------------------------

static std::vector<double> svd_values(std::vector<double> A, int m, int n) {
    std::vector<double> B(A);  // working copy, col-major m x n
    const int max_sweep = 60;
    for (int sweep = 0; sweep < max_sweep; ++sweep) {
        bool rotated = false;
        for (int pp = 0; pp < n - 1; ++pp) {
            for (int q = pp + 1; q < n; ++q) {
                double app = 0.0, aqq = 0.0, apq = 0.0;
                for (int i = 0; i < m; ++i) {
                    double xp = B[i + (size_t)pp * m], xq = B[i + (size_t)q * m];
                    app += xp * xp; aqq += xq * xq; apq += xp * xq;
                }
                if (std::fabs(apq) <= 1e-15 * std::sqrt(app * aqq)) continue;
                rotated = true;
                double tau = (aqq - app) / (2.0 * apq);
                double t = (tau >= 0.0 ? 1.0 : -1.0) /
                    (std::fabs(tau) + std::sqrt(1.0 + tau * tau));
                double cs = 1.0 / std::sqrt(1.0 + t * t);
                double sn = cs * t;
                for (int i = 0; i < m; ++i) {
                    double xp = B[i + (size_t)pp * m], xq = B[i + (size_t)q * m];
                    B[i + (size_t)pp * m] = cs * xp - sn * xq;
                    B[i + (size_t)q * m] = sn * xp + cs * xq;
                }
            }
        }
        if (!rotated) break;
    }
    std::vector<double> vals(n);
    for (int j = 0; j < n; ++j) {
        double s = 0.0;
        for (int i = 0; i < m; ++i) s += B[i + (size_t)j * m] * B[i + (size_t)j * m];
        vals[j] = std::sqrt(s);
    }
    std::sort(vals.begin(), vals.end(), std::greater<double>());
    return vals;
}

// Gauss-Legendre nodes/weights mapped to [0, 1]
static void gauss_legendre01(int n, std::vector<double>& nodes,
                             std::vector<double>& wts) {
    nodes.assign(n, 0.0);
    wts.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double z = std::cos(M_PI * (i + 0.75) / (n + 0.5));
        double pp = 1.0;
        for (int it = 0; it < 100; ++it) {
            double p0 = 1.0, p1 = z;
            for (int j = 2; j <= n; ++j) {
                double p2 = ((2.0 * j - 1.0) * z * p1 - (j - 1.0) * p0) / j;
                p0 = p1; p1 = p2;
            }
            pp = n * (z * p1 - p0) / (z * z - 1.0);
            double dz = p1 / pp;
            z -= dz;
            if (std::fabs(dz) < 1e-15) break;
        }
        nodes[i] = 0.5 * (z + 1.0);
        wts[i] = 1.0 / ((1.0 - z * z) * pp * pp);
    }
}

//' Effects matrix from .lmEffects (contrast coordinate + residual coordinates)
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix st_effects_cpp(Rcpp::NumericMatrix E,
                                   Rcpp::Nullable<Rcpp::NumericMatrix> weights,
                                   Rcpp::NumericMatrix design,
                                   int contrast_idx) {
    int G = E.rows();
    Rcpp::NumericMatrix W, Eff;
    bool wtd = false;
    if (weights.isNotNull()) {
        W = Rcpp::as<Rcpp::NumericMatrix>(weights);
        check_weights(W, G, E.cols());
        wtd = true;
    }
    lm_effects(E, W, wtd, design, contrast_idx, Eff);
    return Eff;
}

//' limma-style squeezeVar for scalar residual df (fitFDist prior)
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List st_squeeze_var_cpp(std::vector<double> var, double df) {
    SVResult sv = squeeze_var(var, df);
    return Rcpp::List::create(
        Rcpp::Named("var.prior") = sv.var_prior,
        Rcpp::Named("df.prior") = sv.df_prior,
        Rcpp::Named("var.post") = sv.var_post);
}

//' Native camera: inter-gene-correlation adjusted mean-rank set test
//'
//' Mirrors \code{limma::camera} for \code{use.ranks = FALSE}.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List st_camera_cpp(Rcpp::NumericMatrix E,
                         Rcpp::Nullable<Rcpp::NumericMatrix> weights,
                         Rcpp::NumericMatrix design, int contrast_idx,
                         Rcpp::List index,
                         Rcpp::Nullable<double> inter_gene_cor,
                         bool allow_neg_cor) {
    int G = E.rows(), n = E.cols();
    if (G < 3) Rcpp::stop("Too few genes in dataset: need at least 3");
    int p = design.cols();
    if (design.rows() != n)
        Rcpp::stop("row dimension of design matrix must match column dimension of data");
    int dfr = n - p;
    if (dfr < 1) Rcpp::stop("No residual df: cannot compute t-tests");
    for (int i = 0; i < E.size(); ++i)
        if (!R_finite(E[i]))
            Rcpp::stop("All y values must be finite and non-NA");
    bool fixed_cor = inter_gene_cor.isNotNull();
    double cor_fixed = fixed_cor ? Rcpp::as<double>(inter_gene_cor) : 0.0;
    int df_camera = fixed_cor ? (G - 2) : std::min(dfr, G - 2);

    Rcpp::NumericMatrix W;
    bool wtd = false;
    if (weights.isNotNull()) {
        W = Rcpp::as<Rcpp::NumericMatrix>(weights);
        check_weights(W, G, n);
        wtd = true;
    }

    // contrast column last, preserving the order of the others
    Rcpp::NumericMatrix X(n, p);
    int c = 0;
    for (int j = 0; j < p; ++j)
        if (j != contrast_idx - 1) X(Rcpp::_, c++) = design(Rcpp::_, j);
    X(Rcpp::_, p - 1) = design(Rcpp::_, contrast_idx - 1);

    std::vector<double> unscaledt(G), sigma2(G);
    std::vector<double> U((size_t)dfr * G);
    std::vector<double> A((size_t)n * p), Y((size_t)n * G), y(n);
    if (!wtd) {
        for (int j = 0; j < p; ++j)
            for (int i = 0; i < n; ++i) A[i + (size_t)j * n] = X(i, j);
        for (int g = 0; g < G; ++g)
            for (int i = 0; i < n; ++i) Y[i + (size_t)g * n] = E(g, i);
        double rpp = qty_householder(n, p, A.data(), G, Y.data(), true);
        double s = (rpp < 0.0) ? -1.0 : 1.0;
        for (int g = 0; g < G; ++g) {
            unscaledt[g] = s * Y[(p - 1) + (size_t)g * n];
            double s2 = 0.0;
            for (int k = 0; k < dfr; ++k) {
                double v = Y[(p + k) + (size_t)g * n];
                U[k * (size_t)G + g] = v;
                s2 += v * v;
            }
            sigma2[g] = s2 / dfr;
        }
    } else {
        for (int g = 0; g < G; ++g) {
            for (int j = 0; j < p; ++j)
                for (int i = 0; i < n; ++i)
                    A[i + (size_t)j * n] = X(i, j) * std::sqrt(W(g, i));
            for (int i = 0; i < n; ++i) y[i] = E(g, i) * std::sqrt(W(g, i));
            double rpp = qty_householder(n, p, A.data(), 1, y.data());
            double s = (rpp < 0.0) ? -1.0 : 1.0;
            unscaledt[g] = s * y[p - 1];
            double s2 = 0.0;
            for (int k = 0; k < dfr; ++k) {
                double v = y[p + k];
                U[k * (size_t)G + g] = v;
                s2 += v * v;
            }
            sigma2[g] = s2 / dfr;
        }
    }

    std::vector<double> Ustd;
    if (!fixed_cor) {
        Ustd.assign(U.size(), 0.0);
        for (int g = 0; g < G; ++g) {
            double sd = std::sqrt(std::max(sigma2[g], 1e-8));
            for (int k = 0; k < dfr; ++k)
                Ustd[k * (size_t)G + g] = U[k * (size_t)G + g] / sd;
        }
    }

    SVResult sv = squeeze_var(sigma2, dfr);
    std::vector<double> modt(G), Stat(G);
    for (int g = 0; g < G; ++g)
        modt[g] = unscaledt[g] / std::sqrt(sv.var_post[g]);
    double df_total = std::min(dfr + sv.df_prior, (double)G * dfr);
    for (int g = 0; g < G; ++g) Stat[g] = zscore_t_hill(modt[g], df_total);
    double meanStat = 0.0;
    for (int g = 0; g < G; ++g) meanStat += Stat[g];
    meanStat /= G;
    double varStat = 0.0;
    for (int g = 0; g < G; ++g) varStat += (Stat[g] - meanStat) * (Stat[g] - meanStat);
    varStat /= (G - 1);

    std::vector<std::vector<int>> sets = as_index(index);
    int nsets = (int)sets.size();
    Rcpp::IntegerVector ngenes(nsets);
    Rcpp::NumericVector correlation(nsets), down(nsets), up(nsets),
        twosided(nsets);
    for (int i = 0; i < nsets; ++i) {
        const std::vector<int>& iset = sets[i];
        int m = (int)iset.size();
        int m2 = G - m;
        ngenes[i] = m;
        double corr_i, vif;
        if (fixed_cor) {
            corr_i = cor_fixed;
            vif = 1.0 + (m - 1) * cor_fixed;
        } else if (m > 1) {
            // vif = m * mean over residual coords of (mean over set genes)^2
            double acc = 0.0;
            for (int k = 0; k < dfr; ++k) {
                double cs = 0.0;
                for (int a = 0; a < m; ++a) cs += Ustd[k * (size_t)G + iset[a]];
                cs /= m;
                acc += cs * cs;
            }
            acc /= dfr;
            vif = m * acc;
            corr_i = (vif - 1.0) / (m - 1);
        } else {
            vif = 1.0;
            corr_i = NA_REAL;
        }
        if (!allow_neg_cor) vif = std::max(1.0, vif);
        double mis = 0.0;
        for (int a = 0; a < m; ++a) mis += Stat[iset[a]];
        mis /= m;
        double delta = (double)G / m2 * (mis - meanStat);
        double varPooled = (((double)G - 1.0) * varStat -
            delta * delta * m * (double)m2 / G) / (G - 2.0);
        double t2 = delta /
            std::sqrt(varPooled * (vif / m + 1.0 / m2));
        down[i] = Rf_pt(t2, df_camera, 1, 0);
        up[i] = Rf_pt(t2, df_camera, 0, 0);
        twosided[i] = 2.0 * std::min(down[i], up[i]);
        correlation[i] = corr_i;
    }
    return Rcpp::List::create(
        Rcpp::Named("ngenes") = ngenes,
        Rcpp::Named("correlation") = correlation,
        Rcpp::Named("down") = down,
        Rcpp::Named("up") = up,
        Rcpp::Named("twosided") = twosided,
        Rcpp::Named("fixed_cor") = fixed_cor);
}

//' Native roast: rotation gene set test on the effects matrix
//'
//' Mirrors \code{limma::mroast}/\code{.roastEffects} with scalar
//' var.prior/df.prior and \code{approx.zscore = TRUE}. Consumes R's
//' RNG in the same order as limma, so \code{set.seed()} reproduces
//' limma's rotation draws exactly.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List st_roast_cpp(Rcpp::NumericMatrix E,
                        Rcpp::Nullable<Rcpp::NumericMatrix> weights,
                        Rcpp::NumericMatrix design, int contrast_idx,
                        Rcpp::List index, std::string set_statistic,
                        Rcpp::Nullable<Rcpp::NumericVector> gene_weights,
                        int nrot) {
    if (set_statistic != "mean" && set_statistic != "floormean" &&
        set_statistic != "mean50" && set_statistic != "msq")
        Rcpp::stop("set.statistic must be one of 'mean', 'floormean', 'mean50', 'msq'");
    if (nrot < 1) Rcpp::stop("nrot must be positive");
    int G = E.rows();
    Rcpp::NumericMatrix W, Eff;
    bool wtd = false;
    if (weights.isNotNull()) {
        W = Rcpp::as<Rcpp::NumericMatrix>(weights);
        check_weights(W, G, E.cols());
        wtd = true;
    }
    lm_effects(E, W, wtd, design, contrast_idx, Eff);
    int ne = Eff.cols();
    int dfr = ne - 1;

    std::vector<double> s2(G);
    for (int g = 0; g < G; ++g) {
        double acc = 0.0;
        for (int k = 1; k < ne; ++k) acc += Eff(g, k) * Eff(g, k);
        s2[g] = acc / dfr;
    }
    SVResult sv = squeeze_var(s2, dfr);
    double var_prior = sv.var_prior, df_prior = sv.df_prior;
    double df_total = df_prior + dfr;
    double df_winsor = std::min(df_total, 10000.0);

    std::vector<double> gw;
    if (gene_weights.isNotNull()) {
        gw = Rcpp::as<std::vector<double>>(gene_weights);
        if ((int)gw.size() != G)
            Rcpp::stop("gene.weights vector should be of length nrow(E)");
    }

    std::vector<std::vector<int>> sets = as_index(index);
    int nsets = (int)sets.size();
    Rcpp::NumericMatrix pv(nsets, 4), active(nsets, 4);
    Rcpp::IntegerVector ngenes(nsets);

    const int chunk = 1000;
    int nchunk = (nrot + chunk - 1) / chunk;
    int nroti0 = (nrot + nchunk - 1) / nchunk;
    int overshoot = nchunk * nroti0 - nrot;
    const double chimed = Rf_qnorm5(0.25, 0.0, 1.0, 0, 0);
    const double sq2 = std::sqrt(2.0);

    for (int i = 0; i < nsets; ++i) {
        const std::vector<int>& iset = sets[i];
        int m = (int)iset.size();
        ngenes[i] = m;
        std::vector<double> ES((size_t)m * ne), re2(m), vp(m);
        for (int a = 0; a < m; ++a) {
            int g = iset[a];
            for (int k = 0; k < ne; ++k) ES[(size_t)a * ne + k] = Eff(g, k);
            double acc = 0.0;
            for (int k = 0; k < ne; ++k)
                acc += ES[(size_t)a * ne + k] * ES[(size_t)a * ne + k];
            re2[a] = acc;
            vp[a] = sv.var_post[g];
        }
        std::vector<double> gws(m, 1.0);
        bool has_gw = !gw.empty();
        if (has_gw)
            for (int a = 0; a < m; ++a) gws[a] = gw[iset[a]];
        // limma's msq chunk branch rescales gene.weights in place
        std::vector<double> gw_acc(m, 0.0);
        if (has_gw)
            for (int a = 0; a < m; ++a) gw_acc[a] = std::fabs(gws[a]);

        std::vector<double> modt(m);
        for (int a = 0; a < m; ++a)
            modt[a] = zscore_t_bailey(ES[(size_t)a * ne] / std::sqrt(vp[a]),
                                      df_winsor);
        double a1, a2;
        if (!has_gw) {
            a1 = 0.0; a2 = 0.0;
            for (int a = 0; a < m; ++a) {
                if (modt[a] > sq2) a1 += 1.0;
                if (modt[a] < -sq2) a2 += 1.0;
            }
            a1 /= m; a2 /= m;
        } else {
            double ss = 0.0; a1 = 0.0; a2 = 0.0;
            for (int a = 0; a < m; ++a) ss += std::fabs(sgn(gws[a]));
            for (int a = 0; a < m; ++a) {
                double sm = sgn(gws[a]) * modt[a];
                if (sm > sq2) a1 += 1.0;
                if (sm < -sq2) a2 += 1.0;
            }
            a1 /= ss; a2 /= ss;  // 0/0 = NaN with all-zero weights, as in limma
        }

        // observed statistics: down, up, upordown, mixed
        double obs[4] = {0.0, 0.0, 0.0, 0.0};
        if (set_statistic == "mean") {
            if (has_gw)
                for (int a = 0; a < m; ++a) modt[a] *= gws[a];
            double mm = 0.0, am = 0.0;
            for (int a = 0; a < m; ++a) { mm += modt[a]; am += std::fabs(modt[a]); }
            mm /= m; am /= m;
            obs[0] = -mm; obs[1] = mm; obs[3] = am;
        } else if (set_statistic == "floormean") {
            std::vector<double> amodt(m);
            for (int a = 0; a < m; ++a)
                amodt[a] = std::max(std::fabs(modt[a]), chimed);
            if (has_gw)
                for (int a = 0; a < m; ++a) {
                    amodt[a] *= gws[a];
                    modt[a] *= gws[a];
                }
            for (int a = 0; a < m; ++a) {
                obs[0] += std::max(-modt[a], 0.0);
                obs[1] += std::max(modt[a], 0.0);
                obs[3] += amodt[a];
            }
            obs[0] /= m; obs[1] /= m; obs[3] /= m;
            obs[2] = std::max(obs[0], obs[1]);
        } else if (set_statistic == "mean50") {
            int half1, half2;
            if (m % 2 == 0) { half1 = m / 2; half2 = half1 + 1; }
            else half1 = half2 = m / 2 + 1;
            if (has_gw)
                for (int a = 0; a < m; ++a) modt[a] *= gws[a];
            std::vector<double> s1(modt);
            std::nth_element(s1.begin(), s1.begin() + half2 - 1, s1.end());
            double dacc = 0.0, uacc = 0.0;
            for (int a = 0; a < half1; ++a) dacc += s1[a];
            for (int a = half2 - 1; a < m; ++a) uacc += s1[a];
            obs[0] = -dacc / half1;
            obs[1] = uacc / (m - half2 + 1);
            obs[2] = std::max(obs[0], obs[1]);
            std::vector<double> s2v(m);
            for (int a = 0; a < m; ++a) s2v[a] = std::fabs(modt[a]);
            std::nth_element(s2v.begin(), s2v.begin() + half2 - 1, s2v.end());
            double macc = 0.0;
            for (int a = half2 - 1; a < m; ++a) macc += s2v[a];
            obs[3] = macc / (m - half2 + 1);
        } else {  // msq
            std::vector<double> modt2(m);
            for (int a = 0; a < m; ++a) modt2[a] = modt[a] * modt[a];
            if (has_gw) {
                for (int a = 0; a < m; ++a) {
                    modt2[a] *= std::fabs(gws[a]);
                    modt[a] *= gws[a];
                }
            }
            for (int a = 0; a < m; ++a) {
                if (modt[a] < 0.0) obs[0] += modt2[a];
                if (modt[a] > 0.0) obs[1] += modt2[a];
                obs[3] += modt2[a];
            }
            obs[0] /= m; obs[1] /= m; obs[3] /= m;
            obs[2] = std::max(obs[0], obs[1]);
        }

        long count[4] = {0, 0, 0, 0};
        int nroti = nroti0;
        std::vector<double> M, rot((size_t)m * nroti0), amodtr;
        for (int ch = 0; ch < nchunk; ++ch) {
            if (ch == nchunk - 1) nroti = nroti0 - overshoot;
            M.assign((size_t)nroti * ne, 0.0);
            for (int j = 0; j < ne; ++j)
                for (int r = 0; r < nroti; ++r)
                    M[r + (size_t)j * nroti] = R::norm_rand();
            for (int r = 0; r < nroti; ++r) {
                double n2 = 0.0;
                for (int j = 0; j < ne; ++j)
                    n2 += M[r + (size_t)j * nroti] * M[r + (size_t)j * nroti];
                n2 = std::sqrt(n2);
                for (int j = 0; j < ne; ++j) M[r + (size_t)j * nroti] /= n2;
            }
            // rotated effects: rot[a, r] = sum_k ES[a, k] * M[r, k]
            std::fill(rot.begin(), rot.end(), 0.0);
            for (int a = 0; a < m; ++a) {
                for (int r = 0; r < nroti; ++r) {
                    double acc = 0.0;
                    for (int k = 0; k < ne; ++k)
                        acc += ES[(size_t)a * ne + k] * M[r + (size_t)k * nroti];
                    rot[(size_t)a * nroti + r] = acc;
                }
            }
            for (int a = 0; a < m; ++a) {
                for (int r = 0; r < nroti; ++r) {
                    double sr = (re2[a] -
                        rot[(size_t)a * nroti + r] * rot[(size_t)a * nroti + r]) / dfr;
                    if (std::isfinite(df_prior))
                        sr = (df_prior * var_prior + dfr * sr) / df_total;
                    else
                        sr = var_prior;
                    rot[(size_t)a * nroti + r] = zscore_t_bailey(
                        rot[(size_t)a * nroti + r] / std::sqrt(sr), df_winsor);
                }
            }

            if (set_statistic == "mean") {
                if (has_gw)
                    for (int a = 0; a < m; ++a)
                        for (int r = 0; r < nroti; ++r)
                            rot[(size_t)a * nroti + r] *= gws[a];
                for (int r = 0; r < nroti; ++r) {
                    double mm = 0.0, am = 0.0;
                    for (int a = 0; a < m; ++a) {
                        double v = rot[(size_t)a * nroti + r];
                        mm += v; am += std::fabs(v);
                    }
                    mm /= m; am /= m;
                    count[0] += (-mm > obs[0]) + (mm > obs[0]);
                    count[1] += (-mm > obs[1]) + (mm > obs[1]);
                    count[3] += am > obs[3];
                }
            } else if (set_statistic == "floormean") {
                // limma floors |raw| at chimed first, then multiplies by gw
                amodtr.assign((size_t)m * nroti, 0.0);
                for (size_t k = 0; k < amodtr.size(); ++k)
                    amodtr[k] = std::max(std::fabs(rot[k]), chimed);
                if (has_gw)
                    for (int a = 0; a < m; ++a)
                        for (int r = 0; r < nroti; ++r) {
                            rot[(size_t)a * nroti + r] *= gws[a];
                            amodtr[(size_t)a * nroti + r] *= gws[a];
                        }
                for (int r = 0; r < nroti; ++r) {
                    double d = 0.0, u = 0.0, mx = 0.0;
                    for (int a = 0; a < m; ++a) {
                        double v = rot[(size_t)a * nroti + r];
                        d += std::max(-v, 0.0);
                        u += std::max(v, 0.0);
                        mx += amodtr[(size_t)a * nroti + r];
                    }
                    d /= m; u /= m; mx /= m;
                    count[0] += (d > obs[0]) + (u > obs[0]);
                    count[1] += (d > obs[1]) + (u > obs[1]);
                    count[2] += std::max(d, u) > obs[2];
                    count[3] += mx > obs[3];
                }
            } else if (set_statistic == "mean50") {
                if (has_gw)
                    for (int a = 0; a < m; ++a)
                        for (int r = 0; r < nroti; ++r)
                            rot[(size_t)a * nroti + r] *= gws[a];
                int half1, half2;
                if (m % 2 == 0) { half1 = m / 2; half2 = half1 + 1; }
                else half1 = half2 = m / 2 + 1;
                std::vector<double> col(m);
                for (int r = 0; r < nroti; ++r) {
                    for (int a = 0; a < m; ++a)
                        col[a] = rot[(size_t)a * nroti + r];
                    std::nth_element(col.begin(), col.begin() + half2 - 1, col.end());
                    double dacc = 0.0, uacc = 0.0;
                    for (int a = 0; a < half1; ++a) dacc += col[a];
                    for (int a = half2 - 1; a < m; ++a) uacc += col[a];
                    double d = -dacc / half1, u = uacc / (m - half2 + 1);
                    for (int a = 0; a < m; ++a)
                        col[a] = std::fabs(rot[(size_t)a * nroti + r]);
                    std::nth_element(col.begin(), col.begin() + half2 - 1, col.end());
                    double macc = 0.0;
                    for (int a = half2 - 1; a < m; ++a) macc += col[a];
                    count[0] += (d > obs[0]) + (u > obs[0]);
                    count[1] += (d > obs[1]) + (u > obs[1]);
                    count[2] += std::max(d, u) > obs[2];
                    count[3] += macc / (m - half2 + 1) > obs[3];
                }
            } else {  // msq
                if (has_gw)
                    // limma rewrites gene.weights <- sqrt(abs(gw)) inside
                    // the chunk loop, so each chunk applies a smaller root
                    for (int a = 0; a < m; ++a) {
                        gw_acc[a] = std::sqrt(gw_acc[a]);
                        for (int r = 0; r < nroti; ++r)
                            rot[(size_t)a * nroti + r] *= gw_acc[a];
                    }
                for (int r = 0; r < nroti; ++r) {
                    double d = 0.0, u = 0.0, mx = 0.0;
                    for (int a = 0; a < m; ++a) {
                        double v = rot[(size_t)a * nroti + r];
                        d += std::max(-v, 0.0) * std::max(-v, 0.0);
                        u += std::max(v, 0.0) * std::max(v, 0.0);
                        mx += v * v;
                    }
                    d /= m; u /= m; mx /= m;
                    count[0] += (d > obs[0]) + (u > obs[0]);
                    count[1] += (d > obs[1]) + (u > obs[1]);
                    count[2] += std::max(d, u) > obs[2];
                    count[3] += mx > obs[3];
                }
            }
        }
        if (set_statistic == "mean")
            count[2] = std::min(count[0], count[1]);
        pv(i, 0) = (count[0] + 1.0) / (2.0 * nrot + 1.0);
        pv(i, 1) = (count[1] + 1.0) / (2.0 * nrot + 1.0);
        pv(i, 2) = (count[2] + 1.0) / (1.0 * nrot + 1.0);
        pv(i, 3) = (count[3] + 1.0) / (1.0 * nrot + 1.0);
        active(i, 0) = a2;
        active(i, 1) = a1;
        active(i, 2) = std::max(a1, a2);
        active(i, 3) = a1 + a2;
    }
    return Rcpp::List::create(
        Rcpp::Named("pv") = pv,
        Rcpp::Named("active") = active,
        Rcpp::Named("ngenes") = ngenes);
}

//' Native fry: fast rotation-free set test on the effects matrix
//'
//' Mirrors \code{limma::fry} (standardize = "posterior.sd" by default).
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List st_fry_cpp(Rcpp::NumericMatrix E,
                      Rcpp::Nullable<Rcpp::NumericMatrix> weights,
                      Rcpp::NumericMatrix design, int contrast_idx,
                      Rcpp::List index,
                      Rcpp::Nullable<Rcpp::NumericVector> gene_weights,
                      std::string standardize) {
    if (standardize != "none" && standardize != "residual.sd" &&
        standardize != "posterior.sd" && standardize != "p2")
        Rcpp::stop("standardize must be one of 'none', 'residual.sd', 'posterior.sd', 'p2'");
    int G = E.rows();
    Rcpp::NumericMatrix W, Eff;
    bool wtd = false;
    if (weights.isNotNull()) {
        W = Rcpp::as<Rcpp::NumericMatrix>(weights);
        check_weights(W, G, E.cols());
        wtd = true;
    }
    lm_effects(E, W, wtd, design, contrast_idx, Eff);
    int ne = Eff.cols();
    int dfr = ne - 1;

    if (standardize != "none") {
        std::vector<double> nodes, wts;
        gauss_legendre01(128, nodes, wts);
        double Eu2max = 0.0;
        for (int i = 0; i < 128; ++i)
            Eu2max += (dfr + 1.0) * std::pow(nodes[i], (double)dfr) *
                Rf_qchisq(nodes[i], 1.0, 1, 0) * wts[i];
        std::vector<double> s2r(G), s2(G);
        for (int g = 0; g < G; ++g) {
            double sum2 = 0.0, mx = 0.0;
            for (int k = 0; k < ne; ++k) {
                double v = Eff(g, k) * Eff(g, k);
                sum2 += v;
                if (v > mx) mx = v;
            }
            s2r[g] = (sum2 - mx) / (dfr + 1.0 - Eu2max);
            double acc = 0.0;
            for (int k = 1; k < ne; ++k) acc += Eff(g, k) * Eff(g, k);
            s2[g] = acc / dfr;
        }
        if (standardize == "p2") {
            SVResult sv = squeeze_var(s2r, 0.92 * dfr);
            s2r = sv.var_post;
        } else if (standardize == "posterior.sd") {
            double scale, df2;
            fit_fdist(s2, dfr, scale, df2);
            if (!std::isfinite(df2) || df2 > 1e100) {
                for (int g = 0; g < G; ++g) s2r[g] = scale;
            } else {
                for (int g = 0; g < G; ++g)
                    s2r[g] = (0.92 * dfr * s2r[g] + df2 * scale) /
                        (0.92 * dfr + df2);
            }
        }
        for (int g = 0; g < G; ++g) {
            double f = std::sqrt(s2r[g]);
            for (int k = 0; k < ne; ++k) Eff(g, k) /= f;
        }
    }

    std::vector<double> gw;
    if (gene_weights.isNotNull()) {
        gw = Rcpp::as<std::vector<double>>(gene_weights);
        if ((int)gw.size() != G)
            Rcpp::stop("gene.weights vector should be of length nrow(E)");
    }

    std::vector<std::vector<int>> sets = as_index(index);
    int nsets = (int)sets.size();
    Rcpp::IntegerVector ngenes(nsets);
    Rcpp::NumericVector tstat(nsets), pmixed(nsets);
    for (int i = 0; i < nsets; ++i) {
        const std::vector<int>& iset = sets[i];
        int m = (int)iset.size();
        ngenes[i] = m;
        std::vector<double> ES((size_t)m * ne);
        for (int a = 0; a < m; ++a)
            for (int k = 0; k < ne; ++k)
                ES[(size_t)a * ne + k] = Eff(iset[a], k);
        if (!gw.empty())
            for (int a = 0; a < m; ++a)
                for (int k = 0; k < ne; ++k)
                    ES[(size_t)a * ne + k] *= gw[iset[a]];
        // set-level t statistic
        std::vector<double> M(ne, 0.0);
        for (int a = 0; a < m; ++a)
            for (int k = 0; k < ne; ++k) M[k] += ES[(size_t)a * ne + k];
        for (int k = 0; k < ne; ++k) M[k] /= m;
        double vacc = 0.0;
        for (int k = 1; k < ne; ++k) vacc += M[k] * M[k];
        tstat[i] = M[0] / std::sqrt(vacc / (ne - 1.0));
        // mixed p-value via the rotated beta distribution
        if (m > 1) {
            std::vector<double> ESc((size_t)m * ne);  // col-major for the SVD
            for (int a = 0; a < m; ++a)
                for (int k = 0; k < ne; ++k)
                    ESc[(size_t)k * m + a] = ES[(size_t)a * ne + k];
            std::vector<double> A = svd_values(ESc, m, ne);
            for (double& v : A) v *= v;  // limma uses SVD$d^2
            int d1 = std::min(m, ne);
            double c0 = 0.0;
            for (int a = 0; a < m; ++a) c0 += ES[(size_t)a * ne] * ES[(size_t)a * ne];
            double a_max = A[0], a_min = A[d1 - 1];
            double Fobs = (c0 - a_min) / (a_max - a_min);
            double asum = 0.0, asq = 0.0;
            for (int j = 0; j < d1; ++j) { asum += A[j]; asq += A[j] * A[j]; }
            double d = d1 - 1.0;
            double bmean = 1.0 / d1;
            double bvar = ((d / d1) / d1) / (d1 / 2.0 + 1.0);
            double Fm = (asum * bmean - a_min) / (a_max - a_min);
            // A^T COV A with COV = bvar*I - (bvar/d)*(J - I):
            // bvar*(1 + 1/d)*sum(A^2) - (bvar/d)*sum(A)^2; this form is
            // guaranteed >= 0, a dropped term makes it go negative when
            // singular values are nearly equal
            double Fv = (bvar * (1.0 + 1.0 / d) * asq -
                (bvar / d) * asum * asum) /
                ((a_max - a_min) * (a_max - a_min));
            double apb = Fm * (1.0 - Fm) / Fv - 1.0;
            double alpha = apb * Fm;
            double beta = apb - alpha;
            pmixed[i] = Rf_pbeta(Fobs, alpha, beta, 0, 0);
        } else {
            pmixed[i] = NA_REAL;
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("t.stat") = tstat,
        Rcpp::Named("ngenes") = ngenes,
        Rcpp::Named("pvalue.mixed") = pmixed,
        Rcpp::Named("df.residual") = dfr);
}
