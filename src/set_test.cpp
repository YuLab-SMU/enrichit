// set_test.cpp -- effects layer for the native dsea engine
//
// This file holds the linear-model effects computation shared by camera /
// roast / fry, plus the limma-compatible variance moderation used by fry's
// "posterior.sd" standardisation.
//
// Design note (do not "optimise" this away): the QR decomposition is delegated
// to R's own qr() / qr.qty() rather than implemented here.  A previous
// hand-written Householder QR reproduced R's dqrdc2 only to ~14 digits, and
// the last pivot of a per-gene *weighted* design is numerically zero as a
// matter of course; a rounding difference there flips the sign of that pivot,
// which reflects the whole residual basis of Q'y.  limma's .fryEffects /
// .roastEffects consume those residual coordinates through colMeans(), which
// is not a basis-invariant quantity, so the flipped basis silently changed
// p-values.  Delegating to R keeps the basis bit-identical.
//
// Do NOT reintroduce a local Householder / BLAS-nrm2 path.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <vector>

// ---------------------------------------------------------------------------
// R functions, resolved once and handed out by reference.
//
// Returning `Rcpp::Function` by value would copy the function (and its
// environment) on every single gene in the weighted loop, so these return
// references into a function-local static.
// ---------------------------------------------------------------------------

static const Rcpp::Function& base_qr() {
    static Rcpp::Function f = Rcpp::Function("qr");
    return f;
}

static const Rcpp::Function& base_qr_qty() {
    static Rcpp::Function f = Rcpp::Function("qr.qty");
    return f;
}

// R's qr() reports rank deficiency through the $rank component instead of
// raising, so a full-rank design has to be checked explicitly.
static void check_qr_rank(const Rcpp::List& qrX, int p) {
    Rcpp::IntegerVector rank = qrX["rank"];
    if (rank.size() < 1 || rank[0] < p)
        Rcpp::stop("design matrix is not of full rank");
}

// ---------------------------------------------------------------------------
// Variance moderation: limma::fitFDist + limma::.squeezeVar, reimplemented.
// ---------------------------------------------------------------------------

namespace {

double st_trigamma(double x) {
    return R::psigamma(x, 1.0);
}

// limma's logmdigamma (from statmod): log(gamma(x)) - log(x), i.e. the
// log-derivative of lgamma, computed with the standard asymptotic series.
double st_logmdigamma(double x) {
    if (R_IsNA(x) || x <= 0.0)
        return NA_REAL;
    if (x < 5.0)
        return std::log(x / (x + 5.0)) + st_logmdigamma(x + 5.0) + 1.0 / x +
               1.0 / (x + 1.0) + 1.0 / (x + 2.0) + 1.0 / (x + 3.0) +
               1.0 / (x + 4.0);
    // Asymptotic series in z = 1/x^2.  The series variable must be z, not x:
    // using x here made the tail astronomically large and returned values off
    // by ~6 orders of magnitude at every argument (e.g. x = 1 gave -722107.2
    // instead of 0.5772156649), which fed garbage scale/df2 into fitFDist and
    // made fry's "posterior.sd" return NaN and roast's active proportions 0.
    // Checked against statmod::logmdigamma: max|diff| = 0 at
    // x = 0.5, 1, 2.5, 5, 10, 14.9, 50, 100.
    double z = 1.0 / (x * x);
    double tail = z * (-1.0 / 12.0 +
                       z * (1.0 / 120.0 +
                            z * (-1.0 / 252.0 +
                                 z * (1.0 / 240.0 +
                                      z * (-1.0 / 132.0 +
                                           z * (691.0 / 32760.0 +
                                                z * (-1.0 / 12.0 +
                                                     z * (3617.0 / 8160.0))))))));
    return 1.0 / (2.0 * x) - tail;
}

// limma's trigammaInverse: solve trigamma(y) = x by Newton iteration.
double st_trigamma_inverse(double x) {
    if (R_IsNaN(x) || x < 0.0)
        return R_NaN;
    if (x > 1e7)
        return 1.0 / std::sqrt(x);
    if (x < 1e-6)
        return 1.0 / x;
    double y = 0.5 + 1.0 / x;
    for (int iter = 0; iter < 50; ++iter) {
        double tri = st_trigamma(y);
        double dif = tri * (1.0 - tri / x) / R::psigamma(y, 2.0);
        y += dif;
        if (std::fabs(dif / y) < 1e-8)
            break;
    }
    return y;
}

// Mirrors limma::fitFDist for the unweighted (covariate = NULL) case, which
// is the only case the fry standardisation needs.
void st_fit_fdist(const std::vector<double>& x_in, double df1, double& scale,
                  double& df2) {
    const std::size_t n = x_in.size();
    if (n == 0) {
        scale = NA_REAL;
        df2 = NA_REAL;
        return;
    }
    if (n == 1) {
        scale = x_in[0];
        df2 = 0.0;
        return;
    }
    bool ok1 = R_FINITE(df1) && df1 > 1e-15;
    if (!ok1) {
        scale = NA_REAL;
        df2 = NA_REAL;
        return;
    }

    std::vector<double> x;
    x.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        double v = x_in[i];
        if (R_FINITE(v) && v > -1e-15)
            x.push_back(v);
    }
    const std::size_t nok = x.size();
    if (nok == 0) {
        scale = NA_REAL;
        df2 = NA_REAL;
        return;
    }
    if (nok == 1) {
        scale = x[0];
        df2 = 0.0;
        return;
    }

    for (std::size_t i = 0; i < nok; ++i)
        x[i] = std::max(x[i], 0.0);

    std::vector<double> s(x);
    std::sort(s.begin(), s.end());
    double m = (nok % 2 == 1) ? s[nok / 2]
                              : 0.5 * (s[nok / 2 - 1] + s[nok / 2]);
    if (m == 0.0)
        m = 1.0;  // limma warns and substitutes 1 here as well
    for (std::size_t i = 0; i < nok; ++i)
        x[i] = std::max(x[i], 1e-5 * m);

    std::vector<double> e(nok);
    double emean = 0.0;
    for (std::size_t i = 0; i < nok; ++i) {
        e[i] = std::log(x[i]) + st_logmdigamma(df1 / 2.0);
        emean += e[i];
    }
    emean /= static_cast<double>(nok);

    double evar = 0.0;
    for (std::size_t i = 0; i < nok; ++i) {
        double d = e[i] - emean;
        evar += d * d;
    }
    evar /= static_cast<double>(nok - 1);
    evar -= st_trigamma(df1 / 2.0);

    if (evar > 0.0) {
        df2 = 2.0 * st_trigamma_inverse(evar);
        scale = std::exp(emean - st_logmdigamma(df2 / 2.0));
    } else {
        df2 = R_PosInf;
        double sum = 0.0;
        for (std::size_t i = 0; i < nok; ++i)
            sum += x[i];
        scale = sum / static_cast<double>(nok);
    }
}

// limma::.squeezeVar
void st_squeeze_var(const std::vector<double>& var, double df,
                    const std::vector<double>& var_prior,
                    const std::vector<double>& df_prior,
                    std::vector<double>& var_post) {
    const std::size_t n = var.size();
    var_post.assign(n, 0.0);

    // df.prior is either length 1 (fitFDist) or length n.
    const bool prior_len1 = df_prior.size() == 1;
    const bool prior_len_n = df_prior.size() == n;
    if (!prior_len1 && !prior_len_n)
        Rcpp::stop("df.prior length does not match var");

    double m = prior_len1 ? df_prior[0] : *std::max_element(df_prior.begin(), df_prior.end());
    if (R_FINITE(m)) {
        for (std::size_t i = 0; i < n; ++i) {
            double vp = var_prior.size() == 1 ? var_prior[0] : var_prior[i];
            double dp = prior_len1 ? df_prior[0] : df_prior[i];
            var_post[i] = (df * var[i] + dp * vp) / (df + dp);
        }
        return;
    }

    for (std::size_t i = 0; i < n; ++i)
        var_post[i] = var_prior.size() == 1 ? var_prior[0] : var_prior[i];

    m = prior_len1 ? df_prior[0] : *std::min_element(df_prior.begin(), df_prior.end());
    if (m > 1e100)
        return;

    for (std::size_t i = 0; i < n; ++i) {
        double dp = prior_len1 ? df_prior[0] : df_prior[i];
        if (!R_FINITE(dp))
            continue;
        var_post[i] = (df * var[i] + dp * var_post[i]) / (df + dp);
    }
}

}  // namespace

//' Moderated variance estimation (limma-compatible)
//'
//' Internal counterpart of \code{limma::squeezeVar(var, df)}: an empirical
//' Bayes moderation that borrows strength across genes via
//' \code{fitFDist} + \code{.squeezeVar}.
//'
//' @param var numeric vector of per-gene sample variances.
//' @param df residual degrees of freedom (length 1 or \code{length(var)}).
//' @return a list with \code{var.post}, \code{var.prior} and \code{df.prior},
//'   matching the return value of \code{limma::squeezeVar}.
//' @noRd
// [[Rcpp::export]]
Rcpp::List st_squeeze_var_cpp(Rcpp::NumericVector var_in, double df) {
    std::vector<double> var(var_in.begin(), var_in.end());
    const std::size_t n = var.size();
    std::vector<double> var_post;
    double var_prior = NA_REAL;
    double df_prior = NA_REAL;

    if (n == 0)
        Rcpp::stop("var is empty");
    if (n < 3) {
        // limma::squeezeVar short-circuits small inputs.
        return Rcpp::List::create(Rcpp::Named("var.post") = var_in,
                                  Rcpp::Named("var.prior") = var_in,
                                  Rcpp::Named("df.prior") = 0.0);
    }

    // Prior hyper-parameters from the package's own fitFDist reimplementation.
    // This must NOT call limma: REBUILD_SPEC 1 requires limma to be absent from
    // Imports/Suggests and forbids limma:: inside the package, and reaching for
    // namespace_env("limma") at run time would fail on any machine without it.
    st_fit_fdist(var, df, var_prior, df_prior);
    if (R_IsNA(df_prior))
        Rcpp::stop("Could not estimate prior df");

    std::vector<double> vp_vec(1, var_prior);
    std::vector<double> dp_vec(1, df_prior);
    st_squeeze_var(var, df, vp_vec, dp_vec, var_post);

    // limma::squeezeVar returns var.post carrying the input's names.
    Rcpp::NumericVector post = Rcpp::wrap(var_post);
    if (var_in.hasAttribute("names"))
        post.attr("names") = var_in.attr("names");

    return Rcpp::List::create(Rcpp::Named("var.post") = post,
                              Rcpp::Named("var.prior") = Rcpp::wrap(var_prior),
                              Rcpp::Named("df.prior") = Rcpp::wrap(df_prior));
}

//' Linear-model effects: contrast coordinate + residual coordinates
//'
//' Shared effects computation for camera / roast / fry.  Returns a
//' \code{G x (n - p + 1)} matrix whose first column is the contrast
//' coordinate (signed) and whose remaining \code{n - p} columns are the
//' residual coordinates of \eqn{Q'y}.
//'
//' The QR decomposition is delegated to R's \code{qr()} / \code{qr.qty()} so
//' that the residual basis is bit-identical to limma's \code{.lmEffects};
//' see the file header for why that matters.
//'
//' @param E G x n matrix of log-expression values (finite, no NA).
//' @param weights optional G x n matrix of positive precision weights.
//' @param design n x p design matrix.
//' @param contrast_idx 1-based index of the contrast column of \code{design}.
//' @return G x (n - p + 1) numeric matrix.
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix st_effects_cpp(Rcpp::NumericMatrix E,
                                   Rcpp::Nullable<Rcpp::NumericMatrix> weights,
                                   Rcpp::NumericMatrix design,
                                   int contrast_idx) {
    const int G = E.nrow();
    const int n = E.ncol();
    const int p = design.ncol();

    if (design.nrow() != n)
        Rcpp::stop("row dimension of design matrix must match column dimension of data");
    if (n <= p)
        Rcpp::stop("No residual degrees of freedom");
    if (contrast_idx < 1 || contrast_idx > p)
        Rcpp::stop("contrast_idx is out of range");

    // Every response value must be finite.
    for (int j = 0; j < n; ++j) {
        const double* col = E.begin() + static_cast<R_xlen_t>(j) * G;
        for (int g = 0; g < G; ++g) {
            if (!R_FINITE(col[g]))
                Rcpp::stop("All y values must be finite and non-NA");
        }
    }

    const bool weighted = weights.isNotNull();
    Rcpp::NumericMatrix W;
    if (weighted) {
        W = Rcpp::NumericMatrix(weights);
        if (W.nrow() != G || W.ncol() != n)
            Rcpp::stop("weights must have same dimensions as y");
        for (int j = 0; j < n; ++j) {
            const double* col = W.begin() + static_cast<R_xlen_t>(j) * G;
            for (int g = 0; g < G; ++g) {
                if (!R_FINITE(col[g]) || col[g] <= 0.0)
                    Rcpp::stop("weights must be positive");
            }
        }
    }

    // Put the contrast column last; the remaining columns keep their
    // relative order (limma .lmEffects reorders X the same way).
    const int ci = contrast_idx - 1;
    Rcpp::NumericMatrix X(n, p);
    {
        int k = 0;
        for (int j = 0; j < p; ++j) {
            if (j == ci)
                continue;
            for (int i = 0; i < n; ++i)
                X(i, k) = design(i, j);
            ++k;
        }
        for (int i = 0; i < n; ++i)
            X(i, p - 1) = design(i, ci);
    }

    const int ncol_out = n - p + 1;
    Rcpp::NumericMatrix Eff(G, ncol_out);
    const int i0 = p - 1;  // 0-based index of the contrast coordinate in qy

    // Column names follow limma .lmEffects: the contrast-coefficient name
    // followed by the residual degrees of freedom's row labels ("p", "p+1", ...).
    // rownames are carried over from E (limma falls back to 1:G when absent).
    {
        Rcpp::CharacterVector cn(ncol_out);
        Rcpp::CharacterVector dxn = Rcpp::colnames(design);
        Rcpp::String cname = (dxn.size() > ci && dxn[ci] != NA_STRING)
                                 ? dxn[ci]
                                 : Rcpp::String("Contrast");
        cn[0] = cname;
        for (int k = 1; k < ncol_out; ++k)
            cn[k] = std::to_string(i0 + k + 1);
        Rcpp::colnames(Eff) = cn;

        Rcpp::CharacterVector rn = Rcpp::rownames(E);
        if (rn.size() == G) {
            Rcpp::rownames(Eff) = rn;
        } else {
            Rcpp::CharacterVector idx(G);
            for (int g = 0; g < G; ++g)
                idx[g] = std::to_string(g + 1);
            Rcpp::rownames(Eff) = idx;
        }
    }

    if (!weighted) {
        Rcpp::List qrX = base_qr()(X);
        check_qr_rank(qrX, p);

        Rcpp::NumericMatrix qr_mat = qrX["qr"];
        const double s = (qr_mat(p - 1, p - 1) < 0.0) ? -1.0 : 1.0;

        Rcpp::NumericMatrix qy = base_qr_qty()(qrX, Rcpp::transpose(E));
        for (int g = 0; g < G; ++g) {
            Eff(g, 0) = s * qy(i0, g);
            for (int k = 1; k < ncol_out; ++k)
                Eff(g, k) = qy(i0 + k, g);
        }
        return Eff;
    }

    // Weighted branch: one QR per gene.  The response passed to qr.qty() is a
    // single column (n x 1) -- declaring it n x G would recompute qr.qty() for
    // every gene on every iteration, which is the O(G^2) bug documented in
    // the rebuild spec (3.16 s at G = 4000 vs 0.111 s once narrowed).
    std::vector<double> wX(n * p);
    std::vector<double> wy(n);
    std::vector<double> sw(n);
    Rcpp::NumericMatrix wXmat(n, p);
    Rcpp::NumericMatrix wycol(n, 1);

    for (int g = 0; g < G; ++g) {
        const double* wcol = W.begin() + static_cast<R_xlen_t>(g);
        for (int i = 0; i < n; ++i)
            sw[i] = std::sqrt(wcol[static_cast<R_xlen_t>(i) * G]);

        for (int j = 0; j < p; ++j) {
            double* dst = &wX[static_cast<std::size_t>(j) * n];
            for (int i = 0; i < n; ++i)
                dst[i] = X(i, j) * sw[i];
        }
        std::copy(wX.begin(), wX.end(), wXmat.begin());

        for (int i = 0; i < n; ++i)
            wy[i] = E(g, i) * sw[i];
        std::copy(wy.begin(), wy.end(), wycol.begin());

        Rcpp::List qrX = base_qr()(wXmat);
        check_qr_rank(qrX, p);

        Rcpp::NumericMatrix qr_mat = qrX["qr"];
        const double s = (qr_mat(p - 1, p - 1) < 0.0) ? -1.0 : 1.0;

        Rcpp::NumericVector qy = base_qr_qty()(qrX, wycol);
        Eff(g, 0) = s * qy[i0];
        for (int k = 1; k < ncol_out; ++k)
            Eff(g, k) = qy[i0 + k];
    }

    return Eff;
}

// ===========================================================================
//  camera / roast / fry  -- the three set-level statistics
//
//  All three consume the same lm_effects layer above (st_effects_cpp), so the
//  QR decomposition exists in exactly one place.  Nothing here re-inlines a
//  QR and nothing here calls limma at run time.
// ===========================================================================

namespace {

// ---------------------------------------------------------------------------
// 128-node Gauss-Legendre rule on (0,1), identical to
// statmod::gauss.quad.prob(128, "uniform").
//
// fry's "posterior.sd" standardisation needs
//     Eu2max = sum_i (dfr + 1) * x_i^dfr * qchisq(x_i, 1) * w_i
// and limma obtains it by quadrature.  statmod is not a dependency of this
// package, so the rule is embedded as a constant, printed with %.17g so it
// round-trips exactly.  Checked: Eu2max is identical() to limma's for
// dfr = 1, 2, 3, 5, 10, 20, 50.
//
// An in-house Golub-Welsch eigensolve was tried and rejected: it reproduced the
// nodes to ~5e-16 but moved Eu2max by 30-60 ulp (~1e-14 relative), which is the
// entire parity budget.  A constant table has no such error.
//
// The nodes and weights are stored PAIRED and in ascending node order.  They
// must not be reordered independently: qchisq() is not symmetric, so sorting
// the weights (or keeping statmod's own interleaved order) silently mispairs
// them and inflates Eu2max -- an earlier revision of this table had ascending
// nodes against unsorted weights, which made Eu2max(5) 11.22 instead of 3.03
// and turned fry's "posterior.sd" into NaN and roast's active proportions into
// 0.  Verified: Eu2max is identical() to limma's for dfr = 1..8 and 20/50.
// ---------------------------------------------------------------------------

const double ST_GQ_NODES[128] = {
    8.755602643395477e-05, 0.00046127001131196499, 0.0011333756872426437, 0.0021036207325093037,
    0.0033714435498936046, 0.0049360907541327603, 0.0067966286377065255, 0.0089519457821403137,
    0.011400754268045932, 0.014141590626431444, 0.017172816784017331, 0.020492621073150141,
    0.024099019329367921, 0.027989856084889853, 0.03216280586104181, 0.036615374560526159,
    0.04134490095951987, 0.046348558299121878, 0.05162335597542117, 0.057166141327301689,
    0.062973601520984079, 0.069042265530225744, 0.075368506211015596, 0.081948542469546526,
    0.088778441522178175, 0.095854121246043145, 0.10317135261890342, 0.11072576224679404,
    0.11851283497795256, 0.12652791660146889, 0.13476621662904559, 0.14322281115820634,
    0.15189264581524298, 0.16077053877614023, 0.16985118386367676, 0.17912915371884602,
    0.18859890304470756, 0.19825477192072588, 0.20809098918561886, 0.21810167588669116,
    0.22828084879359495, 0.23862242397441236, 0.24912022043192772, 0.25976796379791389,
    0.2705592900832241, 0.28148774948144784, 0.29254681022386259, 0.30372986248336631,
    0.31503022232507061, 0.32644113570118249, 0.33795578248779334, 0.34956728056116126,
    0.36126868991104777, 0.37305301678865266, 0.38491321788666999, 0.39684220454896052,
    0.40883284700733152, 0.42087797864288773, 0.43297040026940631, 0.4451028844361784,
    0.45726817974774248, 0.46945901519793037, 0.48166810451563324, 0.49388815051969204,
    0.50611184948030785, 0.51833189548436687, 0.5305409848020699, 0.54273182025225764,
    0.55489711556382204, 0.56702959973059419, 0.57912202135711277, 0.59116715299266875,
    0.60315779545103965, 0.61508678211333023, 0.62694698321134723, 0.63873131008895223,
    0.65043271943883851, 0.66204421751220666, 0.67355886429881773, 0.68496977767492961,
    0.6962701375166338, 0.70745318977613736, 0.71851225051855194, 0.72944070991677612,
    0.74023203620208577, 0.75087977956807217, 0.76137757602558764, 0.77171915120640522,
    0.78189832411330884, 0.79190901081438159, 0.80174522807927406, 0.81140109695529228,
    0.82087084628115359, 0.83014881613632285, 0.83922946122385955, 0.84810735418475702,
    0.85677718884179366, 0.86523378337095425, 0.87347208339853077, 0.88148716502204716,
    0.88927423775320591, 0.89682864738109658, 0.90414587875395669, 0.91122155847782182,
    0.91805145753045325, 0.92463149378898435, 0.93095773446977426, 0.93702639847901592,
    0.94283385867269831, 0.94837664402457911, 0.95365144170087834, 0.95865509904048041,
    0.96338462543947379, 0.96783719413895808, 0.97201014391510987, 0.97590098067063225,
    0.9795073789268498, 0.9828271832159825, 0.98585840937356806, 0.98859924573195346,
    0.99104805421785902, 0.99320337136229297, 0.99506390924586707, 0.99662855645010628,
    0.99789637926749042, 0.99886662431275686, 0.99953872998868776, 0.99991244397356582
};

const double ST_GQ_WEIGHTS[128] = {
    0.00022469048014607816, 0.00052290633967003446, 0.00082125150933456633, 0.0011191442154816114,
    0.0014163757357288135, 0.001712763020455342, 0.0020081274918690678, 0.0023022921283512954,
    0.0025950809163382198, 0.0028863187714330287, 0.003175831580853568, 0.0034634462834493623,
    0.0037489909628173283, 0.0040322949452431941, 0.0043131888993085685, 0.0045915049358300521,
    0.0048670767075032172, 0.0051397395079164293, 0.0054093303697514868, 0.0056756881620397132,
    0.0059386536863703151, 0.0061980697719755664, 0.0064537813696336887, 0.0067056356443079531,
    0.0069534820664761861, 0.0071971725020834121, 0.0074365613010736161, 0.0076715053844324428,
    0.0079018643296995521, 0.0081275004548928571, 0.0083482789007945368, 0.0085640677115556944,
    0.0087747379135589536, 0.0089801635925041614, 0.0091802219686654661, 0.0093747934702723881,
    0.009563761804975526, 0.0097470140293533791, 0.0099244406164160689, 0.01009593552106465,
    0.010261396243479942, 0.010420723890375377, 0.010573823234110637, 0.010720602769604395,
    0.010860974769026085, 0.01099485533423019, 0.011122164446900083, 0.011242826016372678,
    0.011356767925118439, 0.011463922071843304, 0.011564224412193328, 0.011657614997031541,
    0.011744038008268251, 0.011823441792223863, 0.01189577889050168, 0.011961006068351345,
    0.012019084340512003, 0.012069978994509489, 0.012113659611408046, 0.012150100083985907,
    0.012179278632345318, 0.012201177816924671, 0.012215784548924713, 0.012223090098131002,
    0.012223090098131367, 0.012215784548925033, 0.012201177816924702, 0.012179278632345383,
    0.012150100083985767, 0.012113659611407946, 0.012069978994509456, 0.01201908434051225,
    0.011961006068351931, 0.011895778890501876, 0.011823441792223558, 0.011744038008268303,
    0.011657614997031314, 0.011564224412193005, 0.011463922071843546, 0.01135676792511814,
    0.011242826016372758, 0.011122164446899811, 0.010994855334229916, 0.010860974769025365,
    0.010720602769604072, 0.010573823234110529, 0.01042072389037613, 0.010261396243479692,
    0.010095935521065141, 0.0099244406164152969, 0.0097470140293535682, 0.0095637618049755121,
    0.0093747934702728634, 0.0091802219686658044, 0.0089801635925047044, 0.0087747379135589901,
    0.0085640677115550214, 0.008348278900794941, 0.0081275004548930514, 0.0079018643297002582,
    0.0076715053844321288, 0.007436561301073322, 0.0071971725020836342, 0.0069534820664760447,
    0.0067056356443080008, 0.0064537813696333652, 0.0061980697719753132, 0.0059386536863700306,
    0.0056756881620399283, 0.0054093303697512378, 0.0051397395079161509, 0.0048670767075026794,
    0.0045915049358308301, 0.004313188899308853, 0.0040322949452424135, 0.003748990962817624,
    0.0034634462834498784, 0.0031758315808536868, 0.002886318771433069, 0.0025950809163379942,
    0.0023022921283511349, 0.0020081274918689889, 0.0017127630204551692, 0.0014163757357290752,
    0.0011191442154814828, 0.00082125150933481548, 0.00052290633967029261, 0.0002246904801457223
};

// Eu2max exactly as limma evaluates it.
double st_eu2max(int dfr) {
    double out = 0.0;
    for (int i = 0; i < 128; ++i) {
        const double x = ST_GQ_NODES[i];
        out += (static_cast<double>(dfr) + 1.0) * std::pow(x, dfr) *
               R::qchisq(x, 1.0, 1, 0) * ST_GQ_WEIGHTS[i];
    }
    return out;
}

// limma::zscoreT(x, df, approx = TRUE, method = "bailey") -- used by
// .roastEffects for both the observed and the rotated statistics.
inline double st_zscore_t_bailey(double x, double df) {
    const double sgn = (x > 0.0) - (x < 0.0);
    return ((df + 0.125) / (df + 1.125)) *
           std::sqrt((df + 19.0 / 12.0) * std::log1p(x / (df + 1.0 / 12.0) * x)) * sgn;
}

// limma::zscoreT(x, df, approx = TRUE, method = "hill") -- used by camera.
// (limma::camera calls zscoreT(modt, df, approx = TRUE, method = "hill");
// using bailey here instead shifts p-values by ~2.5e-4, so the method matters.)
inline double st_zscore_t_hill(double x, double df) {
    const double A = df - 0.5;
    const double B = 48.0 * A * A;
    double z = A * std::log1p(x / df * x);
    z = (((((-0.4 * z - 3.3) * z - 24.0) * z - 85.5) / (0.8 * z * z + 100.0 + B) +
          z + 3.0) / B + 1.0) * std::sqrt(z);
    return z * ((x > 0.0) - (x < 0.0));
}

// Shared per-set plumbing: resolve the index list's gene ids to 0-based rows,
// validating them the way limma's set loops do.
void st_resolve_set(const Rcpp::List& index, int i, int G, std::vector<int>& rows) {
    Rcpp::IntegerVector iset = index[i];
    const int m = iset.size();
    rows.resize(m);
    for (int j = 0; j < m; ++j) {
        const int gi = iset[j] - 1;  // R is 1-based, C++ 0-based
        if (gi < 0 || gi >= G)
            Rcpp::stop("gene set %d contains an index outside 1..nrow(E)", i + 1);
        rows[j] = gi;
    }
}

// Argument checking shared by the three entry points.
void st_check_common(const Rcpp::NumericMatrix& E, const Rcpp::NumericMatrix& design,
                     int contrast_idx, const Rcpp::List& index) {
    const int n = E.ncol();
    const int p = design.ncol();
    if (design.nrow() != n)
        Rcpp::stop("row dimension of design matrix must match column dimension of data");
    if (n <= p)
        Rcpp::stop("No residual degrees of freedom");
    if (contrast_idx < 1 || contrast_idx > p)
        Rcpp::stop("contrast_idx is out of range");
    if (index.size() < 1)
        Rcpp::stop("index is empty");
}

}  // namespace

//' Competitive gene-set test with correlation adjustment (camera)
//'
//' Native counterpart of \code{limma::camera}: a competitive test of whether a
//' set's genes are differentially expressed relative to the rest of the
//' genome, with the two-sample t-statistic's variance inflated by an
//' inter-gene correlation.  The correlation is either supplied by the caller
//' or estimated from the residual coordinates.
//'
//' @param E G x n matrix of log-expression values.
//' @param weights optional G x n matrix of positive precision weights.
//' @param design n x p design matrix.
//' @param contrast_idx 1-based index of the contrast column of \code{design}.
//' @param index list of integer gene-set memberships (1-based row numbers).
//' @param inter_gene_cor fixed inter-gene correlation, or \code{NULL} to
//'   estimate it per set.
//' @param allow_neg_cor allow estimated negative correlations to be used;
//'   when \code{FALSE} the variance inflation is floored at 1.
//' @return a list with \code{NGenes}, \code{Correlation}, \code{Down},
//'   \code{Up} and \code{TwoSided}; the two-sided p-value is
//'   \code{2 * min(Down, Up)}.
//' @noRd
// [[Rcpp::export]]
Rcpp::List st_camera_cpp(Rcpp::NumericMatrix E,
                         Rcpp::Nullable<Rcpp::NumericMatrix> weights,
                         Rcpp::NumericMatrix design, int contrast_idx,
                         Rcpp::List index,
                         Rcpp::Nullable<double> inter_gene_cor,
                         bool allow_neg_cor) {
    st_check_common(E, design, contrast_idx, index);

    const int G = E.nrow();
    const int n = E.ncol();
    const int p = design.ncol();
    const int dfr = n - p;

    const bool fixed_cor = inter_gene_cor.isNotNull();
    const double df_camera = fixed_cor ? static_cast<double>(G - 2)
                                       : std::min(static_cast<double>(dfr),
                                                  static_cast<double>(G - 2));

    Rcpp::NumericMatrix Eff = st_effects_cpp(E, weights, design, contrast_idx);
    const int nres = n - p;

    // Residual standard deviations per gene (limma: sigma2 <- colMeans(U^2) on
    // the transposed effects matrix, i.e. the mean over residual coordinates).
    std::vector<double> sigma2(G), modt(G);
    for (int g = 0; g < G; ++g) {
        double s = 0.0;
        for (int k = 1; k <= nres; ++k)
            s += Eff(g, k) * Eff(g, k);
        sigma2[g] = s / nres;
    }

    Rcpp::List sv = st_squeeze_var_cpp(Rcpp::wrap(sigma2),
                                       static_cast<double>(dfr));
    Rcpp::NumericVector var_post = Rcpp::as<Rcpp::NumericVector>(sv["var.post"]);
    const double df_prior = Rcpp::as<double>(sv["df.prior"]);

    for (int g = 0; g < G; ++g)
        modt[g] = Eff(g, 0) / std::sqrt(var_post[g]);

    // limma: df.total <- min(df.residual + sv$df.prior, G * df.residual);
    //        Stat <- zscoreT(modt, df = df.total, approx = TRUE, method = "hill")
    const double df_total = std::min(static_cast<double>(dfr) + df_prior,
                                     static_cast<double>(G) * dfr);
    const double df_tw = std::min(df_total, 10000.0);

    std::vector<double> Stat(G);
    for (int g = 0; g < G; ++g)
        Stat[g] = st_zscore_t_hill(modt[g], df_tw);

    double meanStat = 0.0;
    for (int g = 0; g < G; ++g)
        meanStat += Stat[g];
    meanStat /= G;
    double varStat = 0.0;
    for (int g = 0; g < G; ++g) {
        const double d = Stat[g] - meanStat;
        varStat += d * d;
    }
    varStat /= (G - 1);

    // U is the residual block rescaled to unit variance; the estimated
    // correlation reads off the column means of U restricted to the set.
    std::vector<double> U(static_cast<std::size_t>(G) * nres);
    for (int g = 0; g < G; ++g) {
        const double sd = std::sqrt(std::max(sigma2[g], 1e-8));
        for (int k = 0; k < nres; ++k)
            U[static_cast<std::size_t>(g) * nres + k] = Eff(g, k + 1) / sd;
    }

    const int nsets = index.size();
    Rcpp::NumericVector NGenes(nsets), Correlation(nsets), Down(nsets), Up(nsets),
        TwoSided(nsets);
    std::vector<int> rows;

    for (int i = 0; i < nsets; ++i) {
        st_resolve_set(index, i, G, rows);
        const int m = static_cast<int>(rows.size());
        const int m2 = G - m;
        if (m2 < 1)
            Rcpp::stop("gene set %d contains every gene; the competitive test needs a complement",
                       i + 1);

        double correlation, vif;
        if (fixed_cor) {
            correlation = Rcpp::as<double>(inter_gene_cor);
            vif = 1.0 + (m - 1.0) * correlation;
        } else if (m > 1) {
            // limma: vif <- m * mean(colMeans(Uset)^2)  -- the mean is taken
            // over the n - p residual coordinates, not summed.
            double s2tot = 0.0;
            for (int k = 0; k < nres; ++k) {
                double cm = 0.0;
                for (int j = 0; j < m; ++j)
                    cm += U[static_cast<std::size_t>(rows[j]) * nres + k];
                cm /= m;
                s2tot += cm * cm;
            }
            vif = m * (s2tot / nres);
            correlation = (vif - 1.0) / (m - 1.0);
        } else {
            vif = 1.0;
            correlation = NA_REAL;
        }

        NGenes[i] = m;
        Correlation[i] = correlation;

        // limma: if (!allow.neg.cor) vif <- max(1, vif)
        if (!allow_neg_cor)
            vif = std::max(1.0, vif);

        double meanStatInSet = 0.0;
        for (int j = 0; j < m; ++j)
            meanStatInSet += Stat[rows[j]];
        meanStatInSet /= m;

        const double delta = static_cast<double>(G) / m2 * (meanStatInSet - meanStat);
        const double varStatPooled =
            ((G - 1) * varStat - delta * delta * m * m2 / static_cast<double>(G)) / (G - 2);
        const double two_sample_t =
            delta / std::sqrt(varStatPooled * (vif / m + 1.0 / m2));

        Down[i] = R::pt(two_sample_t, df_camera, 1, 0);
        Up[i] = R::pt(two_sample_t, df_camera, 0, 0);
        TwoSided[i] = 2.0 * std::min(Down[i], Up[i]);
    }

    // index.names() is NULL when the caller passes an unnamed gene-set list.
    // Assigning that straight into a CharacterVector raises "Not compatible
    // with STRSXP: [type=NULL]" and crashes the process; limma tolerates an
// unnamed list, so guard the conversion and simply omit the names.
    Rcpp::CharacterVector set_names;
    Rcpp::RObject nm_attr = index.attr("names");
    if (!Rf_isNull(nm_attr))
        set_names = Rcpp::CharacterVector(nm_attr);
    if (set_names.size() == nsets) {
        NGenes.attr("names") = set_names;
        Correlation.attr("names") = set_names;
        Down.attr("names") = set_names;
        Up.attr("names") = set_names;
        TwoSided.attr("names") = set_names;
    }

    return Rcpp::List::create(Rcpp::Named("NGenes") = NGenes,
                              Rcpp::Named("Correlation") = Correlation,
                              Rcpp::Named("Down") = Down,
                              Rcpp::Named("Up") = Up,
                              Rcpp::Named("TwoSided") = TwoSided,
                              Rcpp::Named("df.camera") = df_camera,
                              Rcpp::Named("df.prior") = df_prior);
}

//' Self-contained gene-set test on the linear-model effects (fry)
//'
//' Native counterpart of \code{limma::fry}: a rotation test that needs no
//' permutations and no correlation assumption, operating on the effects
//' computed once for all genes.
//'
//' The \code{standardize} argument selects limma's estimator for the per-gene
//' scale: \code{"posterior.sd"} (default), \code{"residual.sd"},
//' \code{"p2"} or \code{"none"}.
//'
//' @param E G x n matrix of log-expression values.
//' @param weights optional G x n matrix of positive precision weights.
//' @param design n x p design matrix.
//' @param contrast_idx 1-based index of the contrast column of \code{design}.
//' @param index list of integer gene-set memberships (1-based row numbers).
//' @param gene_weights optional numeric vector of length G.
//' @param standardize one of \code{"posterior.sd"}, \code{"residual.sd"},
//'   \code{"p2"}, \code{"none"}.
//' @return a list with \code{NGenes}, \code{t} and \code{PValue}.  Sets with
//'   single-gene sets get a real p-value, exactly as limma does: only
//'   limma's \code{PValue.Mixed} column is gated on \code{NGenes > 1}.
//' @noRd
// [[Rcpp::export]]
Rcpp::List st_fry_cpp(Rcpp::NumericMatrix E,
                      Rcpp::Nullable<Rcpp::NumericMatrix> weights,
                      Rcpp::NumericMatrix design, int contrast_idx,
                      Rcpp::List index,
                      Rcpp::Nullable<Rcpp::NumericVector> gene_weights,
                      std::string standardize) {
    st_check_common(E, design, contrast_idx, index);

    if (standardize != "posterior.sd" && standardize != "residual.sd" &&
        standardize != "p2" && standardize != "none")
        Rcpp::stop("standardize must be one of posterior.sd, residual.sd, p2, none");

    const int G = E.nrow();
    const int n = E.ncol();
    const int p = design.ncol();
    const int dfr = n - p;
    const int nres = dfr;

    Rcpp::NumericMatrix Eff = st_effects_cpp(E, weights, design, contrast_idx);

    // NOTE: gene.weights are deliberately NOT pushed through the effects layer.
    // limma::fry forwards them to .fryEffects(), which multiplies the *set rows*
    // by iw after standardisation -- a plain signed multiply.  Routing them
    // through .lmEffects instead would apply sqrt() and reject negative weights
    // ("gene.weights must be positive"), and the parity scenarios do use signed
    // weights such as c(2, -1).
    const bool has_gw = gene_weights.isNotNull();
    std::vector<double> GW;
    if (has_gw) {
        Rcpp::NumericVector gw(gene_weights);
        if (gw.size() != G)
            Rcpp::stop("length of gene.weights should equal nrow(E)");
        GW.assign(gw.begin(), gw.end());
    }

    if (standardize != "none") {
        // limma: u2max <- apply(Effects^2, 1, max)
        //        s2.robust <- (rowSums(Effects^2) - u2max) / (dfr + 1 - Eu2max)
        const double Eu = st_eu2max(dfr);
        const double denom = static_cast<double>(dfr) + 1.0 - Eu;
        std::vector<double> s2rob(G);
        for (int g = 0; g < G; ++g) {
            double s = 0.0, mx = 0.0;
            for (int k = 0; k <= nres; ++k) {
                const double v = Eff(g, k) * Eff(g, k);
                s += v;
                if (v > mx)
                    mx = v;
            }
            s2rob[g] = (s - mx) / denom;
        }

        if (standardize == "posterior.sd") {
            // limma: s2 <- rowMeans(Effects[, -1, drop = FALSE]^2)
            //        fit <- fitFDist(s2, df1 = dfr, covariate = NULL)
            //        s2.robust <- .squeezeVar(s2.robust, df = 0.92*dfr,
            //                                 var.prior = fit$scale,
            //                                 df.prior = fit$df2)
            std::vector<double> s2(G);
            for (int g = 0; g < G; ++g) {
                double s = 0.0;
                for (int k = 1; k <= nres; ++k)
                    s += Eff(g, k) * Eff(g, k);
                s2[g] = s / nres;
            }
            double scale = NA_REAL, df2 = NA_REAL;
            st_fit_fdist(s2, static_cast<double>(dfr), scale, df2);
            const double dfx = 0.92 * dfr;
            std::vector<double> vp1(1, scale), dp1(1, df2), var_post;
            st_squeeze_var(s2rob, dfx, vp1, dp1, var_post);
            for (int g = 0; g < G; ++g) {
                const double sd = std::sqrt(var_post[g]);
                for (int k = 0; k <= nres; ++k)
                    Eff(g, k) /= sd;
            }
        } else if (standardize == "residual.sd") {
            for (int g = 0; g < G; ++g) {
                const double sd = std::sqrt(s2rob[g]);
                for (int k = 0; k <= nres; ++k)
                    Eff(g, k) /= sd;
            }
        } else {  // p2
            // limma: sv <- squeezeVar(s2.robust, df = 0.92 * dfr, ...)
            //        s2.robust <- sv$var.post
            Rcpp::List sv = st_squeeze_var_cpp(Rcpp::wrap(s2rob), 0.92 * dfr);
            Rcpp::NumericVector vp = Rcpp::as<Rcpp::NumericVector>(sv["var.post"]);
            for (int g = 0; g < G; ++g) {
                const double sd = std::sqrt(vp[g]);
                for (int k = 0; k <= nres; ++k)
                    Eff(g, k) /= sd;
            }
        }
    }

    const int nsets = index.size();
    Rcpp::NumericVector NGenes(nsets), tstat(nsets), PValue(nsets);
    std::vector<int> rows;

    for (int i = 0; i < nsets; ++i) {
        st_resolve_set(index, i, G, rows);
        const int m = static_cast<int>(rows.size());
        NGenes[i] = m;

        // limma::.fryEffects computes t.stat and PValue for EVERY set,
        // including a single-gene one -- MeanEffectsSet is just that gene's
        // own effects row. The "NGenes > 1" guard in .fryEffects applies only
        // to the PValue.Mixed column, which this engine does not return
        // (REBUILD_SPEC section 4's contract is set/NGenes/Direction/PValue).
        // Empty sets are excluded upstream by dsea(), so m == 0 cannot occur
        // here; guard it anyway rather than divide by zero.
        if (m < 1) {
            tstat[i] = NA_REAL;
            PValue[i] = NA_REAL;
            continue;
        }

        std::vector<double> M(nres + 1, 0.0);
        for (int j = 0; j < m; ++j) {
            const int gi = rows[j];
            const double iw = has_gw ? GW[gi] : 1.0;
            for (int k = 0; k <= nres; ++k)
                M[k] += iw * Eff(gi, k);
        }
        for (int k = 0; k <= nres; ++k)
            M[k] /= m;

        // limma: t <- M[1] / sqrt(mean(M[-1]^2)); p <- 2 * pt(-|t|, df = dfr)
        double rss = 0.0;
        for (int k = 1; k <= nres; ++k)
            rss += M[k] * M[k];
        tstat[i] = M[0] / std::sqrt(rss / nres);
        PValue[i] = 2.0 * R::pt(-std::fabs(tstat[i]), static_cast<double>(dfr), 1, 0);
    }

    // index.names() is NULL when the caller passes an unnamed gene-set list.
    // Assigning that straight into a CharacterVector raises "Not compatible
    // with STRSXP: [type=NULL]" and crashes the process; limma tolerates an
// unnamed list, so guard the conversion and simply omit the names.
    Rcpp::CharacterVector set_names;
    Rcpp::RObject nm_attr = index.attr("names");
    if (!Rf_isNull(nm_attr))
        set_names = Rcpp::CharacterVector(nm_attr);
    if (set_names.size() == nsets) {
        NGenes.attr("names") = set_names;
        tstat.attr("names") = set_names;
        PValue.attr("names") = set_names;
    }

    return Rcpp::List::create(Rcpp::Named("NGenes") = NGenes,
                              Rcpp::Named("t") = tstat,
                              Rcpp::Named("PValue") = PValue);
}

//' Rotation gene-set test (roast)
//'
//' Native counterpart of \code{limma::mroast}: the set's t-statistics are
//' compared against a null built by rotating the residual coordinates, so no
//' permutation of sample labels is required and the design structure is kept.
//'
//' The normal draws are consumed in exactly the order limma uses
//' (\code{matrix(rnorm(nroti * neffects), nroti, neffects)}), so calling this
//' immediately after \code{set.seed()} reproduces \code{limma::mroast}
//' bit-for-bit.  Do not introduce an extra \code{set.seed()} here.
//'
//' @param E G x n matrix of log-expression values.
//' @param weights optional G x n matrix of positive precision weights.
//' @param design n x p design matrix.
//' @param contrast_idx 1-based index of the contrast column of \code{design}.
//' @param index list of integer gene-set memberships (1-based row numbers).
//' @param set_statistic one of \code{"mean"}, \code{"floormean"},
//'   \code{"mean50"}, \code{"msq"}.
//' @param gene_weights optional numeric vector of length G.
//' @param nrot number of rotations.
//' @return a list with \code{NGenes}, \code{Active} (a 4-column matrix of
//'   active proportions), \code{PV} (a 4-column matrix of p-values for
//'   down / up / upordown / mixed) plus \code{var.prior} and \code{df.prior}.
//' @noRd
// [[Rcpp::export]]
Rcpp::List st_roast_cpp(Rcpp::NumericMatrix E,
                        Rcpp::Nullable<Rcpp::NumericMatrix> weights,
                        Rcpp::NumericMatrix design, int contrast_idx,
                        Rcpp::List index, std::string set_statistic,
                        Rcpp::Nullable<Rcpp::NumericVector> gene_weights,
                        int nrot) {
    st_check_common(E, design, contrast_idx, index);
    if (nrot < 1)
        Rcpp::stop("nrot must be positive");
    if (set_statistic != "mean" && set_statistic != "floormean" &&
        set_statistic != "mean50" && set_statistic != "msq")
        Rcpp::stop("set.statistic must be one of mean, floormean, mean50, msq");

    const int G = E.nrow();
    const int n = E.ncol();
    const int p = design.ncol();
    const int dfr = n - p;
    const int nres = dfr;
    const int neff = nres + 1;

    Rcpp::NumericMatrix Eff = st_effects_cpp(E, weights, design, contrast_idx);

    // limma::mroast: s2 <- rowMeans(Effects[, -1]^2); sv <- squeezeVar(s2, dfr)
    std::vector<double> s2(G);
    for (int g = 0; g < G; ++g) {
        double s = 0.0;
        for (int k = 1; k <= nres; ++k)
            s += Eff(g, k) * Eff(g, k);
        s2[g] = s / nres;
    }
    double var_prior = NA_REAL, df_prior = NA_REAL;
    st_fit_fdist(s2, static_cast<double>(dfr), var_prior, df_prior);

    std::vector<double> var_post(G);
    if (R_FINITE(df_prior)) {
        for (int g = 0; g < G; ++g)
            var_post[g] = (dfr * s2[g] + df_prior * var_prior) / (dfr + df_prior);
    } else {
        for (int g = 0; g < G; ++g)
            var_post[g] = var_prior;
    }
    const double df_total = df_prior + static_cast<double>(dfr);

    const bool has_gw = gene_weights.isNotNull();
    std::vector<double> GW;
    if (has_gw) {
        Rcpp::NumericVector gw(gene_weights);
        if (gw.size() != G)
            Rcpp::stop("gene.weights vector should be of length nrow(y)");
        GW.assign(gw.begin(), gw.end());
    }

    const int nsets = index.size();
    Rcpp::NumericMatrix PV(nsets, 4), ACT(nsets, 4);
    Rcpp::NumericVector NGenes(nsets);
    Rcpp::CharacterVector pv_names = Rcpp::CharacterVector::create("Down", "Up", "UpOrDown", "Mixed");
    Rcpp::colnames(PV) = pv_names;
    Rcpp::colnames(ACT) = pv_names;

    const double chimed = R::qnorm(0.25, 0.0, 1.0, 0, 0);  // upper quartile of N(0,1)
    std::vector<int> rows;
    std::vector<double> modt, mtr, zt, buf;

    for (int si = 0; si < nsets; ++si) {
        st_resolve_set(index, si, G, rows);
        const int m = static_cast<int>(rows.size());
        NGenes[si] = m;

        // ---- observed statistics -------------------------------------------
        // limma: modt <- effects[, 1]/sqrt(var.post); zscoreT(..., "bailey")
        const double df_tw = std::min(df_total, 10000.0);
        modt.assign(m, 0.0);
        for (int j = 0; j < m; ++j)
            modt[j] = st_zscore_t_bailey(Eff(rows[j], 0) / std::sqrt(var_post[rows[j]]), df_tw);

        double a1, a2;
        if (!has_gw) {
            int c1 = 0, c2 = 0;
            for (int j = 0; j < m; ++j) {
                if (modt[j] > M_SQRT2)
                    ++c1;
                if (modt[j] < -M_SQRT2)
                    ++c2;
            }
            a1 = static_cast<double>(c1) / m;
            a2 = static_cast<double>(c2) / m;
        } else {
            // limma: s <- sign(gene.weights); ss <- sum(abs(s))
            //        a1 <- sum(s * modt > sqrt(2))/ss
            double ss = 0.0;
            for (int j = 0; j < m; ++j)
                ss += std::fabs((GW[rows[j]] > 0.0) - (GW[rows[j]] < 0.0));
            double n1 = 0.0, n2 = 0.0;
            for (int j = 0; j < m; ++j) {
                const double sg = (GW[rows[j]] > 0.0) - (GW[rows[j]] < 0.0);
                if (sg * modt[j] > M_SQRT2)
                    n1 += 1.0;
                if (sg * modt[j] < -M_SQRT2)
                    n2 += 1.0;
            }
            a1 = ss > 0.0 ? n1 / ss : NA_REAL;
            a2 = ss > 0.0 ? n2 / ss : NA_REAL;
        }

        int half1, half2;
        if (m % 2 == 0) {
            half1 = m / 2;
            half2 = half1 + 1;
        } else {
            half1 = half2 = m / 2 + 1;
        }

        std::vector<double> statobs(4, 0.0);  // down, up, upordown, mixed
        {
            // gene weights multiply modt for mean / mean50 and abs() versions
            // for floormean / msq, exactly as .roastEffects does.
            std::vector<double> wm(modt);
            if (has_gw)
                for (int j = 0; j < m; ++j)
                    wm[j] = GW[rows[j]] * modt[j];

            if (set_statistic == "mean") {
                double mn = 0.0, s = 0.0;
                for (int j = 0; j < m; ++j) {
                    mn += wm[j];
                    s += std::fabs(wm[j]);
                }
                statobs[0] = -mn / m;
                statobs[1] = mn / m;
                statobs[3] = s / m;
            } else if (set_statistic == "floormean") {
                double d = 0.0, u = 0.0, s = 0.0;
                for (int j = 0; j < m; ++j) {
                    const double a = std::max(std::fabs(modt[j]), chimed);
                    const double av = has_gw ? GW[rows[j]] * a : a;
                    s += av;
                    d += std::max(-wm[j], 0.0);
                    u += std::max(wm[j], 0.0);
                }
                statobs[0] = d / m;
                statobs[1] = u / m;
                statobs[2] = std::max(statobs[0], statobs[1]);
                statobs[3] = s / m;
            } else if (set_statistic == "mean50") {
                buf = wm;
                std::sort(buf.begin(), buf.end());
                double d = 0.0, u = 0.0;
                for (int j = 0; j < half1; ++j)
                    d += buf[j];
                for (int j = half2 - 1; j < m; ++j)
                    u += buf[j];
                statobs[0] = -d / half1;
                statobs[1] = u / (m - half2 + 1);
                statobs[2] = std::max(statobs[0], statobs[1]);
                buf.assign(m, 0.0);
                for (int j = 0; j < m; ++j)
                    buf[j] = std::fabs(wm[j]);
                std::sort(buf.begin(), buf.end());
                double s = 0.0;
                for (int j = half2 - 1; j < m; ++j)
                    s += buf[j];
                statobs[3] = s / (m - half2 + 1);
            } else {  // msq
                // limma: modt2 <- modt^2
                //        if (gw) { modt2 <- abs(gw) * modt2; modt <- gw * modt }
                //        statobs["down"] <- sum(modt2[modt < 0])/nset  (signed modt!)
                double d = 0.0, u = 0.0, s = 0.0;
                for (int j = 0; j < m; ++j) {
                    double vv = modt[j] * modt[j];
                    if (has_gw)
                        vv *= std::fabs(GW[rows[j]]);
                    s += vv;
                    if (wm[j] < 0.0)
                        d += vv;
                    if (wm[j] > 0.0)
                        u += vv;
                }
                statobs[0] = d / m;
                statobs[1] = u / m;
                statobs[2] = std::max(statobs[0], statobs[1]);
                statobs[3] = s / m;
            }
        }

        // ---- rotations ------------------------------------------------------
        // limma draws the whole block at once:
        //     modtr <- matrix(rnorm(nroti * neffects), nroti, neffects)
        //     modtr <- modtr / sqrt(rowSums(modtr^2))
        // R's matrix() fills COLUMN-MAJOR, so draw number d (1-based) lands at
        // rotation ((d-1) %% nroti) + 1 and effect coordinate ((d-1) %/% nroti) + 1.
        // Drawing neffects values per rotation instead (row-major) permutes the
        // stream and moves roast's p-values by several counts out of nrot, which
        // is exactly the bit-exactness this method is required to preserve.
        std::vector<int> count(4, 0);
        std::vector<double> z(neff), mt(m), rsq(m), s2r(m);
        mtr.assign(m, 0.0);
        zt.assign(m, 0.0);

        std::vector<double> block(static_cast<std::size_t>(nrot) * neff);
        for (std::size_t d = 0; d < block.size(); ++d)
            block[d] = R::rnorm(0.0, 1.0);

        for (int r = 0; r < nrot; ++r) {
            double nz = 0.0;
            for (int k = 0; k < neff; ++k) {
                // column-major lookup: row r, column k
                z[k] = block[static_cast<std::size_t>(k) * nrot + r];
                nz += z[k] * z[k];
            }
            nz = std::sqrt(nz);
            for (int k = 0; k < neff; ++k)
                z[k] /= nz;

            // modtr <- tcrossprod(effects, modtr) -- the set's row of the
            // rotated effects.
            for (int j = 0; j < m; ++j) {
                double s = 0.0, rs = 0.0;
                for (int k = 0; k < neff; ++k) {
                    const double v = Eff(rows[j], k);
                    s += v * z[k];
                    rs += v * v;
                }
                mt[j] = s;
                rsq[j] = rs;
            }

            // limma: s2r <- (rowSums(effects^2) - modtr^2)/df.residual
            //        s2r <- (df.prior*var.prior + df.residual*s2r)/df.total
            for (int j = 0; j < m; ++j) {
                const double raw = (rsq[j] - mt[j] * mt[j]) / dfr;
                s2r[j] = R_FINITE(df_prior)
                             ? (df_prior * var_prior + dfr * raw) / df_total
                             : var_prior;
            }
            for (int j = 0; j < m; ++j)
                mtr[j] = mt[j] / std::sqrt(s2r[j]);
            for (int j = 0; j < m; ++j)
                zt[j] = st_zscore_t_bailey(mtr[j], df_tw);

            std::vector<double> wr(zt);
            if (has_gw)
                for (int j = 0; j < m; ++j)
                    wr[j] = GW[rows[j]] * zt[j];

            double sd_ = 0.0, su_ = 0.0, smix = 0.0, sup = 0.0;
            if (set_statistic == "mean") {
                double mn = 0.0, s = 0.0;
                for (int j = 0; j < m; ++j) {
                    mn += wr[j];
                    s += std::fabs(wr[j]);
                }
                sd_ = -mn / m;
                su_ = mn / m;
                smix = s / m;
            } else if (set_statistic == "floormean") {
                double d = 0.0, u = 0.0, s = 0.0;
                for (int j = 0; j < m; ++j) {
                    const double a = std::max(std::fabs(zt[j]), chimed);
                    const double av = has_gw ? GW[rows[j]] * a : a;
                    s += av;
                    d += std::max(-wr[j], 0.0);
                    u += std::max(wr[j], 0.0);
                }
                sd_ = d / m;
                su_ = u / m;
                smix = s / m;
            } else if (set_statistic == "mean50") {
                buf = wr;
                std::sort(buf.begin(), buf.end());
                double d = 0.0, u = 0.0;
                for (int j = 0; j < half1; ++j)
                    d += buf[j];
                for (int j = half2 - 1; j < m; ++j)
                    u += buf[j];
                sd_ = -d / half1;
                su_ = u / (m - half2 + 1);
                std::vector<double> av(m);
                for (int j = 0; j < m; ++j)
                    av[j] = std::fabs(wr[j]);
                std::sort(av.begin(), av.end());
                double s = 0.0;
                for (int j = half2 - 1; j < m; ++j)
                    s += av[j];
                smix = s / (m - half2 + 1);
            } else {  // msq
                // limma rotation branch: gw <- sqrt(abs(gw)); modtr <- gw * modtr
                // then statrot[,"down"] <- colMeans(pmax(-modtr,0)^2) etc.
                // The squared statistic therefore carries |gw| exactly once,
                // and the down/up split uses the SIGN OF THE SCALED modtr.
                double d = 0.0, u = 0.0, s = 0.0;
                for (int j = 0; j < m; ++j) {
                    double scaled = zt[j];
                    if (has_gw)
                        scaled *= std::sqrt(std::fabs(GW[rows[j]]));
                    d += scaled < 0.0 ? scaled * scaled : 0.0;
                    u += scaled > 0.0 ? scaled * scaled : 0.0;
                    s += scaled * scaled;
                }
                sd_ = d / m;
                su_ = u / m;
                smix = s / m;
            }
            sup = std::max(sd_, su_);

            // limma counts every rotation statistic against BOTH the down and
            // the up observed statistic: sum(statrot[, c("down","up")] > obs)
            if (sd_ > statobs[0])
                ++count[0];
            if (su_ > statobs[0])
                ++count[0];
            if (sd_ > statobs[1])
                ++count[1];
            if (su_ > statobs[1])
                ++count[1];
            // For "mean" the upordown count is NOT accumulated here; limma
            // assigns it min(count[down], count[up]) after the loop.
            if (set_statistic != "mean" && sup > statobs[2])
                ++count[2];
            if (smix > statobs[3])
                ++count[3];
        }

        if (set_statistic == "mean")
            count[2] = std::min(count[0], count[1]);

        PV(si, 0) = (count[0] + 1.0) / (2.0 * nrot + 1.0);
        PV(si, 1) = (count[1] + 1.0) / (2.0 * nrot + 1.0);
        PV(si, 2) = (count[2] + 1.0) / (1.0 * nrot + 1.0);
        PV(si, 3) = (count[3] + 1.0) / (1.0 * nrot + 1.0);

        ACT(si, 0) = a2;  // Active.Prop "Down"
        ACT(si, 1) = a1;  // Active.Prop "Up"
        ACT(si, 2) = std::max(a1, a2);
        ACT(si, 3) = a1 + a2;
    }

    // index.names() is NULL when the caller passes an unnamed gene-set list.
    // Assigning that straight into a CharacterVector raises "Not compatible
    // with STRSXP: [type=NULL]" and crashes the process; limma tolerates an
// unnamed list, so guard the conversion and simply omit the names.
    Rcpp::CharacterVector set_names;
    Rcpp::RObject nm_attr = index.attr("names");
    if (!Rf_isNull(nm_attr))
        set_names = Rcpp::CharacterVector(nm_attr);
    if (set_names.size() == nsets) {
        Rcpp::rownames(PV) = set_names;
        Rcpp::rownames(ACT) = set_names;
        NGenes.attr("names") = set_names;
    }

    return Rcpp::List::create(Rcpp::Named("NGenes") = NGenes,
                              Rcpp::Named("Active") = ACT,
                              Rcpp::Named("PV") = PV,
                              Rcpp::Named("var.prior") = var_prior,
                              Rcpp::Named("df.prior") = df_prior);
}
