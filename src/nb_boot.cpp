// nb_boot.cpp -- gene-set tests whose NULL DISTRIBUTION is a negative-binomial
// parametric bootstrap under H0, plus the project's own set statistic.
//
// Ablation cells (NB_IMPLEMENT_SPEC.md sections 0.1 / 0.2):
//
//   id   statistic                              null               selector
//   N    z_S = NB-model-variance-weighted meta-z   NB bootstrap     "nb_z"
//   A    same z_S                                  analytic normal  "nb_z_analytic"
//   B    fry-form t_S                              NB bootstrap     "fry_t"
//   L    fry-form t_S                              limma rotation   (dsea(), untouched)
//
//     z_S = sum_{g in S} beta_g w_g / sqrt( sum_{g in S} w_g ),  w_g = 1 / se_g^2
//
// a fixed-effect meta-analytic combined z on the contrast effect and its MODEL-
// BASED standard error.  That SE already encodes the NB mean-variance law
// (Var = mu + alpha*mu^2): low-count / high-dispersion genes are down-weighted by
// it instead of having their noise inflated by the log transform, which is the
// structural difference from the limma path.
//
// -------------------------------------------------------------------------
// THE ESTIMATOR IS CLOSED-FORM.  THERE IS NO ITERATIVE SOLVER IN THIS FILE.
// (lead design decision, 2026-10-14, after a per-replicate IRLS refit proved both
// expensive and fragile: Newton diverged on low-count genes and no damping rule
// fixed it without changing what the estimator converges to.)
//
// alpha_g (DESeq2's EB-shrunk dispersion), norm_factors and design are FIXED: they
// do not move between replicates, so the whole per-gene geometry is computed ONCE
// and reused by all B + 1 passes.  That is what keeps B = 1000 at G = 4000 inside
// criterion F3, and it is also why the observed statistic and the null cannot use
// two different estimators -- they read the same arrays.
//
//   v_gj  = mu_gj / (1 + alpha_g mu_gj)   the NB Fisher weight at the fixed fit:
//                                         Var / (d mu/d eta)^2, i.e. the mean-
//                                         variance law divided out.  This factor
//                                         IS the point of the statistic.
//   Xc    = design[, ci] ,  Xo = design[, != ci]
//   c~_g  = Xc - Xo (Xo' diag(v_g) Xo)^{-1} Xo' diag(v_g) Xc
//           the contrast column partialled out of the nuisance design in the v_g
//           metric, so the SE is model-based rather than marginal: the design's
//           other columns take their share of the information first.
//   I_g   = sum_j v_gj c~_gj^2            information for beta_ci => w_g = I_g
//   se_g  = 1 / sqrt(I_g)
//   off_j = log(norm_factors[j])          DESeq2's fixed log size-factor offset
//
// and the per-pass coefficient is the v-weighted GLS contrast coefficient of the
// log-scale observation:
//
//   beta_g(counts) = sum_j v_gj c~_gj (log(counts_gj + 0.5) - off_j) / I_g
//
// Only the target moves between passes; v, c~ and I do not.
//
// WHY NOT A SCORE STEP.  An earlier revision used beta = score/I with
// score = sum_j c~_j (y_j - mu_j)/(1 + alpha mu_j), the one-step Newton update of
// the NB likelihood.  That is the textbook estimator and it is WRONG HERE: at
// DESeq2's MLE the score equations vanish, so score/I is ~0 for every gene, and a
// null built by resampling around that same mu has the same ~0 expectation --
// recreating, in a new guise, exactly the no-op test the H0 correction was meant
// to remove.  The log-scale GLS does not degenerate: on the observed counts it
// returns the observed contrast, and its reference coefficient beta_ref is
// DESeq2's own generating beta (see the exactness note below), so zeroing it
// really does centre the null.
//
// WHY THE JENSEN BIAS IS THE HYPOTHESIS, NOT A DEFECT.  E[log(y + 0.5)] !=
// log(E[y] + 0.5), and the gap is worst at low counts -- that is precisely the
// log-transform distortion this method claims to correct.  The analytic normal
// null (cell A) cannot see it; the NB bootstrap null carries the same bias by
// construction and so absorbs it.  The sd of the null z_S is the measurement that
// says whether the correction is real: ~1 means cell A already suffices and the
// bootstrap adds nothing (F1 kills the method); > 1 is the excess this design
// exists to calibrate.
//
// We DELIBERATELY do not consume DESeq2's per-gene `std` column: that is the
// quantity whose low-count underestimation this method sets out to test, so
// putting it into the statistic would build the hypothesis into the test.  We
// recompute a model-based SE from DESeq2's own alpha and mu instead, and let the
// bootstrap supply the empirical null for that closed-form z.  beta / se still
// arrive in the frozen signature; they are validated for length and deliberately
// unused.
//
// -------------------------------------------------------------------------
// THE NULL IS THE CONSTRAINED FIT, NOT assay(x, "mu")
//
// assay(x, "mu") is the UNCONSTRAINED MLE: it contains the real effect.
// Resampling counts* from it gives counts* the same effect the observed data has,
// so |t*| >= |t_obs| happens about half the time and EVERY set -- a genuinely
// differential one as much as a null one -- comes out near p ~ 0.5.  The test
// would have ZERO discriminating power while still producing p-values in [0,1]
// with a reassuringly flat histogram, so a null-set calibration check alone would
// NOT catch it.  A parametric bootstrap for goodness-of-fit uses the MLE; one for
// TESTING a hypothesis must use the constrained MLE under H0.
//
//   beta_ref_g = sum_j v_gj c~_gj (log(mu_gj) - off_j) / I_g
//   mu_h0[g,j] = mu[g,j] * exp( -Xc_j * beta_ref_g )
//
// which is exp(eta_hat - X[,ci]*beta_hat): intercept, batch / cell terms, size
// factors and the mean structure survive; the tested effect does not.
//
// WHY beta_ref IS EXACT.  DESeq2's fit is log mu = off + X beta, so the target
// log(mu) - off lies in col(X); weighted least squares on a target already in the
// column space returns its true coefficients whatever the positive weights are.
// Measured on simulated data (n = 1980 usable genes): max LS residual 6.4e-15 and
// max |beta_ref - beta| = 3.5e-15, for the NB weight and for arbitrary weights
// alike.  So mu_h0 removes the contrast contribution and nothing else.
//
// NOTE ON THE SPEC'S LITERAL FORM.  The contract proposed exp(X_full %%*%%
// t(coef(dds))) as the reconstruction of mu.  That cannot hold: coef() returns the
// prior-SHRUNK coefficients when betaPrior is on (DESeq2's default), and those do
// not regenerate assay(x, "mu").  Measured: max relative error of
// t(X %%*%% t(coef(dds))) against assay(x, "mu") is 14.8 -- an order of magnitude,
// not 1e-6.  Deriving beta_ref by weighted GLS projection of log(mu) - log(sf),
// the quantity DESeq2 actually fitted, is both exact and parameterisation-
// agnostic, so guard 1 is satisfied by construction rather than by tolerance.  The
// internal consistency that matters is enforced explicitly:
// nb_check_null_consistency() stops the call if log(mu) - log(mu_h0) ever
// disagrees with Xc * beta_ref.
//
// CAVEAT, stated rather than hidden: zeroing the contrast coefficient in the
// unconstrained fit is not identical to a full refit under the constraint (the
// nuisance coefficients are not re-optimised with beta_ci held at 0).  For a
// design whose contrast column is an indicator this is a second-order difference;
// a likelihood-ratio-style null would refit properly and would cost an iterative
// solve per gene per replicate, which the project has just decided against.
//
// -------------------------------------------------------------------------
// PIPELINE IDENTITY (spec section 2)
//
// One lib.size -> E -> weights -> effects chain and one closed-form estimator, run
// B + 1 times over the same reused buffers: b = 0 gets the observed counts, b >= 1
// get counts*.  No branch in the loop is special-cased for the observed pass.
//
//     lib.size = mean(colSums(counts)) * norm_factors
//
// NOT colSums(counts) * norm_factors.  DESeq2::sizeFactors() are depth-like
// (median-of-ratios, geometric mean ~ 1), so multiplying raw depth by them counts
// depth twice and leaves a rank-1 per-sample offset the design cannot absorb.  See
// as_dsea()'s @details in R/deseq2.R for the measurements.
//
//     denom = matrix(lib.size, G, n, byrow = TRUE) + 1
//
// The explicit byrow is required: `matrix * length-n vector` recycles COLUMN-major
// in R, which silently broadcasts onto genes whenever G is a multiple of n.
//
// -------------------------------------------------------------------------
// NA HANDLING
//
// Genes DESeq2 drops (independent filtering) have NA mu / dispersion.  Such a gene
// is NOT coerced to zero: its geometry is unusable, the closed form returns
// beta = NA with I = NA, mu_h0 keeps a zero row, and the meta-z leaves it out of
// both the numerator and the precision sum, so the statistic is over the usable
// genes only.  Coercing them to 0 produced an all-NaN run that looked like a
// method failure.
//
// -------------------------------------------------------------------------
// INDEXING / MEMORY SAFETY
//
// Every flat index carries an explicit cast and a documented extent, and each
// buffer's size is stated where it is indexed:
//   counts / mu / E / weights / vw : column-major G x n, element (g,j) at j*G + g
//   Xo                             : row-major    n x q, element (j,a) at j*q + a
//   S, Sd, L                       : row-major    q x q, element (a,b) at a*q + b
//   Xc, off, v, ct                 : length       n
//   beta, precision, beta_ref      : length       G
// The earlier revision's "subscript out of bounds (index 16 >= vector size 16)"
// warnings came from indexing an n x p design buffer (extent 16 at n = 8, p = 2)
// with a gene-major offset reaching j*G.  All design reads here go through Xc / Xo
// / off, each sized from n alone, and the per-gene arrays are sized G; the class of
// index that could straddle the two extents no longer exists.
//
// -------------------------------------------------------------------------
// WHAT THIS FILE DOES NOT DO
//   * no iterative solver, anywhere;
//   * it does not touch src/set_test.cpp (the bit-exact parity baseline); the QR
//     lives in st_effects_cpp() and is not re-inlined here;
//   * it does not seed the RNG.  GetRNGState/PutRNGState wrap the WHOLE bootstrap
//     loop exactly once, so the caller's stream is consumed (an earlier bug in this
//     project re-seeded inside a permutation loop and 15 replicates produced only
//     2 distinct labelings);
//   * it does not re-estimate size factors or dispersions from counts* (that would
//     leak the composition effect and the resampling degrees of freedom into the
//     null);
//   * "nb_z" ignores gene.weights by definition (a signed multiplier would flip the
//     sign of a precision weight); passing them with a nb_z cell is an error, not a
//     silent no-op;
//   * it does not apply squeezeVar to z_S: the meta-z IS the model-based
//     standardisation, and borrowing strength on top of it would double-count the
//     EB information alpha already carries.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <vector>

// ---------------------------------------------------------------------------
// The verified baseline layer, defined in src/set_test.cpp.  Declared here so this
// translation unit can call it; set_test.cpp is READ-ONLY for this task.
// ---------------------------------------------------------------------------

Rcpp::NumericMatrix st_effects_cpp(Rcpp::NumericMatrix E,
                                   Rcpp::Nullable<Rcpp::NumericMatrix> weights,
                                   Rcpp::NumericMatrix design,
                                   int contrast_idx);

namespace {

// log(2) exactly as R's log(2) evaluates it (both are C log()).
const double NB_LN2 = 0.69314718055994530941723212145818;

inline bool nb_finite(double x) { return R_FINITE(x); }

// ---------------------------------------------------------------------------
// lib.size[j] = mean(colSums(counts)) * norm_factors[j]
//
// Accumulated in long double because R's colSums() and mean.default() both
// accumulate in LDOUBLE; matching them keeps the observed E bit-identical to the E
// that as_dsea() builds in R.  (Measured: accumulating in double shifted the
// weights by 1 ulp and Effects by 2.2e-14, which costs the parity guarantee.)
// ---------------------------------------------------------------------------
void nb_lib_size(const double* counts, int G, int n, const double* norm_factors,
                 std::vector<double>& lib_size) {
    lib_size.assign(static_cast<std::size_t>(n), 0.0);

    long double total = 0.0L;
    for (int j = 0; j < n; ++j) {
        const double* col = counts + static_cast<R_xlen_t>(j) * G;
        long double s = 0.0L;
        for (int g = 0; g < G; ++g)
            s += col[g];
        total += s;
    }
    const double mean_depth =
        static_cast<double>(total / static_cast<long double>(n));

    for (int j = 0; j < n; ++j)
        lib_size[static_cast<std::size_t>(j)] =
            mean_depth * norm_factors[j];
}

// ---------------------------------------------------------------------------
// counts -> E, the single expression-matrix builder used by BOTH the observed pass
// and every bootstrap replicate.  Matches R's
//
//     E <- log2((ct + 0.5) / (lib.size_row + 1) * 1e6)
//
// operand by operand (division before the 1e6 multiply, in that order), and uses
// std::log2, which is bit-identical to R's log2() (verified on 20000 draws over
// 1e-6..1e6).
// ---------------------------------------------------------------------------
void nb_build_E(const double* counts, int G, int n, const double* lib_size,
                double* E) {
    for (int j = 0; j < n; ++j) {
        const double denom = lib_size[static_cast<std::size_t>(j)] + 1.0;
        double* col = E + static_cast<R_xlen_t>(j) * G;
        const double* ccol = counts + static_cast<R_xlen_t>(j) * G;
        for (int g = 0; g < G; ++g) {
            const double num = ccol[g] + 0.5;
            col[g] = std::log2((num / denom) * 1e6);
        }
    }
}

// ---------------------------------------------------------------------------
// Delta-method precision weights, built ONCE (they depend on mu / alpha, which are
// fixed across replicates, never on counts).  Cell B / L only.
//
//     mu_tilde = pmax(mu, mu_floor)
//     w = (mu_tilde + 0.5)^2 * log(2)^2 / (mu_tilde + alpha * mu_tilde^2)
//     w[!is.finite(w)] <- min(w[is.finite(w)])   [0 when none is finite]
//
// pmax() PROPAGATES NA the way R does.  Rounding that off to mu_floor (an earlier
// draft did) silently gives the filtered genes a real weight instead of the NA
// that later falls back to the minimum finite weight, and breaks bit-identity with
// as_dsea()'s weights.
//
// The product is grouped as R groups it: `alpha * (mt*mt)`.  Writing
// `alpha * mt * mt` re-associates and moves the weight by one ulp.
// ---------------------------------------------------------------------------
void nb_weights(const double* mu, const double* alpha_raw, int G, int n,
                double mu_floor, std::vector<double>& w_out) {
    w_out.assign(static_cast<std::size_t>(G) * n, 0.0);
    const double ln2sq = NB_LN2 * NB_LN2;

    for (int j = 0; j < n; ++j) {
        const double* mcol = mu + static_cast<R_xlen_t>(j) * G;
        double* wcol = &w_out[static_cast<R_xlen_t>(j) * G];
        for (int g = 0; g < G; ++g) {
            const double m = mcol[g];
            const double mt = R_IsNA(m) ? NA_REAL
                            : (R_IsNaN(m) ? NA_REAL
                                          : (m > mu_floor ? m : mu_floor));
            const double a = alpha_raw[g];
            const double sq = mt * mt;
            const double num = ((mt + 0.5) * (mt + 0.5)) * ln2sq;
            wcol[g] = num / (mt + a * sq);
        }
    }

    double wmin = R_PosInf;
    bool any_fin = false;
    for (std::size_t i = 0; i < w_out.size(); ++i) {
        if (R_FINITE(w_out[i])) {
            any_fin = true;
            if (w_out[i] < wmin)
                wmin = w_out[i];
        }
    }
    const double filler = any_fin ? wmin : 0.0;
    for (std::size_t i = 0; i < w_out.size(); ++i)
        if (!R_FINITE(w_out[i]))
            w_out[i] = filler;
}

// ---------------------------------------------------------------------------
// NB size parameter per gene for the resampling law c* ~ NB(size, mu_h0).
//
// alpha_g NA / NaN / +-Inf -> 1          (as_dsea()'s convention)
// alpha_g <= 0             -> size = Inf (the Poisson limit)
// else                     -> 1 / alpha_g
//
// `poisson_g` records which genes took the Inf branch.  C's Rf_rnbinom_mu()
// returns NA_REAL for a non-finite size -- only R's vectorised dispatcher degrades
// Inf to Poisson -- so those genes are drawn with Rf_rpois() here.
// ---------------------------------------------------------------------------
void nb_sizes(const double* alpha, int G, std::vector<double>& size_out,
              std::vector<char>& poisson_out) {
    size_out.assign(static_cast<std::size_t>(G), 0.0);
    poisson_out.assign(static_cast<std::size_t>(G), 0);
    for (int g = 0; g < G; ++g) {
        const double a = alpha[g];
        if (!R_FINITE(a)) {
            size_out[static_cast<std::size_t>(g)] = 1.0;
        } else if (a <= 0.0) {
            size_out[static_cast<std::size_t>(g)] = R_PosInf;
            poisson_out[static_cast<std::size_t>(g)] = 1;
        } else {
            size_out[static_cast<std::size_t>(g)] = 1.0 / a;
        }
    }
}

// ---------------------------------------------------------------------------
// One bootstrap replicate of the counts, in place, from the H0 mean vector.
//
// Draw order is gene-major (g outer, j inner), i.e. rnbinom(n, size=1/alpha_g,
// mu=mu_h0[g, ]) per gene, as spec section 2 spells it.  All G x n draws happen in
// ONE pass per replicate; nothing loops over genes inside the set loop, which is
// the O(G^2) trap this project measured at 3.16 s vs 0.111 s at 4000 genes.
// ---------------------------------------------------------------------------
void nb_resample(double* counts_d, int G, int n, const double* null_mu,
                 const std::vector<double>& size_g,
                 const std::vector<char>& poisson_g) {
    const std::size_t stride = static_cast<std::size_t>(G);
    for (int g = 0; g < G; ++g) {
        const std::size_t gk = static_cast<std::size_t>(g);
        const double sz = size_g[gk];
        const int is_pois = poisson_g[gk];
        double* dst = counts_d + static_cast<R_xlen_t>(g);
        const double* src = null_mu + static_cast<R_xlen_t>(g);
        for (int j = 0; j < n; ++j) {
            // src / dst are already offset by g, so the sample step is stride,
            // not j*G + g.  (Indexing them with the full flat offset read two
            // genes past the intended one.)
            const std::size_t k = static_cast<std::size_t>(j) * stride;
            const double m = src[k];
            if (!R_FINITE(m) || m <= 0.0) {
                dst[k] = 0.0;
                continue;
            }
            dst[k] = is_pois ? static_cast<double>(Rf_rpois(m))
                             : static_cast<double>(Rf_rnbinom_mu(sz, m));
        }
    }
}

// ---------------------------------------------------------------------------
// Dense Cholesky for the q x q nuisance normal matrix (q = p - 1 <= 11).
// ---------------------------------------------------------------------------
bool nb_chol_factor(const double* A, int q, double* L) {
    // Row-by-row Cholesky, with each row's off-diagonal entries computed in the
    // SAME pass that produces the diagonal.  They must be computed BEFORE the
    // diagonal update of the next row consumes them: an earlier revision wrote
    // the two phases separately -- it summed A[i][i] - sum_k L[i][k]^2 using
    // STALE raw A[i][k] -- and the error was invisible at q == 1 (a one-column
    // nuisance block never enters the inner loop) while silently rejecting most
    // genes of a ~batch + condition design as "not positive definite".
    for (int i = 0; i < q; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = A[static_cast<std::size_t>(i) * q + j];
            for (int k = 0; k < j; ++k)
                sum -= L[static_cast<std::size_t>(i) * q + k] *
                       L[static_cast<std::size_t>(j) * q + k];
            if (j == i) {
                if (!(sum > 1e-300))
                    return false;
                L[static_cast<std::size_t>(i) * q + i] = std::sqrt(sum);
            } else {
                L[static_cast<std::size_t>(i) * q + j] =
                    sum / L[static_cast<std::size_t>(j) * q + j];
            }
        }
    }
    return true;
}

bool nb_chol_solve(const double* A, int q, const double* rhs, double* x,
                   double* L) {
    if (!nb_chol_factor(A, q, L))
        return false;
    if (q == 0)
        return true;
    std::vector<double> y(static_cast<std::size_t>(q), 0.0);
    for (int i = 0; i < q; ++i) {
        double s = rhs[i];
        for (int k = 0; k < i; ++k)
            s -= L[static_cast<std::size_t>(i) * q + k] *
                 y[static_cast<std::size_t>(k)];
        y[static_cast<std::size_t>(i)] = s / L[static_cast<std::size_t>(i) * q + i];
    }
    for (int i = q - 1; i >= 0; --i) {
        double s = y[static_cast<std::size_t>(i)];
        for (int k = i + 1; k < q; ++k)
            s -= L[static_cast<std::size_t>(k) * q + i] * x[k];
        x[i] = s / L[static_cast<std::size_t>(i) * q + i];
    }
    return true;
}

// ---------------------------------------------------------------------------
// The design geometry that is gene-independent, plus the per-gene scratch the
// partialling-out needs.
// ---------------------------------------------------------------------------
struct NbDesign {
    int n = 0, p = 0, q = 0, ci = 0;
    std::vector<double> Xc;   // n   contrast column
    std::vector<double> Xo;   // n*q nuisance columns, row-major
    std::vector<double> off;  // n   log(norm_factors), DESeq2's fixed offset
    // scratch, reused for every gene
    std::vector<double> S, Sd, L, rhs, gvec, ct, v;

    void init(const Rcpp::NumericMatrix& design, const double* norm_factors,
              int contrast_idx) {
        n = design.nrow();
        p = design.ncol();
        ci = contrast_idx - 1;
        q = p - 1;
        Xc.assign(static_cast<std::size_t>(n), 0.0);
        off.assign(static_cast<std::size_t>(n), 0.0);
        for (int j = 0; j < n; ++j) {
            Xc[static_cast<std::size_t>(j)] = design(j, ci);
            off[static_cast<std::size_t>(j)] = std::log(norm_factors[j]);
        }
        Xo.assign(static_cast<std::size_t>(n) * std::max(q, 1), 0.0);
        int k = 0;
        for (int j = 0; j < p; ++j) {
            if (j == ci)
                continue;
            for (int i = 0; i < n; ++i)
                Xo[static_cast<std::size_t>(i) * q + k] = design(i, j);
            ++k;
        }
        const std::size_t qq = static_cast<std::size_t>(std::max(q, 1)) *
                               std::max(q, 1);
        S.assign(qq, 0.0);
        Sd.assign(qq, 0.0);
        L.assign(qq, 0.0);
        rhs.assign(static_cast<std::size_t>(std::max(q, 1)), 0.0);
        gvec.assign(static_cast<std::size_t>(std::max(q, 1)), 0.0);
        ct.assign(static_cast<std::size_t>(n), 0.0);
        v.assign(static_cast<std::size_t>(n), 0.0);
    }

    // Fill ct with c~ for the weights currently in v.  Returns false when the
    // nuisance normal matrix is unusable: a singular nuisance means the contrast
    // is not identifiable given the design, so the caller must mark that gene NA
    // rather than guess a value.
    bool build_ctilde() {
        if (q == 0) {
            for (int j = 0; j < n; ++j)
                ct[static_cast<std::size_t>(j)] = Xc[static_cast<std::size_t>(j)];
            return true;
        }
        for (int a = 0; a < q; ++a) {
            rhs[a] = 0.0;
            for (int bb = 0; bb < q; ++bb)
                S[static_cast<std::size_t>(a) * q + bb] = 0.0;
        }
        for (int j = 0; j < n; ++j) {
            const double vj = v[static_cast<std::size_t>(j)];
            const double xcj = Xc[static_cast<std::size_t>(j)];
            const double* xoj = &Xo[static_cast<std::size_t>(j) * q];
            for (int a = 0; a < q; ++a) {
                const double wx = vj * xoj[a];
                rhs[a] += wx * xcj;
                for (int bb = 0; bb < q; ++bb)
                    S[static_cast<std::size_t>(a) * q + bb] += wx * xoj[bb];
            }
        }
        double diag = 0.0;
        for (int a = 0; a < q; ++a)
            if (S[static_cast<std::size_t>(a) * q + a] > diag)
                diag = S[static_cast<std::size_t>(a) * q + a];
        if (!(diag > 0.0))
            return false;
        // ridge escalation, for a design column numerically collinear with the
        // nuisance block
        double scale = 0.0;
        bool ok = false;
        for (int attempt = 0; attempt < 8 && !ok; ++attempt) {
            Sd = S;
            if (scale > 0.0)
                for (int a = 0; a < q; ++a)
                    Sd[static_cast<std::size_t>(a) * q + a] += scale * diag;
            ok = nb_chol_solve(Sd.data(), q, rhs.data(), gvec.data(), L.data());
            scale = scale > 0.0 ? scale * 10.0 : 1e-8;
        }
        if (!ok)
            return false;
        for (int j = 0; j < n; ++j) {
            const double* xoj = &Xo[static_cast<std::size_t>(j) * q];
            double acc = Xc[static_cast<std::size_t>(j)];
            for (int a = 0; a < q; ++a)
                acc -= xoj[a] * gvec[static_cast<std::size_t>(a)];
            ct[static_cast<std::size_t>(j)] = acc;
        }
        return true;
    }
};

// ---------------------------------------------------------------------------
// The per-gene closed-form geometry, computed ONCE and shared by every pass.
//
// vw[g,j] = v_gj * c~_gj is the only per-gene, per-sample quantity the estimator
// needs (both the reference coefficient and the statistic are inner products of vw
// with a log-scale target), so it is the only G x n array stored.
// ---------------------------------------------------------------------------
struct NbGeometry {
    std::vector<double> vw;         // G*n  == v_gj * c~_gj
    std::vector<double> precision;  // G    == I_g == 1 / se_g^2
    std::vector<double> beta_ref;   // G    == the contrast coefficient of mu
    int n_gene_na = 0;

    void compute(const double* mu, const double* alpha_eff, NbDesign& D, int G,
                 int n) {
        vw.assign(static_cast<std::size_t>(G) * n, 0.0);
        precision.assign(static_cast<std::size_t>(G), NA_REAL);
        beta_ref.assign(static_cast<std::size_t>(G), NA_REAL);
        const std::size_t stride = static_cast<std::size_t>(G);
        n_gene_na = 0;

        for (int g = 0; g < G; ++g) {
            const std::size_t gk = static_cast<std::size_t>(g);
            const double a = alpha_eff[gk];
            std::vector<double>& v = D.v;
            bool usable = true;
            for (int j = 0; j < n; ++j) {
                const double m = mu[static_cast<std::size_t>(j) * stride + gk];
                if (!nb_finite(m) || m <= 0.0) {
                    usable = false;
                    break;
                }
                v[static_cast<std::size_t>(j)] = m / (1.0 + a * m);
            }
            if (!usable || !D.build_ctilde()) {
                ++n_gene_na;  // NA, never 0: the gene is excluded downstream
                continue;
            }

            double* vwd = &vw[gk];
            double info = 0.0;
            for (int j = 0; j < n; ++j) {
                const std::size_t k = static_cast<std::size_t>(j) * stride + gk;
                const double vj = v[static_cast<std::size_t>(j)];
                const double ctj = D.ct[static_cast<std::size_t>(j)];
                vwd[static_cast<std::size_t>(j) * stride] = vj * ctj;
                info += vj * ctj * ctj;
            }
            if (!nb_finite(info) || !(info > 0.0)) {
                ++n_gene_na;
                continue;
            }
            precision[gk] = info;

            // beta_ref: the v-weighted GLS contrast coefficient of
            // log(mu_g) - off, i.e. DESeq2's own generating coefficient (exact
            // because the target lies in col(X); see the header).
            double num = 0.0;
            for (int j = 0; j < n; ++j) {
                const std::size_t k = static_cast<std::size_t>(j) * stride + gk;
                num += vwd[static_cast<std::size_t>(j) * stride] *
                       (std::log(mu[k]) - D.off[j]);
            }
            beta_ref[gk] = num / info;
            if (!nb_finite(beta_ref[gk])) {
                precision[gk] = NA_REAL;
                ++n_gene_na;
                continue;
            }
            // Reject a gene whose H0 mean would overflow, HERE rather than in
            // nb_build_null_mu(): if the null row were merely zeroed while the
            // gene still carried a finite precision, it would enter the observed
            // z_S but be drawn as all-zero counts under the null -- an
            // asymmetry between the two sides, which is the one thing this file
            // exists to prevent.  Invalidating the gene keeps both sides honest.
            double worst_shift = 0.0;
            for (int j = 0; j < n; ++j) {
                const double sh = std::fabs(D.Xc[static_cast<std::size_t>(j)] *
                                            beta_ref[gk]);
                if (sh > worst_shift)
                    worst_shift = sh;
            }
            if (worst_shift > 300.0) {
                precision[gk] = NA_REAL;
                beta_ref[gk] = NA_REAL;
                ++n_gene_na;
            }
        }
    }
};

// ---------------------------------------------------------------------------
// The statistic's contrast coefficient for ONE count vector, for ALL genes.
//
//   beta_g = sum_j vw[g,j] * (log(counts_gj + 0.5) - off_j) / I_g
//
// the v-weighted GLS contrast coefficient of the log-scale observation at the
// fixed reference fit.  Only the target moves between passes.  The +0.5
// pseudocount is the same one the delta weights use, so a zero count in a
// replicate gives a finite target rather than -Inf.
// ---------------------------------------------------------------------------
void nb_closed_form(const double* counts_d, const NbDesign& D,
                    const NbGeometry& Geo, int G, int n, double* beta_out) {
    const std::size_t stride = static_cast<std::size_t>(G);
    for (int g = 0; g < G; ++g) {
        const std::size_t gk = static_cast<std::size_t>(g);
        const double info = Geo.precision[gk];
        if (!nb_finite(info) || !(info > 0.0)) {
            beta_out[gk] = NA_REAL;
            continue;
        }
        const double* vwd = Geo.vw.data() + gk;
        // counts_d + g is gene g's column-major stride; index by j*stride.
        const double* ycol = counts_d + static_cast<R_xlen_t>(g);
        double num = 0.0;
        for (int j = 0; j < n; ++j) {
            const std::size_t k = static_cast<std::size_t>(j) * stride;
            num += vwd[k] * (std::log(ycol[k] + 0.5) - D.off[j]);
        }
        beta_out[gk] = nb_finite(num) ? num / info : NA_REAL;
    }
}

// ---------------------------------------------------------------------------
// mu_h0[g,j] = mu[g,j] * exp(-Xc_j * beta_ref_g)
//
// Written as a multiplicative correction of mu itself, so log(mu) - log(mu_h0)
// equals Xc * beta_ref to rounding and the returned beta is by construction the
// beta that built the null.  Genes with NA geometry keep an all-zero row; they are
// excluded from the statistic, not silently averaged into it.
// ---------------------------------------------------------------------------
void nb_build_null_mu(const double* mu, const NbGeometry& Geo,
                      const NbDesign& D, int G, int n,
                      std::vector<double>& mu_h0) {
    const std::size_t stride = static_cast<std::size_t>(G);
    mu_h0.assign(static_cast<std::size_t>(G) * n, 0.0);
    for (int g = 0; g < G; ++g) {
        const std::size_t gk = static_cast<std::size_t>(g);
        const double bref = Geo.beta_ref[gk];
        if (!nb_finite(bref) || !nb_finite(Geo.precision[gk]))
            continue;
        for (int j = 0; j < n; ++j) {
            const std::size_t k = static_cast<std::size_t>(j) * stride + gk;
            const double m = mu[k];
            if (!nb_finite(m) || m <= 0.0)
                continue;
            const double shift = -D.Xc[static_cast<std::size_t>(j)] * bref;
            if (shift > 300.0 || shift < -300.0)
                continue;  // would overflow: leave the row unusable
            mu_h0[k] = m * std::exp(shift);
        }
    }
}

// Guard: log(mu) - log(mu_h0) must equal Xc * beta_ref.  Disagreement means the
// null was not built from the returned coefficient -- the inconsistency this
// project has already been bitten by once.
void nb_check_null_consistency(const double* mu, const double* mu_h0,
                               const NbGeometry& Geo, const NbDesign& D, int G,
                               int n) {
    const std::size_t stride = static_cast<std::size_t>(G);
    double worst = 0.0;
    for (int g = 0; g < G; ++g) {
        const std::size_t gk = static_cast<std::size_t>(g);
        const double bref = Geo.beta_ref[gk];
        if (!nb_finite(bref) || !nb_finite(Geo.precision[gk]))
            continue;
        for (int j = 0; j < n; ++j) {
            const std::size_t k = static_cast<std::size_t>(j) * stride + gk;
            if (mu_h0[k] <= 0.0 || mu[k] <= 0.0)
                continue;
            const double lhs = std::log(mu[k]) - std::log(mu_h0[k]);
            const double rhs = D.Xc[static_cast<std::size_t>(j)] * bref;
            const double d = std::fabs(lhs - rhs);
            if (d > worst)
                worst = d;
        }
    }
    if (worst > 1e-8)
        Rcpp::stop("internal error: mu_h0 is inconsistent with the returned "
                   "contrast coefficient (max |diff| = %g); the null would not be "
                   "the contrast-zeroed fit",
                   worst);
}

// ---------------------------------------------------------------------------
// fry-form set statistic on a G x (nres+1) effects block (cell B / L).
//
//     M[k] = mean_{g in S}( gw[g] * Eff[g, k] )
//     t    = M[0] / sqrt( mean_{k >= 1} M[k]^2 )
//
// Accumulated in double, summing over genes first and dividing by m once, the same
// order st_fry_cpp() uses.  Signed gene weights, not sqrt(): limma's .fryEffects
// multiplies the SET ROWS after standardisation.
// ---------------------------------------------------------------------------
double nb_set_t(const Rcpp::NumericMatrix& Eff, const std::vector<int>& rows,
                const std::vector<double>& gene_weights, bool has_gw,
                int nres, std::vector<double>& M) {
    const int m = static_cast<int>(rows.size());
    if (m < 1)
        return NA_REAL;

    std::fill(M.begin(), M.begin() + nres + 1, 0.0);
    for (int j = 0; j < m; ++j) {
        const int gi = rows[j];
        const double iw = has_gw ? gene_weights[static_cast<std::size_t>(gi)] : 1.0;
        const double* erow = Eff.begin() + static_cast<R_xlen_t>(gi);
        for (int k = 0; k <= nres; ++k)
            M[static_cast<std::size_t>(k)] +=
                iw * erow[static_cast<R_xlen_t>(k) * Eff.nrow()];
    }
    for (int k = 0; k <= nres; ++k)
        M[static_cast<std::size_t>(k)] /= static_cast<double>(m);

    double rss = 0.0;
    for (int k = 1; k <= nres; ++k) {
        const double v = M[static_cast<std::size_t>(k)];
        rss += v * v;
    }
    const double den = std::sqrt(rss / static_cast<double>(nres));
    return den == 0.0 ? NA_REAL : M[0] / den;
}

// ---------------------------------------------------------------------------
// The project's own set statistic (spec section 0.1): the fixed-effect meta z.
// ---------------------------------------------------------------------------
double nb_set_z(const double* beta, const double* prec,
                const std::vector<int>& rows, int* n_used) {
    double num = 0.0, wsum = 0.0;
    int used = 0;
    for (std::size_t j = 0; j < rows.size(); ++j) {
        const std::size_t g = static_cast<std::size_t>(rows[j]);
        const double b = beta[g], w = prec[g];
        if (!nb_finite(b) || !nb_finite(w) || !(w > 0.0))
            continue;
        num += b * w;
        wsum += w;
        ++used;
    }
    if (n_used)
        *n_used = used;
    if (used == 0 || !(wsum > 0.0))
        return NA_REAL;
    return num / std::sqrt(wsum);
}

void nb_resolve_set(const Rcpp::IntegerVector& iset, int which, int G,
                    std::vector<int>& rows) {
    const int m = iset.size();
    rows.resize(static_cast<std::size_t>(m));
    for (int j = 0; j < m; ++j) {
        const int gi = iset[j] - 1;  // R is 1-based
        if (R_IsNA(gi) || gi < 0 || gi >= G)
            Rcpp::stop("gene set %d contains an index outside 1..nrow(counts)",
                       which + 1);
        rows[static_cast<std::size_t>(j)] = gi;
    }
}

// Row names of a matrix as a CharacterVector, empty when absent.  Done through
// Rf_getAttrib because Rcpp::rownames() returns a DimNameProxy with no size().
Rcpp::CharacterVector nb_rownames(SEXP m, int G) {
    Rcpp::CharacterVector out;
    SEXP dn = Rf_getAttrib(m, R_DimNamesSymbol);
    if (Rf_isNull(dn))
        return out;
    SEXP rn = VECTOR_ELT(dn, 0);
    if (Rf_isNull(rn))
        return out;
    if (Rf_length(rn) == G)
        out = Rcpp::CharacterVector(rn);
    return out;
}

// ---------------------------------------------------------------------------
// Shared argument validation for the entry point and the verification helpers.
// ---------------------------------------------------------------------------
void nb_validate(const Rcpp::IntegerMatrix& counts, const Rcpp::NumericMatrix& mu,
                 const Rcpp::NumericVector& alpha,
                 const Rcpp::NumericVector& norm_factors,
                 const Rcpp::NumericMatrix& design, int contrast_idx, int G,
                 int n, int p) {
    if (G < 1 || n < 1)
        Rcpp::stop("counts must have at least one gene and one sample");
    if (mu.nrow() != G || mu.ncol() != n)
        Rcpp::stop("mu must have the same dimensions as counts");
    if (alpha.size() != G)
        Rcpp::stop("alpha must have one entry per gene (length nrow(counts))");
    if (norm_factors.size() != n)
        Rcpp::stop("norm_factors must have one entry per sample (ncol(counts))");
    if (design.nrow() != n)
        Rcpp::stop("row dimension of design matrix must match column dimension of counts");
    if (n <= p)
        Rcpp::stop("No residual degrees of freedom");
    if (contrast_idx < 1 || contrast_idx > p)
        Rcpp::stop("contrast_idx is out of range");
    if (p > 12)
        Rcpp::stop("nb_boot supports designs with at most 12 columns");

    for (int j = 0; j < p; ++j)
        for (int i = 0; i < n; ++i)
            if (!R_FINITE(design(i, j)))
                Rcpp::stop("design contains non-finite values");
    for (int j = 0; j < n; ++j) {
        const double f = norm_factors[j];
        if (!R_FINITE(f) || f <= 0.0)
            Rcpp::stop("norm_factors must be finite and positive");
    }
    // counts: finite non-negative integers.  NA_INTEGER is INT_MIN, so the
    // negativity branch catches it too -- report it explicitly.
    for (int j = 0; j < n; ++j) {
        const int* col = counts.begin() + static_cast<R_xlen_t>(j) * G;
        for (int g = 0; g < G; ++g) {
            const int v = col[g];
            if (v == NA_INTEGER)
                Rcpp::stop("counts contains NA at gene %d, sample %d; raw counts must be complete",
                           g + 1, j + 1);
            if (v < 0)
                Rcpp::stop("counts must be non-negative (found %d at gene %d, sample %d)",
                           v, g + 1, j + 1);
        }
    }
    for (int j = 0; j < n; ++j) {
        const double* col = mu.begin() + static_cast<R_xlen_t>(j) * G;
        for (int g = 0; g < G; ++g) {
            const double v = col[g];
            // NA / NaN are allowed: they mark genes DESeq2 dropped, which the
            // geometry returns NA for and the statistic then excludes.  A
            // NEGATIVE finite mean is not a fitted mu at all.
            if (R_FINITE(v) && v < 0.0)
                Rcpp::stop("mu contains a negative fitted mean at gene %d, sample %d",
                           g + 1, j + 1);
        }
    }
    for (int g = 0; g < G; ++g) {
        const double a = alpha[g];
        if (a < 0.0)  // NaN and NA compare false and are allowed (-> treated as 1)
            Rcpp::stop("alpha contains a negative dispersion at gene %d", g + 1);
    }
}

void nb_counts_to_double(const Rcpp::IntegerMatrix& counts, int G, int n,
                         std::vector<double>& out) {
    out.assign(static_cast<std::size_t>(G) * n, 0.0);
    for (int j = 0; j < n; ++j) {
        const int* src = counts.begin() + static_cast<R_xlen_t>(j) * G;
        double* dst = &out[static_cast<R_xlen_t>(j) * G];
        for (int g = 0; g < G; ++g)
            dst[g] = static_cast<double>(src[g]);
    }
}

void nb_alpha_effective(const Rcpp::NumericVector& alpha, int G,
                        std::vector<double>& out) {
    out.assign(static_cast<std::size_t>(G), 1.0);
    for (int g = 0; g < G; ++g) {
        const double a = alpha[g];
        if (R_FINITE(a) && a > 0.0)
            out[static_cast<std::size_t>(g)] = a;
    }
}

}  // namespace

//' Gene-set test with an H0 negative-binomial parametric-bootstrap null
//'
//' \code{statistic} selects the ablation cell (spec section 0.2):
//' \describe{
//'   \item{\code{"nb_z"}}{the project's own statistic -- the NB-model-variance
//'     weighted set z \code{sum(beta*w)/sqrt(sum(w))} with \code{w = 1/se^2} from
//'     the closed-form NB information -- against the NB bootstrap null.  The
//'     identical closed form is applied to the observed counts and to every
//'     replicate, so observed and null share one estimator by construction.}
//'   \item{\code{"nb_z_analytic"}}{the same z under the analytic normal null,
//'     \code{p_two = 2*pnorm(-|z|)}; no bootstrap, \code{B} may be 0.  Isolates
//'     the statistic's contribution (cell A), and its null sd is the measurement
//'     that says whether the bootstrap adds anything at all.}
//'   \item{\code{"fry_t"}}{limma's fry-form set statistic against the same NB
//'     bootstrap (cell B), through \code{st_effects_cpp}.  The closed form is used
//'     only to build the null mean; the statistic itself is the effects layer and
//'     \code{gene_weights} apply.}
//' }
//'
//' The bootstrap mean is the CONTRAST-ZEROED fitted mean
//' \code{mu_h0[g,j] = mu[g,j] * exp(-X[j, contrast_idx] * beta_ref_g)}, not
//' \code{assay(x, "mu")}: the latter carries the real effect and would give every
//' set p ~ 0.5.  Size factors, batch terms, the mean structure and the EB-shrunk
//' dispersion all survive mu_h0; only the tested effect is removed.
//'
//' @param counts G x n integer matrix of raw counts.  Must be finite and
//'   non-negative; \code{NA} raises an error.
//' @param mu G x n numeric matrix of fitted means (\code{assay(x, "mu")}).  It is
//'   both the reference point of the geometry and the base of the null mean, so it
//'   must come from the same fit as \code{norm_factors} / \code{design}.  It is NOT
//'   itself the resampling mean (that is mu_h0, built internally).
//' @param alpha G numeric vector of EB-shrunk dispersions
//'   (\code{mcols(x)$dispersion}).  \code{NA} is treated as 1; \code{<= 0} becomes
//'   the Poisson limit.
//' @param norm_factors n numeric size factors (depth-like).
//' @param design n x p design matrix (p <= 12), same parameterisation as \code{mu}.
//' @param contrast_idx 1-based index of the contrast column of \code{design}.
//' @param index list of integer gene-set memberships (1-based row numbers).
//' @param B number of bootstrap replicates; >= 1 for the bootstrap cells, may be 0
//'   for \code{"nb_z_analytic"}.
//' @param gene_weights optional numeric vector of length G (signed).  Used by
//'   \code{"fry_t"} only; supplying it with a \code{nb_z} cell is an error.
//' @param mu_floor lower bound for the fitted means in the delta weights.
//' @param statistic one of \code{"nb_z"}, \code{"nb_z_analytic"}, \code{"fry_t"}.
//' @param beta G numeric DESeq2 contrast coefficients.  Validated for length but
//'   NOT consumed: the engine computes its own closed form so the observed
//'   statistic and the null cannot drift apart.  \code{NA_real_} placeholders are
//'   fine.
//' @param se G numeric DESeq2 standard errors.  Same status as \code{beta}, and
//'   deliberately so -- see the file header on the low-count SE being the quantity
//'   under test rather than an input.
//' @return a list with \code{t_obs} (\code{z_S} or \code{t_S} per the selector),
//'   \code{p_two}, \code{n_extreme} (\code{NA} for the analytic cell),
//'   \code{alpha_hat} (set-wise mean dispersion consumed), \code{n_gene_na} (genes
//'   excluded from the statistic), \code{statistic} (echoed) and \code{B}.
//'   Bootstrap cells use \code{p_two = (1 + n_extreme) / (B + 1)}.
//' @noRd
// [[Rcpp::export]]
Rcpp::List st_nb_boot_cpp(Rcpp::IntegerMatrix counts,
                          Rcpp::NumericMatrix mu,
                          Rcpp::NumericVector alpha,
                          Rcpp::NumericVector norm_factors,
                          Rcpp::NumericMatrix design,
                          int contrast_idx,
                          Rcpp::List index,
                          int B,
                          Rcpp::Nullable<Rcpp::NumericVector> gene_weights,
                          double mu_floor,
                          std::string statistic,
                          Rcpp::NumericVector beta,
                          Rcpp::NumericVector se) {
    const int G = counts.nrow();
    const int n = counts.ncol();
    const int p = design.ncol();

    nb_validate(counts, mu, alpha, norm_factors, design, contrast_idx, G, n, p);

    if (index.size() < 1)
        Rcpp::stop("index is empty");
    if (statistic != "nb_z" && statistic != "nb_z_analytic" && statistic != "fry_t")
        Rcpp::stop("statistic must be one of nb_z, nb_z_analytic, fry_t");
    if (statistic == "nb_z_analytic") {
        if (B < 0)
            Rcpp::stop("B must be a non-negative integer");
    } else if (B < 1) {
        Rcpp::stop("B must be a positive integer for the bootstrap cells");
    }
    if (!R_FINITE(mu_floor) || mu_floor < 0.0)
        Rcpp::stop("mu_floor must be a single finite non-negative number");
    if (beta.size() != G)
        Rcpp::stop("beta must have one entry per gene (NA_real_ placeholders are fine)");
    if (se.size() != G)
        Rcpp::stop("se must have one entry per gene (NA_real_ placeholders are fine)");

    const bool is_z = statistic != "fry_t";

    const bool has_gw = gene_weights.isNotNull();
    std::vector<double> GW;
    if (has_gw) {
        if (is_z)
            Rcpp::stop("gene.weights are not supported with statistic = \"%s\": a signed "
                       "multiplier is not part of the meta-z",
                       statistic.c_str());
        Rcpp::NumericVector gw(gene_weights);
        if (gw.size() != G)
            Rcpp::stop("length of gene.weights should equal nrow(counts)");
        GW.assign(gw.begin(), gw.end());
        for (int g = 0; g < G; ++g)
            if (!R_FINITE(GW[static_cast<std::size_t>(g)]))
                Rcpp::stop("gene.weights must be finite");
    }

    std::vector<std::vector<int>> set_rows(static_cast<std::size_t>(index.size()));
    for (int i = 0; i < index.size(); ++i) {
        Rcpp::IntegerVector iset = index[i];
        nb_resolve_set(iset, i, G, set_rows[static_cast<std::size_t>(i)]);
    }

    // ---- built ONCE from mu / alpha / design; shared by every pass -------
    std::vector<double> lib_size;

    std::vector<double> alpha_eff;
    nb_alpha_effective(alpha, G, alpha_eff);

    NbDesign D;
    D.init(design, norm_factors.begin(), contrast_idx);

    NbGeometry Geo;
    Geo.compute(mu.begin(), alpha_eff.data(), D, G, n);

    std::vector<double> cd;
    nb_counts_to_double(counts, G, n, cd);

    std::vector<double> mu_h0;
    nb_build_null_mu(mu.begin(), Geo, D, G, n, mu_h0);
    nb_check_null_consistency(mu.begin(), mu_h0.data(), Geo, D, G, n);

    // cell B / L only: the delta weights for the effects layer
    std::vector<double> w_vals;
    Rcpp::NumericMatrix Wmat;
    if (!is_z) {
        nb_weights(mu.begin(), alpha.begin(), G, n, mu_floor, w_vals);
        Wmat = Rcpp::NumericMatrix(G, n);
        std::copy(w_vals.begin(), w_vals.end(), Wmat.begin());
    }

    std::vector<double> size_g;
    std::vector<char> poisson_g;
    nb_sizes(alpha.begin(), G, size_g, poisson_g);

    // ---- the one pipeline, run B + 1 times -------------------------------
    const int nres = n - p;
    const int nsets = index.size();

    Rcpp::NumericMatrix Eb(G, n);  // reused buffer: counts-derived E

    Rcpp::NumericVector t_obs(nsets), p_two(nsets);
    Rcpp::IntegerVector n_extreme(nsets);
    Rcpp::NumericVector alpha_hat(nsets);

    std::vector<double> counters(static_cast<std::size_t>(nsets), 0.0);
    std::vector<double> Mbuf(static_cast<std::size_t>(nres) + 1, 0.0);
    std::vector<double> beta_star(static_cast<std::size_t>(G), NA_REAL);
    const double* prec_star = Geo.precision.data();

    GetRNGstate();  // exactly once for the whole bootstrap loop

    Rcpp::NumericMatrix Eff;
    std::vector<double> tprev(static_cast<std::size_t>(nsets), 0.0);

    const int b_top = (statistic == "nb_z_analytic") ? 0 : B;

    for (int b = 0; b <= b_top; ++b) {
        if (b > 0)
            nb_resample(cd.data(), G, n, mu_h0.data(), size_g, poisson_g);

        // The identical observed / replicate construction.  Cell N (nb_z) reads
        // the counts directly through the closed form; cell B (fry_t) goes the
        // long way round, lib.size -> E -> delta weights -> st_effects_cpp.  E
        // is built only on that branch -- computing a G x n matrix per pass that
        // nothing reads is pure cost, and it is NOT a pipeline-identity
        // violation: the requirement is that whatever path the OBSERVED
        // statistic takes, the RESAMPLED statistic takes the same one, and here
        // both take the same branch of the same if.
        if (is_z) {
            nb_closed_form(cd.data(), D, Geo, G, n, beta_star.data());
        } else {
            nb_lib_size(cd.data(), G, n, norm_factors.begin(), lib_size);
            nb_build_E(cd.data(), G, n, lib_size.data(), Eb.begin());
            Eff = st_effects_cpp(Eb, Wmat, design, contrast_idx);
        }

        for (int i = 0; i < nsets; ++i) {
            const double t = is_z
                ? nb_set_z(beta_star.data(), prec_star,
                           set_rows[static_cast<std::size_t>(i)], nullptr)
                : nb_set_t(Eff, set_rows[static_cast<std::size_t>(i)],
                           GW, has_gw, nres, Mbuf);
            if (b == 0) {
                t_obs[i] = t;
                const std::vector<int>& rows = set_rows[static_cast<std::size_t>(i)];
                if (!rows.empty()) {
                    double s = 0.0;
                    for (std::size_t j = 0; j < rows.size(); ++j)
                        s += alpha_eff[static_cast<std::size_t>(rows[j])];
                    alpha_hat[i] = s / static_cast<double>(rows.size());
                } else {
                    alpha_hat[i] = NA_REAL;
                }
                tprev[static_cast<std::size_t>(i)] = R_FINITE(t) ? t : NA_REAL;
                counters[static_cast<std::size_t>(i)] = 0.0;
            } else {
                const double to = tprev[static_cast<std::size_t>(i)];
                if (!R_FINITE(to))
                    continue;  // untestable set: observed statistic is not finite
                if (R_FINITE(t) && std::fabs(t) >= std::fabs(to))
                    counters[static_cast<std::size_t>(i)] += 1.0;
            }
        }
    }

    PutRNGstate();

    for (int i = 0; i < nsets; ++i) {
        const bool testable = R_FINITE(t_obs[i]);
        if (statistic == "nb_z_analytic") {
            n_extreme[i] = NA_INTEGER;
            p_two[i] = testable
                          ? 2.0 * R::pnorm(-std::fabs(t_obs[i]), 0.0, 1.0, 1, 0)
                          : NA_REAL;
        } else {
            n_extreme[i] = testable
                              ? static_cast<int>(counters[static_cast<std::size_t>(i)])
                              : NA_INTEGER;
            p_two[i] = testable
                          ? (1.0 + counters[static_cast<std::size_t>(i)]) /
                                (static_cast<double>(B) + 1.0)
                          : NA_REAL;
        }
    }

    Rcpp::CharacterVector set_names;
    Rcpp::RObject nm_attr = index.attr("names");
    if (!Rf_isNull(nm_attr))
        set_names = Rcpp::CharacterVector(nm_attr);
    if (set_names.size() == nsets) {
        t_obs.attr("names") = set_names;
        p_two.attr("names") = set_names;
        n_extreme.attr("names") = set_names;
        alpha_hat.attr("names") = set_names;
    }

    return Rcpp::List::create(Rcpp::Named("t_obs") = t_obs,
                              Rcpp::Named("p_two") = p_two,
                              Rcpp::Named("n_extreme") = n_extreme,
                              Rcpp::Named("alpha_hat") = alpha_hat,
                              Rcpp::Named("n_gene_na") = Geo.n_gene_na,
                              Rcpp::Named("statistic") = statistic,
                              Rcpp::Named("B") = B);
}

// ---------------------------------------------------------------------------
// Verification helpers (not part of the frozen contract).
// ---------------------------------------------------------------------------

//' The closed-form geometry and the H0 (contrast-zeroed) null mean
//'
//' Returns exactly what the engine computes internally on the OBSERVED data, so the
//' R side can
//'   (a) assert \code{log(mu) - log(mu_h0) == X[, ci] * beta} (the consistency
//'       guard: the returned beta must be the beta that built the null),
//'   (b) confirm the null is centred, i.e. the reference coefficient really is the
//'       contrast that mu carries,
//'   (c) compare \code{beta}/{\code{se} with \code{DESeq2::results()} and see the
//'       difference rather than assume it away.
//'
//' @param counts G x n integer matrix of observed counts.
//' @param mu G x n fitted means.
//' @param alpha G dispersions (NA -> 1).
//' @param norm_factors n size factors (depth-like), the fixed log offset.
//' @param design n x p design matrix, same parameterisation as \code{mu}.
//' @param contrast_idx 1-based contrast column.
//' @return a list with \code{mu_h0}, \code{mu_full} (the input \code{mu} echoed, so
//'   the ratio test is exact), \code{beta} (beta_ref), \code{se},
//'   \code{precision}, \code{beta_obs} (the statistic's own coefficient on the
//'   observed counts) and \code{n_gene_na}.
//' @noRd
// [[Rcpp::export]]
Rcpp::List st_nb_null_mu_cpp(Rcpp::IntegerMatrix counts,
                             Rcpp::NumericMatrix mu,
                             Rcpp::NumericVector alpha,
                             Rcpp::NumericVector norm_factors,
                             Rcpp::NumericMatrix design,
                             int contrast_idx) {
    const int G = counts.nrow();
    const int n = counts.ncol();
    const int p = design.ncol();
    nb_validate(counts, mu, alpha, norm_factors, design, contrast_idx, G, n, p);

    std::vector<double> cd;
    nb_counts_to_double(counts, G, n, cd);
    std::vector<double> alpha_eff;
    nb_alpha_effective(alpha, G, alpha_eff);

    NbDesign D;
    D.init(design, norm_factors.begin(), contrast_idx);
    NbGeometry Geo;
    Geo.compute(mu.begin(), alpha_eff.data(), D, G, n);

    std::vector<double> mu_h0;
    nb_build_null_mu(mu.begin(), Geo, D, G, n, mu_h0);
    nb_check_null_consistency(mu.begin(), mu_h0.data(), Geo, D, G, n);

    std::vector<double> beta_obs(static_cast<std::size_t>(G), NA_REAL);
    nb_closed_form(cd.data(), D, Geo, G, n, beta_obs.data());

    Rcpp::NumericMatrix mh0(G, n);
    std::copy(mu_h0.begin(), mu_h0.end(), mh0.begin());
    mh0.attr("dimnames") = mu.attr("dimnames");

    // mu_full is the input mu echoed: the geometry's reference point IS
    // assay(x, "mu"), so no reconstruction can drift from it.
    Rcpp::NumericMatrix mfull = Rcpp::NumericMatrix(mu);
    mfull.attr("dimnames") = mu.attr("dimnames");

    Rcpp::NumericVector bh(G), sh(G), bo(G), pr(G);
    for (int g = 0; g < G; ++g) {
        const std::size_t gk = static_cast<std::size_t>(g);
        bh[g] = Geo.beta_ref[gk];
        pr[g] = Geo.precision[gk];
        sh[g] = (nb_finite(Geo.precision[gk]) && Geo.precision[gk] > 0.0)
                    ? 1.0 / std::sqrt(Geo.precision[gk])
                    : NA_REAL;
        bo[g] = beta_obs[gk];
    }
    {
        Rcpp::CharacterVector rn = nb_rownames(mu, G);
        if (rn.size() == G) {
            bh.attr("names") = rn;
            sh.attr("names") = rn;
            bo.attr("names") = rn;
            pr.attr("names") = rn;
        }
    }

    return Rcpp::List::create(Rcpp::Named("mu_h0") = mh0,
                              Rcpp::Named("mu_full") = mfull,
                              Rcpp::Named("beta") = bh,
                              Rcpp::Named("se") = sh,
                              Rcpp::Named("precision") = pr,
                              Rcpp::Named("beta_obs") = bo,
                              Rcpp::Named("n_gene_na") = Geo.n_gene_na);
}

//' Resample bootstrap replicates of the counts from a given null mean
//'
//' Runs the engine's own resampling loop \code{b} times from the caller's RNG
//' stream so the draw order and the number of draws per replicate can be checked
//' against vectorised R \code{rnbinom()} calls under an identical seed.
//'
//' @param null_mu G x n mean matrix to draw from (pass \code{mu_h0}).
//' @param alpha G dispersions (NA -> 1, <= 0 -> Poisson).
//' @param b replicate number to return (>= 1).
//' @return a G x n integer matrix of resampled counts.
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerMatrix st_nb_resample_cpp(Rcpp::NumericMatrix null_mu,
                                       Rcpp::NumericVector alpha,
                                       int b) {
    const int G = null_mu.nrow();
    const int n = null_mu.ncol();
    if (b < 1)
        Rcpp::stop("b must be >= 1");
    if (alpha.size() != G)
        Rcpp::stop("alpha must have one entry per gene");

    std::vector<double> size_g;
    std::vector<char> poisson_g;
    nb_sizes(alpha.begin(), G, size_g, poisson_g);

    std::vector<double> cd(static_cast<std::size_t>(G) * n, 0.0);
    Rcpp::IntegerMatrix out(G, n);

    GetRNGstate();
    for (int rep = 0; rep < b; ++rep)
        nb_resample(cd.data(), G, n, null_mu.begin(), size_g, poisson_g);
    PutRNGstate();

    for (int j = 0; j < n; ++j)
        for (int g = 0; g < G; ++g)
            out(g, j) = static_cast<int>(cd[static_cast<R_xlen_t>(j) * G + g]);
    out.attr("dimnames") = null_mu.attr("dimnames");
    return out;
}

//' Engine pipeline for arbitrary counts -- the pipeline-identity probe
//'
//' Runs the SAME internal construction the bootstrap loop uses (lib.size -> E ->
//' delta weights -> st_effects_cpp) on whatever counts it is handed, and returns
//' every intermediate.  Feeding it the OBSERVED counts and comparing against
//' \code{as_dsea()}'s own \code{E}/\code{weights} proves the observed pass is
//' bit-identical; feeding it a RESAMPLED replicate and comparing against the R
//' formulas proves the null traverses the identical path.
//'
//' @param counts G x n integer matrix (observed or a replicate).
//' @param mu G x n fitted means (for the weights only).
//' @param alpha G dispersions (for the weights only).
//' @param norm_factors n size factors.
//' @param design n x p design matrix.
//' @param contrast_idx 1-based contrast column.
//' @param mu_floor lower bound for the weights.
//' @return a list with \code{lib_size}, \code{E}, \code{weights}, \code{Effects}.
//' @noRd
// [[Rcpp::export]]
Rcpp::List st_nb_pipeline_cpp(Rcpp::IntegerMatrix counts,
                              Rcpp::NumericMatrix mu,
                              Rcpp::NumericVector alpha,
                              Rcpp::NumericVector norm_factors,
                              Rcpp::NumericMatrix design,
                              int contrast_idx,
                              double mu_floor) {
    const int G = counts.nrow();
    const int n = counts.ncol();
    const int p = design.ncol();
    nb_validate(counts, mu, alpha, norm_factors, design, contrast_idx, G, n, p);
    if (!R_FINITE(mu_floor) || mu_floor < 0.0)
        Rcpp::stop("mu_floor must be a single finite non-negative number");

    std::vector<double> cd;
    nb_counts_to_double(counts, G, n, cd);

    std::vector<double> lib_size;
    nb_lib_size(cd.data(), G, n, norm_factors.begin(), lib_size);

    Rcpp::NumericMatrix E(G, n);
    nb_build_E(cd.data(), G, n, lib_size.data(), E.begin());
    E.attr("dimnames") = counts.attr("dimnames");

    std::vector<double> w_vals;
    nb_weights(mu.begin(), alpha.begin(), G, n, mu_floor, w_vals);
    Rcpp::NumericMatrix W(G, n);
    std::copy(w_vals.begin(), w_vals.end(), W.begin());
    W.attr("dimnames") = counts.attr("dimnames");

    Rcpp::NumericMatrix Eff = st_effects_cpp(E, W, design, contrast_idx);

    return Rcpp::List::create(
        Rcpp::Named("lib_size") = Rcpp::wrap(lib_size),
        Rcpp::Named("E") = E,
        Rcpp::Named("weights") = W,
        Rcpp::Named("Effects") = Eff);
}

//' The H0 null distribution of the set statistic
//'
//' Runs the same pipeline and returns the B replicate values instead of only the
//' exceedance counts, so the R side can measure the shape of the null.  The sd of
//' z_S under this null IS the hypothesis: ~1 means the analytic normal null
//' already suffices and the bootstrap adds nothing (criterion F1), while sd > 1 is
//' the excess the design exists to calibrate.
//'
//' @inheritParams st_nb_boot_cpp
//' @param statistic \code{"nb_z"} or \code{"fry_t"}.
//' @return a list with \code{t_obs}, \code{null} (nsets x B matrix of replicate
//'   statistics), \code{p_two}, \code{statistic} and \code{B}.
//' @noRd
// [[Rcpp::export]]
Rcpp::List st_nb_null_draws_cpp(Rcpp::IntegerMatrix counts,
                                Rcpp::NumericMatrix mu,
                                Rcpp::NumericVector alpha,
                                Rcpp::NumericVector norm_factors,
                                Rcpp::NumericMatrix design,
                                int contrast_idx,
                                Rcpp::List index,
                                int B,
                                Rcpp::Nullable<Rcpp::NumericVector> gene_weights,
                                double mu_floor,
                                std::string statistic,
                                Rcpp::NumericVector beta,
                                Rcpp::NumericVector se) {
    const int G = counts.nrow();
    const int n = counts.ncol();
    const int p = design.ncol();
    nb_validate(counts, mu, alpha, norm_factors, design, contrast_idx, G, n, p);
    if (index.size() < 1)
        Rcpp::stop("index is empty");
    if (statistic != "nb_z" && statistic != "fry_t")
        Rcpp::stop("statistic must be nb_z or fry_t");
    if (B < 1)
        Rcpp::stop("B must be a positive integer");
    if (!R_FINITE(mu_floor) || mu_floor < 0.0)
        Rcpp::stop("mu_floor must be a single finite non-negative number");
    if (beta.size() != G || se.size() != G)
        Rcpp::stop("beta and se must have one entry per gene");

    const bool is_z = statistic == "nb_z";
    std::vector<std::vector<int>> set_rows(static_cast<std::size_t>(index.size()));
    for (int i = 0; i < index.size(); ++i)
        nb_resolve_set(Rcpp::IntegerVector(index[i]), i, G,
                       set_rows[static_cast<std::size_t>(i)]);

    std::vector<double> lib_size, alpha_eff;
    nb_alpha_effective(alpha, G, alpha_eff);
    NbDesign D;
    D.init(design, norm_factors.begin(), contrast_idx);
    NbGeometry Geo;
    Geo.compute(mu.begin(), alpha_eff.data(), D, G, n);

    std::vector<double> cd;
    nb_counts_to_double(counts, G, n, cd);
    std::vector<double> mu_h0;
    nb_build_null_mu(mu.begin(), Geo, D, G, n, mu_h0);
    nb_check_null_consistency(mu.begin(), mu_h0.data(), Geo, D, G, n);

    std::vector<double> w_vals;
    Rcpp::NumericMatrix Wmat;
    if (!is_z) {
        nb_weights(mu.begin(), alpha.begin(), G, n, mu_floor, w_vals);
        Wmat = Rcpp::NumericMatrix(G, n);
        std::copy(w_vals.begin(), w_vals.end(), Wmat.begin());
    }
    std::vector<double> size_g;
    std::vector<char> poisson_g;
    nb_sizes(alpha.begin(), G, size_g, poisson_g);

    const int nres = n - p;
    const int nsets = index.size();
    Rcpp::NumericMatrix null_draws(nsets, B);
    Rcpp::NumericMatrix Eb(G, n);
    std::vector<double> Mbuf(static_cast<std::size_t>(nres) + 1, 0.0);
    std::vector<double> beta_star(static_cast<std::size_t>(G), NA_REAL);
    std::vector<double> t_obs(static_cast<std::size_t>(nsets), NA_REAL);
    std::vector<double> counters(static_cast<std::size_t>(nsets), 0.0);

    GetRNGstate();
    Rcpp::NumericMatrix Eff;
    for (int b = 0; b <= B; ++b) {
        if (b > 0)
            nb_resample(cd.data(), G, n, mu_h0.data(), size_g, poisson_g);
        if (is_z) {
            nb_closed_form(cd.data(), D, Geo, G, n, beta_star.data());
        } else {
            nb_lib_size(cd.data(), G, n, norm_factors.begin(), lib_size);
            nb_build_E(cd.data(), G, n, lib_size.data(), Eb.begin());
            Eff = st_effects_cpp(Eb, Wmat, design, contrast_idx);
        }
        for (int i = 0; i < nsets; ++i) {
            const double t = is_z
                ? nb_set_z(beta_star.data(), Geo.precision.data(),
                           set_rows[static_cast<std::size_t>(i)], nullptr)
                : nb_set_t(Eff, set_rows[static_cast<std::size_t>(i)],
                           std::vector<double>(), false, nres, Mbuf);
            if (b == 0) {
                t_obs[static_cast<std::size_t>(i)] = t;
            } else {
                null_draws(i, b - 1) = t;
                if (R_FINITE(t_obs[static_cast<std::size_t>(i)]) && R_FINITE(t) &&
                    std::fabs(t) >= std::fabs(t_obs[static_cast<std::size_t>(i)]))
                    counters[static_cast<std::size_t>(i)] += 1.0;
            }
        }
    }
    PutRNGstate();

    Rcpp::NumericVector p_two(nsets);
    Rcpp::NumericVector tob(nsets);
    for (int i = 0; i < nsets; ++i) {
        tob[i] = t_obs[static_cast<std::size_t>(i)];
        p_two[i] = R_FINITE(t_obs[static_cast<std::size_t>(i)])
                       ? (1.0 + counters[static_cast<std::size_t>(i)]) /
                             (static_cast<double>(B) + 1.0)
                       : NA_REAL;
    }
    colnames(null_draws) = Rcpp::CharacterVector::create();

    return Rcpp::List::create(Rcpp::Named("t_obs") = tob,
                              Rcpp::Named("null") = null_draws,
                              Rcpp::Named("p_two") = p_two,
                              Rcpp::Named("statistic") = statistic,
                              Rcpp::Named("B") = B);
}
