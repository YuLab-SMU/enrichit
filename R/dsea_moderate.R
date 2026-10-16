## Design fitting, empirical Bayes variance moderation and consensus correlation.
##
## WHY MODERATE AT ALL.  With a handful of samples the per-pathway residual
## variance carries about n - rank(X) degrees of freedom -- five, in the case
## this package is built for.  A variance estimated from five residuals is
## mostly noise, and a t statistic built on it has a much heavier tail than its
## nominal distribution.  Moderation borrows strength across pathways: the
## residual variances are shrunk towards a common fitted prior, which is what
## makes the test usable at small n.
##
## The estimator is Smyth (2004)'s scaled-inverse-chi-square prior on the log
## variances.  It is re-implemented here in base R so that the package stays
## free of limma; `tests/testthat/test-dsea-pathway.R` pins it against
## limma::eBayes on the same input, so the re-implementation is checked rather
## than merely claimed.

## ---- design ----------------------------------------------------------------

.dsea_check_design <- function(design, samples) {
    if (is.null(design)) stop("'design' is required.", call. = FALSE)
    design <- as.matrix(design)
    if (!is.numeric(design)) {
        stop("'design' must be a numeric matrix. Build it with model.matrix(), ",
             "e.g. model.matrix(~ condition, data = colData); got ", class(design)[1], ".",
             call. = FALSE)
    }
    if (any(!is.finite(design))) stop("'design' contains non-finite values.", call. = FALSE)
    if (nrow(design) != length(samples)) {
        stop("'design' has ", nrow(design), " rows but 'x' has ", length(samples),
             " samples.", call. = FALSE)
    }
    if (is.null(colnames(design))) {
        colnames(design) <- paste0("x", seq_len(ncol(design)))
        warning("'design' had no column names; referring to columns by index only.",
                call. = FALSE)
    }
    if (!is.null(rownames(design)) && !identical(rownames(design), samples)) {
        stop("rownames(design) do not match colnames(x); reorder the design to match the samples.",
             call. = FALSE)
    }
    rank <- qr(design)$rank
    if (rank < ncol(design)) {
        ## Checked BEFORE any fitting, and named explicitly: a rank-deficient
        ## design is the airway dex/cell trap, and it must not surface as a
        ## downstream numerical error.
        stop("'design' is not full rank (rank ", rank, " of ", ncol(design),
             " columns). One or more design terms are linear combinations of ",
             "the others and cannot be estimated; drop the redundant terms first.",
             call. = FALSE)
    }
    design
}

.dsea_resolve_contrasts <- function(contrast, design) {
    if (is.null(contrast)) {
        ## Everything except a leading intercept.  An intercept is only skipped
        ## when the design actually has one (first column is constant); a design
        ## without an intercept has no column to skip, and skipping its first
        ## column would silently drop a real contrast.
        has_int <- ncol(design) > 1L &&
            isTRUE(all.equal(as.numeric(design[, 1L]), rep(1, nrow(design))))
        idx <- if (has_int) seq_len(ncol(design))[-1L] else seq_len(ncol(design))
        if (length(idx) == 0L) {
            stop("'design' has only an intercept; there is nothing to test.", call. = FALSE)
        }
        return(idx)
    }
    if (length(contrast) > 1L) {
        return(vapply(contrast, function(c0) .dsea_resolve_one(c0, design), integer(1)))
    }
    .dsea_resolve_one(contrast, design)
}

.dsea_resolve_one <- function(contrast, design) {
    if (is.character(contrast)) {
        i <- match(contrast, colnames(design))
        if (is.na(i)) {
            stop("contrast '", contrast, "' not found in design; available: ",
                 paste(colnames(design), collapse = ", "), call. = FALSE)
        }
        return(i)
    }
    i <- as.integer(contrast)
    if (length(i) != 1L || is.na(i) || i < 1L || i > ncol(design)) {
        stop("'contrast' must be a single column name or index in 1..", ncol(design),
             call. = FALSE)
    }
    i
}

## ---- standardisation -------------------------------------------------------

.dsea_standardize <- function(S) {
    mu <- rowMeans(S)
    sdv <- apply(S, 1L, stats::sd)
    bad <- !is.finite(sdv) | sdv <= 0
    Z <- (S - mu) / sdv
    if (any(bad)) Z[bad, ] <- NA_real_
    attr(Z, "dropped_constant") <- rownames(S)[bad]
    Z
}

## ---- fitting ---------------------------------------------------------------

.dsea_fit <- function(Z, design, ci, moderate) {
    n <- nrow(design)
    Y <- t(Z)                              # samples x pathways
    ## LAPACK = TRUE is essential, not stylistic.  qr() pivots columns by
    ## default, and qr.R() comes back with those pivot-induced sign flips, so
    ## R[ci, ci] and solve(R) refer to a different parameterisation than the
    ## design columns the caller named.  The result is a coefficient of the
    ## wrong sign that looks entirely plausible.  Without pivoting, X = Q R
    ## exactly and column i of R is design column i.
    qrX <- qr(design, LAPACK = TRUE)
    rank <- ncol(design)
    dof <- n - rank
    ## qr.Q on a LAPACK factorisation returns the full n x n Q, so the first
    ## `rank` columns are the ones belonging to the design.
    Q <- qr.Q(qrX)[, seq_len(rank), drop = FALSE]
    R <- qr.R(qrX)[seq_len(rank), seq_len(rank), drop = FALSE]
    ## solve(), not chol2inv().  chol2inv() inverts a Cholesky factor, i.e. an
    ## upper-triangular matrix with a POSITIVE diagonal; qr.R() has whatever
    ## sign the QR produced (negative here), and chol2inv() silently returns
    ## R^-1 R^-T for that matrix instead of R^-1.  Every coefficient came out
    ## wrong; only the golden test against limma exposed it.
    Rinv <- solve(R)
    ## The contrast variance multiplier is [(X'X)^-1]_{ci,ci} = [R^-1 R^-T]_{ci,ci},
    ## NOT [R^-1]_{ci,ci}.  Those are different numbers (0.5 vs -0.707 on a
    ## balanced two-group design) and using the wrong one makes every se a NaN.
    RinvXtRinv <- Rinv %*% t(Rinv)
    v <- RinvXtRinv[ci, ci]                # contrast variance multiplier
    QtY <- crossprod(Q, Y)                 # t(Q) Y, rank x pathways
    coef <- Rinv %*% QtY                   # rank x pathways
    fitted <- Q %*% coef
    beta <- coef[ci, ]
    resid <- Y - fitted
    s2 <- colSums(resid^2) / dof
    mod <- .dsea_moderate(s2, v, dof, moderate)
    se <- sqrt(mod$se2 * v)
    tstat <- beta / se
    names(beta) <- names(se) <- names(tstat) <- rownames(Z)
    list(logFC = beta, se = se, t = tstat, df.total = mod$df.total,
         dof.residual = dof, moderation = mod)
}

## ---- empirical Bayes moderation (Smyth 2004) --------------------------------

## log(x) - digamma(x), the "mean digamma" term.  statmod computes it with a
## shift for small arguments; at the arguments we use here (dof/2 >= 1.5) the
## direct form is accurate to double precision, and the golden test against
## limma::eBayes is what actually holds us to that.
.dsea_logmdigamma <- function(x) log(x) - digamma(x)

## Inverse of the trigamma function, by Newton iteration (as in Smyth 2004).
.dsea_trigamma_inverse <- function(x) {
    if (anyNA(x)) return(ifelse(is.na(x), NA_real_, .dsea_trigamma_inverse(x)))
    if (any(x < 0)) stop("'x' must be non-negative.", call. = FALSE)
    x <- x[x > 0]
    if (length(x) == 0L) return(numeric(0))
    y <- ifelse(x > 1e7, 1 / sqrt(x), ifelse(x < 1e-6, 1 / x, 0.5 + 1 / x))
    for (i in seq_len(100L)) {
        tri <- trigamma(y)
        dif <- tri * (1 - tri / x) / psigamma(y, deriv = 2)
        y <- y + dif
        if (max(-dif / y) < 1e-8) break
    }
    y
}

## Moment estimator for the scaled-inverse-chi-square prior on log(s2).
## Reproduces limma::fitFDist(x, df1 = dof) with no trend covariate.
.dsea_fit_prior <- function(s2, dof) {
    n <- length(s2)
    ok <- is.finite(s2) & s2 > -1e-15
    nok <- sum(ok)
    if (nok == 0L) return(list(scale = NA_real_, df2 = NA_real_))
    if (nok == 1L) return(list(scale = s2[ok], df2 = 0))
    x <- pmax(s2[ok], 0)
    m <- stats::median(x)
    if (m == 0) {
        warning("More than half of the residual variances are exactly zero: moderation is unreliable.",
                call. = FALSE)
        m <- 1
    }
    x <- pmax(x, 1e-5 * m)
    e <- log(x) + .dsea_logmdigamma(dof / 2)
    emean <- mean(e)
    evar <- sum((e - emean)^2) / (nok - 1)
    evar <- evar - trigamma(dof / 2)
    if (evar > 0) {
        df2 <- 2 * .dsea_trigamma_inverse(evar)
        s20 <- exp(emean - .dsea_logmdigamma(df2 / 2))
    } else {
        df2 <- Inf
        s20 <- mean(x)
    }
    list(scale = s20, df2 = df2)
}

.dsea_moderate <- function(s2, v, dof, moderate) {
    if (dof <= 0) {
        stop("the design leaves no residual degrees of freedom (n = p); ",
             "nothing can be tested. Add samples or simplify the design.",
             call. = FALSE)
    }
    if (moderate == "none") {
        return(list(se2 = s2, df.total = dof, prior_scale = NA_real_, prior_df = Inf,
                    method = "none"))
    }
    prior <- .dsea_fit_prior(s2, dof)
    d0 <- prior$df2
    s0 <- prior$scale
    if (!is.finite(d0)) {
        se2 <- rep(s0, length(s2))
        df.post <- dof + d0
    } else {
        ## (df * s2 + df.prior * s0) / (df + df.prior), exactly as limma's
        ## .squeezeVar does.  An earlier revision wrote (d0 * s2 + s0) here:
        ## the weights were swapped, so the poorly-estimated variance got the
        ## LARGE weight and the prior was added with no degree-of-freedom
        ## weighting at all.  The shrinkage pulled every pathway towards the
        ## prior less when the pathway's own estimate was least trustworthy --
        ## the exact opposite of the intent -- and no test covered it, because
        ## the golden tests pinned fitFDist and the OLS coefficients but not
        ## this line.
        se2 <- (dof * s2 + d0 * s0) / (dof + d0)
        df.post <- dof + d0
    }
    se2[!is.finite(se2) | se2 < 0] <- NA_real_
    list(se2 = se2, df.total = df.post, prior_scale = s0, prior_df = d0,
         method = "eb")
}

## ---- consensus correlation -------------------------------------------------

## Closed-form equicorrelation estimate.  Under a consensus model the residual
## matrix R (pathways x residual df) has row covariance
## sigma2 * ((1 - rho) I + rho J), whose eigenvalues are
## sigma2 (1 + (m-1) rho) once and sigma2 (1 - rho) m-1 times.  The singular
## values of R estimate those eigenvalues (times the number of columns), so:
##     sigma2 = sum_{i>=2} d_i^2 / ((m-1) * ncol)
##     u      = (d_1^2 / ncol) / sigma2
##     rho    = (u - 1) / (u + m - 1)
##
## NOTE this is NOT limma::duplicateCorrelation, which fits per-pathway
## variance components iteratively.  This is our own estimator; the benchmark
## reports limma's value alongside it rather than claiming identity (see
## PATHWAY_DSEA_SPEC.md A11).
.dsea_consensus_correlation <- function(Z, design, block = NULL) {
    m <- nrow(Z)
    nres <- ncol(design) - qr(design)$rank
    if (m < 3L || nres < 1L) return(0)
    X <- design
    if (!is.null(block)) {
        if (max(table(block)) == 1L) {
            warning("Blocks are all of size 1: intrablock correlation set to zero.", call. = FALSE)
            return(0)
        }
        X <- cbind(X, model.matrix(~ factor(block))[, -1, drop = FALSE])
    }
    qrX <- qr(X)
    if (qrX$rank < ncol(X)) {
        ## Block already encoded in the design: nothing left to estimate.
        return(0)
    }
    R <- t(apply(Z, 1L, function(y) qr.resid(qrX, y)))
    R <- R[stats::complete.cases(R), , drop = FALSE]
    if (nrow(R) < 3L) return(0)
    d <- svd(R, nu = 0L, nv = 0L)$d
    nres2 <- ncol(R)
    denom <- (length(d) - 1) * nres2
    if (denom <= 0) return(0)
    sigma2 <- (sum(d^2) - d[1]^2) / denom
    if (!is.finite(sigma2) || sigma2 <= 0) return(0)
    u <- (d[1]^2 / nres2) / sigma2
    rho <- (u - 1) / (u + length(d) - 1)
    if (!is.finite(rho)) return(0)
    ## Negative rho is not meaningful for a consensus that will only ever
    ## inflate significance; report zero rather than a negative variance factor.
    min(max(rho, 0), 0.99)
}

## ---- scale guard -----------------------------------------------------------

## The previous dsea() shipped a z whose magnitude was ~1e4 because the
## numerator was a weighted SUM and the denominator the variance of the
## weighted MEAN.  It changed no p-value, but it made every |z| threshold and
## every effect size meaningless.
##
## The threshold here is deliberately far above anything real.  A standardised
## t is O(1) under the null, but with strong genuine signal |t| can legitimately
## reach tens or hundreds, and a guard that fires on those would reject correct
## results.  1e3 sits above every plausible biological t and below the 1e4 that
## the old scale error produced, so it catches that class of bug and nothing
## else.  The tight O(1) claim is checked under the null in the unit tests
## (PATHWAY_DSEA_SPEC.md A2), where it belongs.
.dsea_guard_scale <- function(contrasts, tol = 1e3) {
    for (r in contrasts) {
        mt <- stats::median(abs(r$t[is.finite(r$t)]))
        ml <- stats::median(abs(r$logFC[is.finite(r$logFC)]))
        if (!is.finite(mt) || mt > tol) {
            stop("internal scale guard: median |t| = ", signif(mt, 4),
                 " exceeds ", tol, ". The statistic is not on a t scale; ",
                 "no result is returned because it would be uninterpretable.",
                 call. = FALSE)
        }
        if (!is.finite(ml) || ml > tol) {
            stop("internal scale guard: median |logFC| = ", signif(ml, 4),
                 " exceeds ", tol, ". The effect size is not standardised; ",
                 "no result is returned because it would be uninterpretable.",
                 call. = FALSE)
        }
    }
    invisible(TRUE)
}