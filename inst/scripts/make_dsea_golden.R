#!/usr/bin/env Rscript
## ===========================================================================
## make_dsea_golden.R -- generate tests/refdata/dsea_golden.rds
##
## WHAT THIS IS
##   The golden reference that pins the numerical output of dsea() for the
##   camera / roast / fry engines.  tests/testthat/test-dsea.R reads the .rds
##   and asserts that a freshly computed run reproduces every locked p-value.
##   limma was removed from the package (REBUILD_SPEC section 1), so this file
##   is the ONLY thing standing between a silent engine change and a wrong
##   p-value: there is no external reference implementation to fall back on.
##
## HOW TO REGENERATE  (deliberately, never casually)
##   1. Be certain the engine is FINAL.  REBUILD_SPEC section 5: the golden may
##      only be produced by the frozen engine, and only after the parity
##      scenarios pass.  Regenerating it while an engine bug is still in the
##      tree re-freezes the bug into the reference.
##   2. Run, from the package root (or from anywhere -- the output path is
##      resolved from this script's own location):
##          Rscript inst/scripts/make_dsea_golden.R
##   3. Inspect the printed table against the previous one.  If ANY p-value
##      moved, that change is the thing you must justify in the commit message
##      -- do not regenerate to make a failing test go away.
##   4. Re-run the full test suite:
##          Rscript -e 'suppressPackageStartupMessages(devtools::load_all(".", quiet=TRUE)); library(testthat); test_dir("tests/testthat")'
##
## WHAT IS LOCKED
##   One p-value vector per method (camera, fry, roast), for six gene sets of
##   increasing size, on the same fitted DESeqDataSet every time:
##
##       set.seed(42)
##       dds <- DESeq2::makeExampleDESeqDataSet(n = 2000, m = 8)
##       dds <- DESeq2::DESeq(dds)
##
##   The contrast is taken from the LAST column of as_dsea()$design rather than
##   hard-coded: makeExampleDESeqDataSet() names it "condition_B_vs_A", whereas
##   the airway example in the documentation names it "dex_trt_vs_untrt".  Both
##   are spelled by base model.matrix()/DESeq2, not by this package, so the name
##   must be discovered at run time (REBUILD_SPEC trap 8.4).
##
##   The gene sets are slices of rownames(E) in a fixed order, so building them
##   consumes no RNG.  The only RNG consumer is roast's rotation loop, and each
##   dsea() call is preceded by its own set.seed(42); the engine deliberately
##   does not seed itself (see dsea()'s "Resampling methods" section).
##
## SIDE EFFECTS
##   Writes exactly one file: <package>/tests/refdata/dsea_golden.rds.
##   No other file is created or modified.  Safe to re-run.
## ===========================================================================

local({
    ## ---- locate the package root by walking up from this script -------------
    args <- commandArgs(trailingOnly = FALSE)
    fpath <- sub("^--file=", "", args[grep("^--file=", args)][1])
    here <- if (nzchar(fpath)) dirname(normalizePath(fpath)) else getwd()
    root <- here
    for (i in 1:6) {
        if (file.exists(file.path(root, "DESCRIPTION"))) break
        root <- dirname(root)
    }
    if (!file.exists(file.path(root, "DESCRIPTION"))) {
        stop("could not locate the package root from ", here,
             "; run this script from inside the enrichit source tree.",
             call. = FALSE)
    }
    out_dir <- file.path(root, "tests", "refdata")
    out_file <- file.path(out_dir, "dsea_golden.rds")

    cat(sprintf("package root : %s\n", root))
    cat(sprintf("output file  : %s\n", out_file))

    ## ---- the engine --------------------------------------------------------
    suppressPackageStartupMessages({
        library(devtools)
        library(DESeq2)
    })
    ## load_all() so that the .so in src/ is what gets exercised; installing the
    ## released package here would silently test a different build.
    devtools::load_all(root, quiet = TRUE)

    ## ---- the data (fixed seed, fixed shape) --------------------------------
    SEED <- 42L
    N_GENES <- 2000L
    M_SAMPLES <- 8L
    set.seed(SEED)
    dds <- DESeq2::makeExampleDESeqDataSet(n = N_GENES, m = M_SAMPLES)
    dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))

    dd <- enrichit::as_dsea(dds)
    genes <- rownames(dd$E)
    contrast <- colnames(dd$design)[ncol(dd$design)]

    ## Six sets, sizes 5 / 20 / 50 / 100 / 200 / 400.  Contiguous slices of the
    ## gene order, so construction is deterministic and RNG-free.
    gene_sets <- lapply(c(5, 20, 50, 100, 200, 400), function(k) {
        head(genes, k)
    })
    names(gene_sets) <- sprintf("set%02d", seq_along(gene_sets))

    cat(sprintf("data         : n=%d genes, m=%d samples, contrast = '%s'\n",
                N_GENES, M_SAMPLES, contrast))
    cat(sprintf("gene sets    : %s (genes per set: %s)\n",
                paste(names(gene_sets), collapse = ", "),
                paste(vapply(gene_sets, length, integer(1)), collapse = ", ")))

    ## ---- the three locked p-value vectors ----------------------------------
    METHODS <- c("camera", "roast", "fry")
    locked <- vector("list", length(METHODS))
    names(locked) <- METHODS
    locked_dir <- vector("list", length(METHODS))
    names(locked_dir) <- METHODS

    for (m in METHODS) {
        set.seed(SEED)
        ## sort = FALSE: the locked vectors are indexed by set NAME, not by
        ## p-value rank, so a change in the ranking of two sets must not look
        ## like a change in a p-value.  (dsea() sorts by PValue by default, and
        ## camera/fry/roast do not rank the same sets the same way.)
        res <- dsea(dds, gene_sets, contrast = contrast, method = m, sort = FALSE)
        stopifnot(
            identical(rownames(res), res$set),
            identical(res$set, names(gene_sets)),
            all(res$PValue >= 0 & res$PValue <= 1)
        )
        locked[[m]] <- res$PValue
        locked_dir[[m]] <- res$Direction
        cat(sprintf("  %-6s : %s\n", m,
                    paste(sprintf("%.10g", res$PValue), collapse = " ")))
    }

    golden <- list(
        ## provenance: everything needed to reproduce this file
        provenance = list(
            script        = "inst/scripts/make_dsea_golden.R",
            generated     = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
            r_version     = R.version.string,
            enrichit      = as.character(utils::packageVersion("enrichit")),
            deseq2        = as.character(utils::packageVersion("DESeq2")),
            seed          = SEED,
            n_genes       = N_GENES,
            m_samples     = M_SAMPLES,
            contrast      = contrast,
            generator     = "DESeq2::makeExampleDESeqDataSet + DESeq2::DESeq + as_dsea + dsea"
        ),
        set_names = names(gene_sets),
        n_genes_per_set = vapply(gene_sets, length, integer(1)),
        pvalue = locked,
        direction = locked_dir
    )

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(golden, out_file)
    cat(sprintf("\nwrote %s (%.1f kB)\n", out_file,
                file.info(out_file)$size / 1024))
    cat("Locked one p-value vector per method: camera, fry, roast.\n")
    cat("Re-read this file and re-run the test suite before committing.\n")
})
