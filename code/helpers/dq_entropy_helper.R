# ─────────────────────────────────────────────────────────────────────────────
# dq_entropy_helper.R
# Utilities to compute entropy-based data-quality indices by county
#
# Exported:
#   • compute_entropy_county(ds, county_var)
#
# Returns a tibble with one row per county present in `ds` and columns:
#   - <county_var>
#   - DQ_overall    : normalized Shannon entropy in [0,1]
#   - DQ_H          : raw Shannon entropy (nats; ln)
#   - DQ_K          : number of distinct underlying-cause categories (root-3)
#   - DQ_eff_causes : exp(H) — effective number of causes
#
# Design notes
# • This implementation is glmnet-free and cannot trigger x/y length mismatches.
# • Works whether `uc3` already exists in `ds` or not (will derive if missing).
# • Safe on sparse counties (K=0 or K=1 handled gracefully).
# • Independent of `year` so it can be joined by county only (your script joins
#   the result of a single-year call to county-year metrics by county).
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")
    if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required")
    if (!requireNamespace("stringr", quietly = TRUE)) stop("stringr required")
    if (!requireNamespace("tibble", quietly = TRUE)) stop("tibble required")
    if (!requireNamespace("purrr", quietly = TRUE)) stop("purrr required")
})

# ---- helpers ----------------------------------------------------------------

# Canonicalize ICD-10 codes: uppercase, remove non [A-Z0-9]
canonical_icd_simple <- function(x) {
    stringr::str_remove_all(stringr::str_to_upper(x), "[^A-Z0-9]")
}

# Extract "root-3" category (letter + 2 digits) from ICD-10
# If input already looks like a root (length >= 3), will take first 3
icd_root3 <- function(x) {
    x <- canonical_icd_simple(x)
    stringr::str_sub(x, 1L, 3L)
}

# Shannon entropy (natural log); vector p must sum to 1
H_shannon <- function(p) {
    p <- p[p > 0 & is.finite(p)]
    if (length(p) == 0) return(0)
    -sum(p * log(p))
}

# ---- main API ----------------------------------------------------------------

#' Compute entropy-based county data-quality index (per-call year scope)
#'
#' @param ds data.frame (or tibble) containing at least:
#'        - `ucod` (underlying cause, string) OR `uc3` (precomputed root-3)
#'        - county id column named by `county_var`
#' @param county_var string name of the county identifier column in `ds`
#'
#' @return tibble with columns: county_var, DQ_overall, DQ_H, DQ_K, DQ_eff_causes
#'
compute_entropy_county <- function(ds, county_var) {
    stopifnot(is.data.frame(ds), is.character(county_var), length(county_var) == 1L)
    if (!county_var %in% names(ds)) {
        stop("compute_entropy_county: county_var '", county_var, "' not found in ds.")
    }
    
    # Prefer existing uc3 if present; else derive from ucod
    if ("uc3" %in% names(ds)) {
        uc3 <- ds[["uc3"]]
    } else if ("ucod" %in% names(ds)) {
        uc3 <- icd_root3(ds[["ucod"]])
    } else {
        stop("compute_entropy_county: require 'uc3' or 'ucod' in ds.")
    }
    
    # Build a compact table: county, uc3 (drop NA/empty)
    comp <- tibble::tibble(
        county = ds[[county_var]],
        uc3 = uc3
    ) %>%
        dplyr::filter(!is.na(county), !is.na(uc3), uc3 != "")
    
    if (nrow(comp) == 0L) {
        # No valid rows: return empty tibble with expected cols
        out <- tibble::tibble(
            !!county_var := character(),
            DQ_overall = numeric(),
            DQ_H = numeric(),
            DQ_K = integer(),
            DQ_eff_causes = numeric()
        )
        return(out)
    }
    
    # Count per county x uc3, then compute entropy per county
    by_county <- comp %>%
        dplyr::count(.data$county, .data$uc3, name = "n_uc3") %>%
        dplyr::group_by(.data$county) %>%
        dplyr::summarise(
            N      = sum(.data$n_uc3),
            K      = dplyr::n_distinct(.data$uc3),
            H      = {
                p <- .data$n_uc3 / sum(.data$n_uc3)
                H_shannon(p)
            },
            .groups = "drop"
        ) %>%
        dplyr::mutate(
            # Normalize by H_max = ln(K); define 0 when K <= 1
            DQ_overall   = dplyr::if_else(K > 1, H / log(K), 0),
            DQ_H         = H,
            DQ_K         = as.integer(K),
            DQ_eff_causes = exp(H)
        ) %>%
        dplyr::select(.data$county, .data$DQ_overall, .data$DQ_H, .data$DQ_K, .data$DQ_eff_causes)
    
    # Rename back to county_var expected by caller
    names(by_county)[names(by_county) == "county"] <- county_var
    by_county
}

# ─────────────────────────────────────────────────────────────────────────────
# Optional: multinomial/glmnet debug stub (disabled by default)
# If you ever restore a glmnet-based entropy model, call `align_xy_complete_cases()`
# to guarantee X and y are aligned, preventing the classic mismatch.
# ─────────────────────────────────────────────────────────────────────────────
# align_xy_complete_cases <- function(X, y) {
#   if (!is.matrix(X)) X <- as.matrix(X)
#   keep <- stats::complete.cases(cbind(y, X))
#   list(X = X[keep, , drop = FALSE], y = y[keep])
# }
#
# dq_glmnet_debug <- function(x, y, ...) {
#   if (!is.matrix(x)) x <- as.matrix(x)
#   ok <- stats::complete.cases(cbind(y, x))
#   if (sum(ok) < length(y)) {
#     cat("[dq_glmnet_debug] Dropping", sum(!ok), "rows with NA/Inf in X or y\n")
#   }
#   glmnet::glmnet(x = x[ok, , drop = FALSE], y = y[ok], ...)
# }

