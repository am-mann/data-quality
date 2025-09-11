# ───────────────── dq_entropy_helper.R — Foreman bins + RI (county×year) ─────────────────
# Inputs:
#   - ds: data frame with 'ucod', county key, 'year', optional 'age_years' or 'ager27',
#         and contributing causes in record_1..record_20 (or record1..record20)
#   - county_var: name of the county column in ds (e.g., "county_ihme")
#   - dict_dir: directory containing:
#       * foreman-icd10-mapping.csv  (ICD → USCOD codes, e.g. A_5, B_3_1, G_1)
#       * foreman-table2-map.csv     (wide: target_cause + G_1..G_9 TRUE/FALSE)
#
# Output columns (joined by county×year):
#   total, DQ_entropy, DQ_K, DQ_overall, DQ_expH_over_K,
#   RI, DQ_rec_expH_over_K_mean, DQ_rec_ig_abs_mean, DQ_rec_ig_frac_mean_garbage,
#   foreman_garbage, foreman_garbage_adj
#
# Options you can set (examples):
#   options(DQ_VERBOSE=TRUE, DQ_DIAG_DIR="output/dq_diag", DQ_SEED=123, DQ_POOL_YEARS=2)
# ──────────────────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr); library(Matrix)
})

# --------------------------- options / switches -------------------------------------------
.pool_years_default       <- function() as.integer(getOption("DQ_POOL_YEARS", 0L))
.dirichlet_prior_default  <- function() as.numeric(getOption("DQ_DIRICHLET_PRIOR", 0.5))
.min_per_class_default    <- function() as.integer(getOption("DQ_MIN_PER_CLASS", 12L))
.max_train_per_class_def  <- function() as.integer(getOption("DQ_MAX_TRAIN_PER_CLASS", 8000L))
.seed_default             <- function() as.integer(getOption("DQ_SEED", NA_integer_))
.verbose_default          <- function() isTRUE(getOption("DQ_VERBOSE", FALSE))
.diag_dir_default         <- function() { x <- getOption("DQ_DIAG_DIR", NULL); if (is.character(x) && nzchar(x)) x else NULL }
.has_glmnet <- function() requireNamespace("glmnet", quietly = TRUE)

# --------------------------- tiny helpers -------------------------------------------------
`%||%`        <- function(a,b) if (!is.null(a)) a else b
canonical_icd <- function(x) stringr::str_remove_all(stringr::str_to_upper(x), "[^A-Z0-9]")
icd3          <- function(x) substr(canonical_icd(x), 1, 3)
icd4          <- function(x) substr(canonical_icd(x), 1, 4)
norm_code     <- function(x) tolower(gsub("[^A-Za-z0-9]+", "", as.character(x)))  # "B_3_1"→"b31", "G_1"→"g1"
.shannon_H    <- function(p) { p <- p[p > 0]; if (!length(p)) 0 else -sum(p * log(p)) }

.rename_record_cols <- function(df) {
    if (!any(grepl("^record_", names(df))) && any(grepl("^record[0-9]+$", names(df)))) {
        dplyr::rename_with(df, ~ sub("^record([0-9]+)$","record_\\1", .x), .cols = dplyr::matches("^record[0-9]+$"))
    } else df
}

.read_icd_map <- function(path) {
    df <- readr::read_csv(
        path, show_col_types = FALSE,
        col_types = readr::cols(.default = readr::col_character(),
                                ICD10 = readr::col_character(),
                                USCOD = readr::col_character())
    )
    if (!all(c("ICD10","USCOD") %in% names(df))) stop("foreman-icd10-mapping.csv must have ICD10 and USCOD (e.g., A_5, B_3_1, G_1).")
    out <- df %>%
        transmute(icd = canonical_icd(ICD10), bin = norm_code(USCOD)) %>%
        filter(icd != "", bin != "") %>%
        distinct(icd, .keep_all = TRUE)
    if (any(grepl("^\\d+$", out$bin))) stop("USCOD parsed as digits-only; ensure alphanumeric codes like A_5, B_3_1, G_1.")
    out
}

.read_table2_map <- function(path) {
    df <- readr::read_csv(path, show_col_types = FALSE,
                          col_types = readr::cols(.default = readr::col_character()))
    if (!"target_cause" %in% names(df)) stop("foreman-table2-map.csv must contain 'target_cause'.")
    g_cols <- grep("^G_[0-9]+$", names(df), value = TRUE)
    if (!length(g_cols)) stop("foreman-table2-map.csv must contain columns G_1..G_9.")
    
    df[g_cols] <- lapply(df[g_cols], function(x){
        y <- tolower(as.character(x)); y %in% c("1","true","t","y","yes","x","✓","check","checked")
    })
    t2 <- df %>% mutate(target = norm_code(target_cause)) %>%
        pivot_longer(all_of(g_cols), names_to = "garbage_raw", values_to = "flag") %>%
        filter(flag) %>%
        transmute(garbage = norm_code(garbage_raw), target = norm_code(target)) %>%
        filter(garbage != "", target != "") %>% distinct(garbage, target)
    miss <- setdiff(paste0("g",1:9), unique(t2$garbage))
    if (length(miss)) stop("Missing garbage bins in Table 2: ", paste(miss, collapse=", "))
    t2
}

.make_icd_lookup <- function(icd_map) {
    list(
        map4 = icd_map %>% filter(nchar(icd) == 4) %>% select(icd4 = icd, bin4 = bin),
        map3 = icd_map %>% filter(nchar(icd) == 3) %>% select(icd3 = icd, bin3 = bin)
    )
}
map_icd_to_bin <- function(x, lu) {
    x4 <- icd4(x); x3 <- icd3(x)
    out <- rep(NA_character_, length(x))
    m4 <- match(x4, lu$map4$icd4); ok4 <- !is.na(m4); out[ok4] <- lu$map4$bin4[m4[ok4]]
    m3 <- match(x3, lu$map3$icd3); ok3 <- !is.na(m3) & !ok4; out[ok3] <- lu$map3$bin3[m3[ok3]]
    out
}

# Align columns to target classes (add zeros for missing; reorder; drop extras)
.align_cols <- function(M, classes) {
    have <- colnames(M) %||% character(0)
    keep <- intersect(classes, have)
    miss <- setdiff(classes, have)
    base <- if (length(keep)) M[, keep, drop = FALSE] else matrix(, nrow = nrow(M), ncol = 0)
    if (length(miss)) {
        base <- cbind(base, matrix(0, nrow = nrow(M), ncol = length(miss), dimnames = list(NULL, miss)))
    }
    base[, classes, drop = FALSE]
}

# Robustly extract a softmax probability matrix (n_te × k) from glmnet
.extract_multinomial_probs <- function(fit, X_te, classes, s = "lambda.1se") {
    out <- tryCatch(predict(fit, newx = X_te, type = "response", s = s), error = function(e) NULL)
    if (is.null(out)) return(NULL)
    
    # Case A: 3D array [n_te, k, nlambdas]
    if (is.array(out) && length(dim(out)) == 3L) {
        M <- out[, , 1, drop = TRUE]  # => matrix n_te × k
        if (!is.matrix(M)) M <- as.matrix(M)
        if (is.null(colnames(M)) || anyNA(colnames(M))) {
            cls <- tryCatch(fit$glmnet.fit$classnames, error = function(e) NULL)
            if (!is.null(cls) && length(cls) == ncol(M)) colnames(M) <- cls
        }
        return(.align_cols(M, classes))
    }
    
    # Case B: 2D matrix
    if (is.matrix(out) && length(dim(out)) == 2L) {
        if (is.null(colnames(out)) || anyNA(colnames(out))) {
            cls <- tryCatch(fit$glmnet.fit$classnames, error = function(e) NULL)
            if (!is.null(cls) && length(cls) == ncol(out)) colnames(out) <- cls
        }
        return(.align_cols(out, classes))
    }
    
    # Case C: list of per-class matrices/vectors
    if (is.list(out)) {
        mats <- lapply(out, function(x) {
            if (is.matrix(x)) x[, 1, drop = TRUE] else as.numeric(x)
        })
        M <- do.call(cbind, mats)
        colnames(M) <- names(out) %||% paste0("cls", seq_len(ncol(M)))
        return(.align_cols(as.matrix(M), classes))
    }
    
    # Fallback
    return(.align_cols(as.matrix(out), classes))
}

# Force/repair probability shape to (n_te × k). Try to do so diagnostically.
.fix_prob_shape <- function(M, n_te, classes, context = NULL, prefer_row_sums_close_to_1 = TRUE) {
    k <- length(classes)
    # already fine
    if (is.matrix(M) && nrow(M) == n_te && ncol(M) == k) {
        if (is.null(colnames(M)) || anyNA(colnames(M))) colnames(M) <- classes
        return(list(M = M, note = "ok"))
    }
    
    dims <- dim(M)
    # transposed case
    if (is.matrix(M) && all(dims == c(k, n_te))) {
        M2 <- t(M); colnames(M2) <- classes
        return(list(M = M2, note = "transpose[k×n_te]→[n_te×k]"))
    }
    
    # flattened row-blowup: (n_te*k) × k   ← common glmnet quirk
    if (is.matrix(M) && dims[1] == n_te * k && dims[2] == k) {
        v <- as.vector(M)  # column-major
        # Try both reshape orientations and pick the one closer to row-sums=1
        cand1 <- matrix(v, nrow = n_te, ncol = k, byrow = FALSE)
        cand2 <- matrix(v, nrow = n_te, ncol = k, byrow = TRUE)
        r1 <- rowSums(cand1); r2 <- rowSums(cand2)
        dev1 <- mean(abs(r1 - 1), na.rm = TRUE); dev2 <- mean(abs(r2 - 1), na.rm = TRUE)
        if (prefer_row_sums_close_to_1 && dev2 < dev1) {
            colnames(cand2) <- classes
            return(list(M = cand2, note = sprintf("reshape[(n_te*k)×k] byrow=TRUE (dev=%.4g)", dev2)))
        } else {
            colnames(cand1) <- classes
            return(list(M = cand1, note = sprintf("reshape[(n_te*k)×k] byrow=FALSE (dev=%.4g)", dev1)))
        }
    }
    
    # completely flattened vector with length n_te*k
    if (is.null(dims) && length(M) == n_te * k) {
        M2 <- matrix(as.numeric(M), nrow = n_te, ncol = k, byrow = TRUE)
        colnames(M2) <- classes
        return(list(M = M2, note = "vector→[n_te×k] byrow=TRUE"))
    }
    
    stop("Probability matrix has dims ", paste(dims, collapse = "×"),
         " incompatible with n_te×k=", n_te, "×", k, if (!is.null(context)) paste0(" (", context, ")") else "")
}

# ------------------------------- main function -------------------------------------------
compute_entropy_county_foreman <- function(
        ds,
        county_var,
        dict_dir      = NULL,
        icd_map_path  = if (!is.null(dict_dir)) file.path(dict_dir, "foreman-icd10-mapping.csv") else NULL,
        code_map_path = if (!is.null(dict_dir)) file.path(dict_dir, "foreman-table2-map.csv") else NULL,
        age_breaks    = c(-Inf, 1, 5, seq(10, 85, by = 5), Inf)
) {
    if (is.null(icd_map_path) || is.null(code_map_path))
        stop("Provide icd_map_path and code_map_path (or dict_dir).")
    
    seed     <- .seed_default(); if (!is.na(seed)) set.seed(seed)
    alpha0   <- .dirichlet_prior_default()
    minpc    <- .min_per_class_default()
    capc     <- .max_train_per_class_def()
    pool_y   <- .pool_years_default()
    verbose  <- .verbose_default()
    diag_dir <- .diag_dir_default()
    diag_rows <- list()
    
    # ---------- dictionaries ----------
    icd_map  <- .read_icd_map(icd_map_path)
    t2_map   <- .read_table2_map(code_map_path)
    icd_lu   <- .make_icd_lookup(icd_map)
    
    garbage_bins     <- sort(unique(t2_map$garbage))
    all_bins_in_map  <- sort(unique(icd_map$bin))
    valid_bins       <- setdiff(all_bins_in_map, garbage_bins)
    K <- length(valid_bins)
    if (K < 10) warning("Few valid bins detected (K=", K, "). Check dictionaries.")
    
    bad_targets <- setdiff(unique(t2_map$target), all_bins_in_map)
    if (length(bad_targets))
        warning("Table-2 targets not in ICD map: ", paste(head(bad_targets,12), collapse=", "),
                if (length(bad_targets) > 12) paste0(" … +", length(bad_targets)-12, " more"))
    
    targets_by_g <- split(t2_map$target, t2_map$garbage, drop = TRUE) %>% lapply(function(x) sort(unique(x)))
    
    # ---------- standardize DS ----------
    if (!county_var %in% names(ds)) {
        if ("county" %in% names(ds)) ds <- dplyr::rename(ds, !!county_var := .data$county)
        else stop("county_var '", county_var, "' not found (nor 'county').")
    }
    if (!"ucod" %in% names(ds)) stop("ds must contain 'ucod'")
    if (!"year" %in% names(ds) || all(is.na(ds$year))) stop("ds must contain 'year'.")
    
    ds  <- .rename_record_cols(ds)
    rec_cols <- intersect(paste0("record_", 1:20), names(ds))
    
    keep <- c("ucod", county_var, "year", "age_years", "ager27", rec_cols)
    keep <- intersect(keep, names(ds))
    dat  <- ds %>%
        select(all_of(keep)) %>%
        mutate(
            uc4  = icd4(ucod),
            uc3  = icd3(ucod),
            ubin = map_icd_to_bin(ucod, icd_lu)
        )
    
    # ---------- CC long (exclude UCOD; include ALL CC bins) ----------
    cc_long <- if (length(rec_cols) == 0) {
        tibble::tibble(.row = integer(), cc_bin = character())
    } else {
        dat %>%
            mutate(.row = row_number()) %>%
            select(.row, all_of(rec_cols), uc4, uc3) %>%
            pivot_longer(cols = all_of(rec_cols), names_to = "rec", values_to = "cc_raw") %>%
            mutate(cc_icd4 = icd4(cc_raw)) %>%
            filter(!(substr(cc_icd4, 1, 4) == uc4 | substr(cc_icd4, 1, 3) == uc3)) %>%
            transmute(.row, cc_bin = map_icd_to_bin(cc_icd4, icd_lu)) %>%
            filter(!is.na(cc_bin) & cc_bin != "")
    }
    
    row_is_garbage <- !is.na(dat$ubin) & dat$ubin %in% garbage_bins
    
    # ---------- counts ----------
    base_counts_agg <- dat %>%
        filter(!row_is_garbage & !is.na(ubin) & ubin %in% valid_bins) %>%
        count(.data[[county_var]], year, bin = ubin, name = "k", .drop = FALSE)
    counts_cy <- dat %>% count(.data[[county_var]], year, name = "N_total", .drop = FALSE)
    
    # ---------- priors ----------
    prior_by_g <- lapply(garbage_bins, function(g) {
        targs <- targets_by_g[[g]]
        if (is.null(targs) || !length(targs)) return(NULL)
        if (!nrow(cc_long)) {
            setNames(rep(1/length(targs), length(targs)), targs)
        } else {
            rows_with_g_cc_all <- unique(cc_long$.row[cc_long$cc_bin == g])
            tr_all <- intersect(rows_with_g_cc_all, which(!is.na(dat$ubin) & dat$ubin %in% targs))
            if (!length(tr_all)) setNames(rep(1/length(targs), length(targs)), targs) else {
                tbl <- table(dat$ubin[tr_all])
                vec <- rep(0, length(targs)); names(vec) <- targs
                vec[names(tbl)] <- as.numeric(tbl)
                if (sum(vec) == 0) rep(1/length(targs), length(targs)) else (vec / sum(vec))
            }
        }
    })
    names(prior_by_g) <- garbage_bins
    
    # ---------- holders ----------
    prob_counts_list <- list()
    metrics_list     <- list()
    
    # ---------- feature builder ----------
    build_X_for_rows <- function(keep_idx) {
        # CC (all bins)
        if (nrow(cc_long) == 0 || length(keep_idx) == 0) {
            X_cc <- Matrix::Matrix(0, nrow = length(keep_idx), ncol = 0, sparse = TRUE)
            colnames(X_cc) <- character(0)
        } else {
            cc_sub <- cc_long %>% filter(.row %in% keep_idx) %>% distinct(.row, cc_bin)
            if (!nrow(cc_sub)) {
                X_cc <- Matrix::Matrix(0, nrow = length(keep_idx), ncol = 0, sparse = TRUE)
                colnames(X_cc) <- character(0)
            } else {
                feat_bins_local <- sort(unique(cc_sub$cc_bin))
                j_lookup <- setNames(seq_along(feat_bins_local), feat_bins_local)
                i_vec <- match(cc_sub$.row, keep_idx)
                j_vec <- unname(j_lookup[cc_sub$cc_bin])
                ok <- !is.na(i_vec) & !is.na(j_vec)
                i_vec <- i_vec[ok]; j_vec <- j_vec[ok]
                if (length(i_vec)) {
                    key <- paste0(i_vec, "_", j_vec); keepu <- !duplicated(key)
                    i_vec <- i_vec[keepu]; j_vec <- j_vec[keepu]
                }
                X_cc <- Matrix::sparseMatrix(i=i_vec, j=j_vec, x=1,
                                             dims=c(length(keep_idx), length(feat_bins_local)),
                                             dimnames=list(NULL, paste0("cc_", feat_bins_local)))
            }
        }
        # Age: bucket age_years else one-hot ager27
        age_mm <- Matrix::Matrix(0, nrow = length(keep_idx), ncol = 0, sparse = TRUE)
        if ("age_years" %in% names(dat) && any(!is.na(dat$age_years))) {
            age_keep <- suppressWarnings(as.numeric(dat$age_years[keep_idx]))
            age_bucket <- cut(age_keep, breaks = age_breaks, right = FALSE, include.lowest = TRUE)
            levs <- levels(age_bucket); P <- length(levs)
            if (P > 0) {
                ok <- !is.na(age_bucket); i <- which(ok); j <- match(age_bucket[ok], levs)
                age_mm <- Matrix::sparseMatrix(i=i,j=j,x=1,
                                               dims=c(length(keep_idx), P),
                                               dimnames=list(NULL, paste0("age_", levs)))
            }
        } else if ("ager27" %in% names(dat) && any(!is.na(dat$ager27))) {
            a <- suppressWarnings(as.integer(dat$ager27[keep_idx])); ok <- !is.na(a) & a>=1 & a<=27
            if (any(ok)) {
                i <- which(ok); j <- a[ok]
                age_mm <- Matrix::sparseMatrix(i=i,j=j,x=1,
                                               dims=c(length(keep_idx), 27),
                                               dimnames=list(NULL, paste0("ager27_",1:27)))
            }
        }
        if (ncol(age_mm) > 0 && ncol(X_cc) > 0) cbind(X_cc, age_mm)
        else if (ncol(age_mm) > 0) age_mm else X_cc
    }
    
    # ============================ main loop per garbage bin ============================
    for (g in garbage_bins) {
        targets_raw <- targets_by_g[[g]]
        if (is.null(targets_raw) || !length(targets_raw)) next
        
        classes_all <- sort(intersect(targets_raw, valid_bins))
        if (!length(classes_all)) next
        
        rows_with_g_cc <- if (nrow(cc_long)) unique(cc_long$.row[cc_long$cc_bin == g]) else integer(0)
        
        tr_all <- which(!is.na(dat$ubin) & dat$ubin %in% classes_all)
        if (pool_y > 0) {
            te_years <- unique(dat$year[row_is_garbage & dat$ubin == g])
            if (length(te_years)) {
                tr_mask <- dat$year >= min(te_years) - pool_y & dat$year <= max(te_years) + pool_y
                tr_all <- intersect(tr_all, which(tr_mask))
            }
        }
        
        tr_idx <- intersect(rows_with_g_cc, tr_all)
        te_idx <- intersect(which(row_is_garbage), which(dat$ubin == g))
        if (!length(te_idx)) next
        
        # ---------------------- fit multinomial (grouped softmax) -----------------------
        model_ok <- FALSE; probs <- NULL
        n_tr <- length(tr_idx); n_te <- length(te_idx)
        
        if (.has_glmnet() && n_tr >= 30) {
            cls_counts <- table(dat$ubin[tr_idx])
            ok_classes <- names(cls_counts)[cls_counts >= minpc]
            if (length(ok_classes) >= 2) {
                take <- unlist(lapply(split(tr_idx, dat$ubin[tr_idx]),
                                      function(ix) if (length(ix) > capc) sample(ix, capc) else ix),
                               use.names = FALSE)
                tr_idx2 <- sort(take)
                keep_idx <- sort(unique(c(tr_idx2, te_idx)))
                X_local  <- build_X_for_rows(keep_idx)
                nz <- if (ncol(X_local)) which(Matrix::colSums(X_local[match(tr_idx2, keep_idx), , drop = FALSE]) > 0) else integer(0)
                
                if (length(nz)) {
                    X_tr <- X_local[match(tr_idx2, keep_idx), nz, drop = FALSE]
                    X_te <- X_local[match(te_idx,  keep_idx), nz, drop = FALSE]
                    y_tr <- factor(dat$ubin[tr_idx2], levels = sort(unique(ok_classes)))
                    
                    fit <- suppressWarnings(tryCatch(
                        glmnet::cv.glmnet(X_tr, y_tr, family = "multinomial", alpha = 0,
                                          type.multinomial = "grouped", nfolds = 3, parallel = FALSE),
                        error = function(e) NULL
                    ))
                    if (!is.null(fit)) {
                        # 1) extract robustly (no accidental flattening)
                        M_raw <- .extract_multinomial_probs(fit, X_te, classes_all, s = "lambda.1se")
                        # 2) enforce shape and diagnose
                        shp <- .fix_prob_shape(M_raw, n_te = n_te, classes = classes_all,
                                               context = paste0("g=", g, ", dims_in=", paste(dim(M_raw), collapse="×")))
                        probs <- shp$M
                        note  <- shp$note
                        
                        # diagnostics
                        rs0 <- rowSums(probs); dev <- abs(rs0 - 1)
                        frac_bad <- mean(dev > 1e-6)
                        if (verbose) {
                            message(sprintf("[diag] g=%s n_tr=%d n_te=%d k=%d  dims_in=%s  rs_min=%.6f rs_max=%.6f frac_bad=%.3f %s",
                                            g, n_tr, n_te, length(classes_all),
                                            paste(dim(M_raw), collapse="×"),
                                            min(rs0), max(rs0), frac_bad,
                                            if (nzchar(note)) paste0(" (", note, ")") else ""))
                        }
                        diag_rows[[length(diag_rows)+1L]] <- tibble::tibble(
                            garbage = g, n_tr = n_tr, n_te = n_te, k_cand = length(classes_all),
                            dims_in = paste(dim(M_raw), collapse="×"),
                            rs_min = min(rs0), rs_max = max(rs0),
                            frac_rows_ne1e6 = frac_bad,
                            nz_feat = if (exists("nz")) length(nz) else NA_integer_,
                            glmnet_grouped = TRUE, lambda_1se = tryCatch(fit$lambda.1se, error=function(e) NA_real_),
                            reshape_note = note
                        )
                        
                        # 3) final renormalization (harmless if already ≈1)
                        rs0[rs0 == 0] <- 1
                        probs <- probs / rs0
                        model_ok <- TRUE
                    }
                }
            }
        }
        
        # ---------------------- Dirichlet fallback -----------------------
        if (!model_ok) {
            counts <- as.numeric(table(factor(dat$ubin[tr_idx], levels = classes_all)))
            prior_vec <- prior_by_g[[g]]
            prior_aligned <- prior_vec[match(classes_all, names(prior_vec))]
            if (!length(prior_aligned) || any(is.na(prior_aligned))) prior_aligned <- rep(1/length(classes_all), length(classes_all))
            prior_aligned <- prior_aligned / sum(prior_aligned)
            probs <- matrix(rep((counts + alpha0 * prior_aligned) / (sum(counts) + alpha0),
                                each = length(te_idx)),
                            nrow = length(te_idx), ncol = length(classes_all),
                            dimnames = list(NULL, classes_all))
            diag_rows[[length(diag_rows)+1L]] <- tibble::tibble(
                garbage = g, n_tr = n_tr, n_te = n_te, k_cand = length(classes_all),
                dims_in = sprintf("%d×%d", nrow(probs), ncol(probs)),
                rs_min = 1, rs_max = 1, frac_rows_ne1e6 = 0,
                nz_feat = if (exists("nz")) length(nz) else NA_integer_,
                glmnet_grouped = FALSE, lambda_1se = NA_real_, reshape_note = "dirichlet_fallback"
            )
        }
        
        # --------- HARD ASSERT before aggregation: dimensions must match ----------
        if (!is.matrix(probs) || nrow(probs) != length(te_idx) || ncol(probs) != length(classes_all)) {
            stop(sprintf("Bad probs shape: got [%s×%s], expected [%s×%s] (g=%s).",
                         nrow(probs), ncol(probs), length(te_idx), length(classes_all), g))
        }
        if (any(!is.finite(probs))) warning("Non-finite probabilities detected; clamping.")
        probs[!is.finite(probs)] <- 0
        
        # ---- aggregate expected counts for this group's garbage rows ----
        te_meta <- dat[te_idx, c(county_var, "year")]
        colnames(te_meta) <- c("cnty", "year")
        cls_summaries <- lapply(seq_along(classes_all), function(j) {
            tibble::tibble(cnty = te_meta$cnty, year = te_meta$year, k = probs[, j]) %>%
                group_by(cnty, year) %>% summarise(k = sum(k), .groups = "drop") %>%
                mutate(bin = classes_all[j])
        })
        prob_counts_list[[length(prob_counts_list)+1L]] <- dplyr::bind_rows(cls_summaries)
        
        # ---- per-record metrics ----
        eps <- .Machine$double.eps
        H_post <- rowSums(-probs * log(pmax(probs, eps)))
        k_cand <- length(classes_all)
        H0     <- log(k_cand)
        
        Hnorm_gc_by_kcand <- if (k_cand > 1) H_post / log(k_cand) else rep(0, length(H_post))
        Hnorm_post_K      <- if (K > 0) H_post / log(K) else rep(NA_real_, length(H_post))
        expH_over_K_vec   <- if (K > 0) exp(H_post) / K else rep(NA_real_, length(H_post))
        IG_abs_vec        <- H0 - H_post
        IG_frac_vec       <- if (H0 > 0) (H0 - H_post) / H0 else rep(NA_real_, length(H_post))
        
        metr <- tibble::tibble(
            cnty = te_meta$cnty, year = te_meta$year,
            sum_Hnorm_post_K      = Hnorm_post_K,
            sum_expH_over_K       = expH_over_K_vec,
            sum_IG_abs            = IG_abs_vec,
            sum_IG_frac           = IG_frac_vec,
            sum_Hnorm_gc_by_kcand = Hnorm_gc_by_kcand,
            one = 1L
        ) %>%
            group_by(cnty, year) %>%
            summarise(
                sum_Hnorm_post_K      = sum(sum_Hnorm_post_K,      na.rm = TRUE),
                sum_expH_over_K       = sum(sum_expH_over_K,       na.rm = TRUE),
                sum_IG_abs            = sum(sum_IG_abs,            na.rm = TRUE),
                sum_IG_frac           = sum(sum_IG_frac,           na.rm = TRUE),
                sum_Hnorm_gc_by_kcand = sum(sum_Hnorm_gc_by_kcand, na.rm = TRUE),
                N_garbage             = sum(one),
                .groups = "drop"
            )
        metrics_list[[length(metrics_list)+1L]] <- metr
    } # end g loop
    
    prob_counts_all <- if (length(prob_counts_list)) dplyr::bind_rows(prob_counts_list) else
        tibble::tibble(cnty = character(), year = numeric(), k = numeric(), bin = character())
    base_counts_all <- base_counts_agg %>% dplyr::rename(cnty = dplyr::all_of(county_var))
    
    cty_year_bin <- dplyr::bind_rows(
        base_counts_all %>% dplyr::select(cnty, year, bin, k),
        prob_counts_all  %>% dplyr::select(cnty, year, bin, k)
    ) %>%
        dplyr::group_by(cnty, year, bin) %>%
        dplyr::summarise(k = sum(k), .groups = "drop")
    
    out_aggregate <- cty_year_bin %>%
        dplyr::group_by(cnty, year) %>%
        dplyr::summarise(
            total          = sum(k),
            DQ_entropy     = { p <- k / sum(k); .shannon_H(p) },
            DQ_K           = K,
            DQ_overall     = ifelse(DQ_K > 0, DQ_entropy / log(DQ_K), NA_real_),
            DQ_expH_over_K = ifelse(DQ_K > 0, exp(DQ_entropy) / DQ_K, NA_real_),
            .groups = "drop"
        ) %>%
        mutate(
            DQ_overall     = pmax(0, pmin(1, DQ_overall)),
            DQ_expH_over_K = pmax(0, pmin(1, DQ_expH_over_K))
        )
    
    metrics_sum <- if (length(metrics_list)) {
        dplyr::bind_rows(metrics_list) %>%
            dplyr::group_by(cnty, year) %>%
            dplyr::summarise(
                sum_Hnorm_post_K      = sum(sum_Hnorm_post_K,      na.rm = TRUE),
                sum_expH_over_K       = sum(sum_expH_over_K,       na.rm = TRUE),
                sum_IG_abs            = sum(sum_IG_abs,            na.rm = TRUE),
                sum_IG_frac           = sum(sum_IG_frac,           na.rm = TRUE),
                sum_Hnorm_gc_by_kcand = sum(sum_Hnorm_gc_by_kcand, na.rm = TRUE),
                N_garbage             = sum(N_garbage),
                .groups = "drop"
            )
    } else {
        tibble::tibble(cnty = character(), year = numeric(),
                       sum_Hnorm_post_K = numeric(), sum_expH_over_K = numeric(),
                       sum_IG_abs = numeric(), sum_IG_frac = numeric(),
                       sum_Hnorm_gc_by_kcand = numeric(), N_garbage = integer())
    }
    
    counts_cy2 <- counts_cy %>% dplyr::rename(cnty = dplyr::all_of(county_var))
    
    out_perrecord <- counts_cy2 %>%
        dplyr::left_join(metrics_sum,  by = c("cnty", "year")) %>%
        dplyr::mutate(
            across(c(sum_Hnorm_post_K, sum_expH_over_K, sum_IG_abs, sum_IG_frac, sum_Hnorm_gc_by_kcand), ~ dplyr::coalesce(.x, 0)),
            N_garbage = dplyr::coalesce(N_garbage, 0L),
            
            RI = dplyr::if_else(N_garbage > 0, 1 - (sum_Hnorm_gc_by_kcand / N_garbage), NA_real_),
            
            DQ_rec_expH_over_K_mean      = dplyr::if_else(N_garbage > 0, sum_expH_over_K / N_garbage, NA_real_),
            DQ_rec_ig_abs_mean           = dplyr::if_else(N_garbage > 0, sum_IG_abs      / N_garbage, NA_real_),
            DQ_rec_ig_frac_mean_garbage  = dplyr::if_else(N_garbage > 0, sum_IG_frac     / N_garbage, NA_real_)
        ) %>%
        dplyr::select(cnty, year, RI, DQ_rec_expH_over_K_mean, DQ_rec_ig_abs_mean, DQ_rec_ig_frac_mean_garbage)
    
    fg_share <- dat %>%
        dplyr::mutate(is_fg = !is.na(ubin) & ubin %in% garbage_bins) %>%
        dplyr::count(.data[[county_var]], year, is_fg, name = "k") %>%
        tidyr::pivot_wider(names_from = is_fg, values_from = k, values_fill = 0) %>%
        dplyr::transmute(
            cnty = .data[[county_var]], year,
            N_total = `FALSE` + `TRUE`,
            foreman_garbage = dplyr::if_else(N_total > 0, `TRUE` / N_total, NA_real_)
        )
    
    out <- out_aggregate %>%
        dplyr::left_join(out_perrecord, by = c("cnty", "year")) %>%
        dplyr::left_join(fg_share %>% dplyr::select(cnty, year, foreman_garbage),
                         by = c("cnty", "year")) %>%
        dplyr::mutate(
            foreman_garbage_adj = foreman_garbage * (1 - dplyr::coalesce(DQ_rec_ig_frac_mean_garbage, 0)),
            RI = ifelse(is.na(RI), RI, pmax(0, pmin(1, RI)))
        ) %>%
        dplyr::rename(!!county_var := cnty)
    
    # -------------------------- write diagnostics if requested ---------------------------
    if (length(diag_rows)) {
        diag_df <- dplyr::bind_rows(diag_rows)
        if (!is.null(diag_dir)) {
            if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
            ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
            f <- file.path(diag_dir, paste0("glmnet_prob_diag_", ts, ".csv"))
            readr::write_csv(diag_df, f)
            if (verbose) message("[diag] wrote probability diagnostics to: ", f)
        }
        attr(out, "prob_diag") <- diag_df
    }
    
    # final sanity notes
    if (any(out$DQ_overall > 1 + 1e-9, na.rm = TRUE) || any(out$DQ_overall < -1e-9, na.rm = TRUE))
        warning("DQ_overall outside [0,1]; values were clamped.")
    if (any(out$DQ_expH_over_K > 1 + 1e-9, na.rm = TRUE))
        warning("DQ_expH_over_K exceeded 1; values were clamped.")
    
    out
}
