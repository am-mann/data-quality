# ───────────────── dq_entropy_helper.R — Foreman bins + RI (county×year) ─────────────────
#
# Outputs (county×year):
#   total, DQ_overall, RI, RI_post_only, N_garbage
# ──────────────────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(readr); library(stringr)
    library(purrr); library(Matrix)
})

# ------- Options & helpers -------
.pool_years_default       <- function() as.integer(getOption("DQ_POOL_YEARS", 2L))
.dirichlet_prior_default  <- function() as.numeric(getOption("DQ_DIRICHLET_PRIOR", 0.25))
.min_per_class_default    <- function() as.integer(getOption("DQ_MIN_PER_CLASS", 5L))
.max_train_per_class_def  <- function() as.integer(getOption("DQ_MAX_TRAIN_PER_CLASS", 3000L))
.seed_default             <- function() as.integer(getOption("DQ_SEED", NA_integer_))
.verbose_default          <- function() isTRUE(getOption("DQ_VERBOSE", TRUE))
.diag_dir_default         <- function() { x <- getOption("DQ_DIAG_DIR", NULL); if (is.character(x) && nzchar(x)) x else NULL }
.has_glmnet               <- function() requireNamespace("glmnet", quietly = TRUE)
.cv_folds_default         <- function() as.integer(getOption("DQ_CV_FOLDS", 3L))
.nlambda_default          <- function() as.integer(getOption("DQ_NLAMBDA", 60L))
.lmr_default              <- function() as.numeric(getOption("DQ_LAMBDA_MIN_RATIO", 0.05))

# NEW knobs (safe defaults)
.alpha_default            <- function() as.numeric(getOption("DQ_ALPHA", 0.7))           # sparsity 0..1
.lambda_choice_default    <- function() as.character(getOption("DQ_LAMBDA_CHOICE","lambda.min")) # "lambda.min"|"lambda.1se"
.tau_default              <- function() as.numeric(getOption("DQ_POST_TAU", 0.75))       # <1 sharpens metrics-only
.min_train_support_def    <- function() as.integer(getOption("DQ_MIN_TRAIN_SUPPORT", 20))# prune low-support targets
.max_k_cand_def           <- function() { x <- getOption("DQ_MAX_K_CAND", NA_integer_); if (is.na(x)) NA_integer_ else as.integer(x) }
.target_per_class_def     <- function() as.integer(getOption("DQ_TARGET_PER_CLASS", 1500)) # balance training

`%||%`        <- function(a,b) if (!is.null(a)) a else b
canonical_icd <- function(x) stringr::str_remove_all(stringr::str_to_upper(x), "[^A-Z0-9]")
icd3          <- function(x) substr(canonical_icd(x), 1, 3)
icd4          <- function(x) substr(canonical_icd(x), 1, 4)
norm_code     <- function(x) tolower(gsub("[^A-Za-z0-9]+", "", as.character(x)))
.shannon_H    <- function(p) { p <- p[p > 0]; if (!length(p)) 0 else -sum(p * log(p)) }

.rename_record_cols <- function(df) {
    if (!any(grepl("^record_", names(df))) && any(grepl("^record[0-9]+$", names(df)))) {
        dplyr::rename_with(df, ~ sub("^record([0-9]+)$","record_\\1", .x), .cols = dplyr::matches("^record[0-9]+$"))
    } else df
}

.read_icd_map <- function(path) {
    df <- readr::read_csv(path, show_col_types = FALSE,
                          col_types = readr::cols(.default = readr::col_character(),
                                                  ICD10 = readr::col_character(),
                                                  USCOD = readr::col_character()))
    if (!all(c("ICD10","USCOD") %in% names(df))) stop("foreman-icd10-mapping.csv must have ICD10 and USCOD.")
    out <- df %>% transmute(icd = canonical_icd(ICD10), bin = norm_code(USCOD)) %>%
        filter(icd != "", bin != "") %>% distinct(icd, .keep_all = TRUE)
    if (any(grepl("^\\d+$", out$bin))) stop("USCOD parsed as digits-only; ensure alphanumeric codes.")
    out
}

.read_table2_map <- function(path) {
    df <- readr::read_csv(path, show_col_types = FALSE, col_types = readr::cols(.default = readr::col_character()))
    if (!"target_cause" %in% names(df)) stop("foreman-table2-map.csv must contain 'target_cause'.")
    g_cols <- grep("^G_[0-9]+$", names(df), value = TRUE); if (!length(g_cols)) stop("Need columns G_1..G_9.")
    df[g_cols] <- lapply(df[g_cols], function(x){ y <- tolower(as.character(x)); y %in% c("1","true","t","y","yes","x","✓","check","checked") })
    t2 <- df %>% mutate(target = norm_code(target_cause)) %>%
        pivot_longer(all_of(g_cols), names_to = "garbage_raw", values_to = "flag") %>% filter(flag) %>%
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

.align_cols <- function(M, classes) {
    have <- colnames(M) %||% character(0)
    keep <- intersect(classes, have); miss <- setdiff(classes, have)
    base <- if (length(keep)) M[, keep, drop = FALSE] else matrix(, nrow = nrow(M), ncol = 0)
    if (length(miss)) base <- cbind(base, matrix(0, nrow = nrow(M), ncol = length(miss), dimnames = list(NULL, miss)))
    base[, classes, drop = FALSE]
}

.extract_multinomial_probs <- function(fit, X_te, classes, s = "lambda.1se") {
    out <- tryCatch(predict(fit, newx = X_te, type = "response", s = s), error = function(e) NULL)
    if (is.null(out)) return(NULL)
    if (is.array(out) && length(dim(out)) == 3L) {
        M <- out[, , 1, drop = TRUE]; if (!is.matrix(M)) M <- as.matrix(M)
        if (is.null(colnames(M)) || anyNA(colnames(M))) {
            cls <- tryCatch(fit$glmnet.fit$classnames, error = function(e) NULL)
            if (!is.null(cls) && length(cls) == ncol(M)) colnames(M) <- cls
        }
        return(.align_cols(M, classes))
    }
    if (is.matrix(out) && length(dim(out)) == 2L) {
        if (is.null(colnames(out)) || anyNA(colnames(out))) {
            cls <- tryCatch(fit$glmnet.fit$classnames, error = function(e) NULL)
            if (!is.null(cls) && length(cls) == ncol(out)) colnames(out) <- cls
        }
        return(.align_cols(out, classes))
    }
    if (is.list(out)) {
        mats <- lapply(out, function(x) if (is.matrix(x)) x[, 1, drop = TRUE] else as.numeric(x))
        M <- do.call(cbind, mats); colnames(M) <- names(out) %||% paste0("cls", seq_len(ncol(M)))
        return(.align_cols(as.matrix(M), classes))
    }
    .align_cols(as.matrix(out), classes)
}

.fix_prob_shape <- function(M, n_te, classes, context = NULL) {
    k <- length(classes)
    if (is.matrix(M) && nrow(M) == n_te && ncol(M) == k) {
        if (is.null(colnames(M)) || anyNA(colnames(M))) colnames(M) <- classes
        return(list(M=M, note="ok"))
    }
    dims <- dim(M)
    if (is.matrix(M) && all(dims == c(k, n_te))) {
        M2 <- t(M); colnames(M2) <- classes
        return(list(M=M2, note="transpose"))
    }
    if (is.matrix(M) && nrow(M) == n_te * k && ncol(M) == k) {
        grp <- rep(seq_len(n_te), each = k)
        M2  <- rowsum(M, grp)
        rs  <- rowSums(M2); rs[rs == 0] <- 1
        M2  <- M2 / rs
        colnames(M2) <- classes
        return(list(M = M2, note = "collapse[(n_te*k)×k] via rowsum"))
    }
    if (is.null(dims) && length(M) == n_te * k) {
        M2 <- matrix(as.numeric(M), nrow = n_te, ncol = k, byrow = TRUE); colnames(M2) <- classes
        return(list(M=M2, note="vector reshape"))
    }
    stop("Bad prob shape: ", paste(dims, collapse="×"), if (!is.null(context)) paste0(" (", context, ")") else "")
}

# ============================== MAIN FUNCTION ==============================
compute_entropy_county_foreman <- function(
        ds, county_var,
        dict_dir      = NULL,
        icd_map_path  = if (!is.null(dict_dir)) file.path(dict_dir, "foreman-icd10-mapping.csv") else NULL,
        code_map_path = if (!is.null(dict_dir)) file.path(dict_dir, "foreman-table2-map.csv") else NULL,
        age_breaks    = c(-Inf, 1, 5, seq(10, 85, by = 5), Inf)
) {
    if (is.null(icd_map_path) || is.null(code_map_path)) stop("Provide icd_map_path and code_map_path (or dict_dir).")
    if (!.has_glmnet()) stop("glmnet package is required.")
    
    seed     <- .seed_default(); if (!is.na(seed)) set.seed(seed)
    alpha0   <- .dirichlet_prior_default()
    minpc    <- .min_per_class_default()
    capc     <- .max_train_per_class_def()
    pool_y   <- .pool_years_default()
    verbose  <- .verbose_default()
    diag_dir <- .diag_dir_default()
    diag_rows <- list(); on.exit({gc()}, add = TRUE)
    
    if (verbose) {
        message(sprintf("[settings] cap/class=%d, CV folds=%d, nlambda=%d, pool years=%d, alpha=%.2f, lambda=%s, tau=%.2f",
                        capc, .cv_folds_default(), .nlambda_default(), pool_y,
                        .alpha_default(), .lambda_choice_default(), .tau_default()))
    }
    
    icd_map  <- .read_icd_map(icd_map_path)
    t2_map   <- .read_table2_map(code_map_path)
    icd_lu   <- .make_icd_lookup(icd_map)
    
    garbage_bins     <- sort(unique(t2_map$garbage))
    all_bins_in_map  <- sort(unique(icd_map$bin))
    valid_bins       <- setdiff(all_bins_in_map, garbage_bins)
    K <- length(valid_bins)
    if (K < 2) stop("Too few valid bins found in dictionaries.")
    
    targets_by_g <- split(t2_map$target, t2_map$garbage, drop = TRUE) |>
        lapply(function(x) sort(unique(x)))
    
    if (!county_var %in% names(ds)) {
        if ("county" %in% names(ds)) ds <- dplyr::rename(ds, !!county_var := .data$county)
        else stop("county_var '", county_var, "' not found.")
    }
    if (!"ucod" %in% names(ds)) stop("ds must contain 'ucod'")
    if (!"year" %in% names(ds) || all(is.na(ds$year))) stop("ds must contain 'year'.")
    
    ds  <- .rename_record_cols(ds)
    rec_cols <- intersect(paste0("record_", 1:20), names(ds))
    keep <- c("ucod", county_var, "year", "age_years", "ager27", "sex", rec_cols)
    keep <- intersect(keep, names(ds))
    
    dat <- ds %>%
        dplyr::select(dplyr::all_of(keep)) %>%
        dplyr::mutate(
            uc4  = icd4(ucod),
            uc3  = icd3(ucod),
            ubin = map_icd_to_bin(ucod, icd_lu)
        )
    
    row_is_garbage <- !is.na(dat$ubin) & dat$ubin %in% garbage_bins
    
    base_counts_agg <- dat %>%
        dplyr::filter(!is.na(ubin) & ubin %in% valid_bins) %>%
        dplyr::count(.data[[county_var]], year, bin = ubin, name = "k", .drop = FALSE)
    base_counts_all <- base_counts_agg %>% dplyr::rename(cnty = dplyr::all_of(county_var))
    
    # age bucket + sex keys (for prior only)
    age_bucket_vec <- dplyr::case_when(
        "age_years" %in% names(dat) & !is.na(dat$age_years) ~ as.character(cut(
            suppressWarnings(as.numeric(dat$age_years)),
            breaks = age_breaks, right = FALSE, include.lowest = TRUE
        )),
        "ager27" %in% names(dat) & !is.na(dat$ager27) ~ paste0("ager27_", pmax(1L, pmin(27L, suppressWarnings(as.integer(dat$ager27))))),
        TRUE ~ NA_character_
    )
    sex_vec <- dplyr::case_when(
        "sex" %in% names(dat) & dat$sex %in% c("M","F") ~ dat$sex,
        TRUE ~ NA_character_
    )
    
    # Sparse feature builder: record_* counts (+ optional age dummies)
    build_X_for_rows <- function(keep_idx) {
        rec_cols_local <- intersect(paste0("record_", 1:20), names(dat))
        n_keep <- length(keep_idx)
        
        # contributing cause (USCOD-bin) *COUNTS* (no 0/1 override)
        if (n_keep == 0 || length(rec_cols_local) == 0) {
            X_cc <- Matrix(0, nrow = n_keep, ncol = 0, sparse = TRUE)
            colnames(X_cc) <- character(0)
        } else {
            bin_to_j <- new.env(parent = emptyenv())
            j_to_bin <- character(0)
            i_vec <- integer(0)
            j_vec <- integer(0)
            
            for (col in rec_cols_local) {
                v <- dat[[col]][keep_idx]
                v <- canonical_icd(v)
                b <- map_icd_to_bin(v, icd_lu)
                ok <- !is.na(b) & nzchar(b)
                if (!any(ok)) next
                
                b_sub <- b[ok]
                new_bins <- unique(b_sub[!(b_sub %in% j_to_bin)])
                if (length(new_bins)) {
                    start <- length(j_to_bin)
                    j_to_bin <- c(j_to_bin, new_bins)
                    for (k2 in seq_along(new_bins)) assign(new_bins[k2], start + k2, envir = bin_to_j)
                }
                j_sub <- unname(vapply(b_sub, function(x) get(x, envir = bin_to_j, inherits = FALSE), 1L))
                i_sub <- which(ok)  # duplicates across columns accumulate -> counts
                
                i_vec <- c(i_vec, i_sub)
                j_vec <- c(j_vec, j_sub)
                
                rm(v, b, ok, b_sub, j_sub, i_sub); gc(FALSE)
            }
            
            if (length(i_vec)) {
                X_cc <- sparseMatrix(i = i_vec, j = j_vec, x = 1,
                                     dims = c(n_keep, length(j_to_bin)),
                                     dimnames = list(NULL, paste0("cc_", j_to_bin)))
                # IMPORTANT: do NOT force X_cc@x <- 1; duplicates already sum to counts
                # Matrix::sparseMatrix sums duplicates by default.
            } else {
                X_cc <- Matrix(0, nrow = n_keep, ncol = 0, sparse = TRUE)
                colnames(X_cc) <- character(0)
            }
        }
        
        # Optional age dummies (help classification; not used in prior construction)
        age_mm <- Matrix(0, nrow = n_keep, ncol = 0, sparse = TRUE)
        if ("age_years" %in% names(dat) && any(!is.na(dat$age_years[keep_idx]))) {
            age_keep   <- suppressWarnings(as.numeric(dat$age_years[keep_idx]))
            age_bucket <- cut(age_keep, breaks = age_breaks, right = FALSE, include.lowest = TRUE)
            levs <- levels(age_bucket); P <- length(levs)
            if (P > 0) {
                ok <- !is.na(age_bucket); i <- which(ok); j <- match(age_bucket[ok], levs)
                age_mm <- sparseMatrix(i = i, j = j, x = 1, dims = c(n_keep, P), dimnames = list(NULL, paste0("age_", levs)))
            }
        } else if ("ager27" %in% names(dat) && any(!is.na(dat$ager27[keep_idx]))) {
            a <- suppressWarnings(as.integer(dat$ager27[keep_idx])); ok <- !is.na(a) & a >= 1 & a <= 27
            if (any(ok)) {
                i <- which(ok); j <- a[ok]
                age_mm <- sparseMatrix(i = i, j = j, x = 1, dims = c(n_keep, 27), dimnames = list(NULL, paste0("ager27_", 1:27)))
            }
        }
        
        if (ncol(age_mm) > 0L && ncol(X_cc) > 0L) cbind(X_cc, age_mm)
        else if (ncol(age_mm) > 0L) age_mm
        else X_cc
    }
    
    # Prior constructor from TRAINING rows by age×sex (Dirichlet-smoothed)
    .make_age_sex_prior_getter <- function(classes_all, tr_idx2, age_bucket_vec, sex_vec, alpha0) {
        cls_train <- as.character(dat$ubin[tr_idx2])
        age_tr    <- age_bucket_vec[tr_idx2]
        sex_tr    <- sex_vec[tr_idx2]
        
        df <- tibble::tibble(
            cls = cls_train,
            age = age_tr,
            sex = sex_tr
        ) %>%
            dplyr::filter(!is.na(cls), cls %in% classes_all) %>%
            dplyr::mutate(
                age = ifelse(is.na(age), "_NA_", age),
                sex = ifelse(is.na(sex), "_NA_", sex)
            )
        
        tabs <- list(
            by_age_sex = df %>% dplyr::count(cls, age, sex, name = "k"),
            by_age     = df %>% dplyr::count(cls, age, name = "k"),
            by_sex     = df %>% dplyr::count(cls, sex, name = "k"),
            by_all     = df %>% dplyr::count(cls, name = "k")
        )
        
        smoothed_norm <- function(named_counts) {
            v <- setNames(numeric(length(classes_all)), classes_all)
            if (!is.null(named_counts) && length(named_counts)) v[names(named_counts)] <- named_counts
            v <- v + alpha0
            v / sum(v)
        }
        
        function(age_key, sex_key) {
            a <- ifelse(is.na(age_key), "_NA_", age_key)
            s <- ifelse(is.na(sex_key), "_NA_", sex_key)
            
            vec <- tabs$by_age_sex %>% dplyr::filter(age == a, sex == s)
            if (nrow(vec)) return(smoothed_norm(stats::setNames(vec$k, vec$cls)))
            
            vec <- tabs$by_age %>% dplyr::filter(age == a)
            if (nrow(vec)) return(smoothed_norm(stats::setNames(vec$k, vec$cls)))
            
            vec <- tabs$by_sex %>% dplyr::filter(sex == s)
            if (nrow(vec)) return(smoothed_norm(stats::setNames(vec$k, vec$cls)))
            
            if (nrow(tabs$by_all)) return(smoothed_norm(stats::setNames(tabs$by_all$k, tabs$by_all$cls)))
            
            rep(1/length(classes_all), length(classes_all))  # uniform fallback
        }
    }
    
    prob_counts_list <- list()
    metrics_list     <- list()
    
    # ============================ per garbage bin ============================
    for (g in garbage_bins) {
        targets_raw <- targets_by_g[[g]]; if (is.null(targets_raw) || !length(targets_raw)) next
        classes_all <- sort(intersect(targets_raw, valid_bins)); if (!length(classes_all)) next
        
        tr_all <- which(!is.na(dat$ubin) & dat$ubin %in% classes_all)
        
        # Temporal pooling of training around test years in this bin
        if (pool_y > 0) {
            te_years_probe <- unique(dat$year[row_is_garbage & dat$ubin == g])
            if (length(te_years_probe)) {
                tr_mask <- dat$year >= min(te_years_probe) - pool_y & dat$year <= max(te_years_probe) + pool_y
                tr_all <- intersect(tr_all, which(tr_mask))
            }
        }
        
        te_idx <- intersect(which(row_is_garbage), which(dat$ubin == g))
        if (!length(te_idx)) next
        
        n_tr <- length(tr_all); n_te <- length(te_idx)
        if (n_tr < 30) stop(sprintf("No-fallback: g=%s too few training rows (n_tr=%d < 30).", g, n_tr))
        
        # Compute per-class counts on pooled training
        cls_counts_all <- table(dat$ubin[tr_all])
        
        # Minimum per-class present for model viability
        ok_classes <- names(cls_counts_all)[cls_counts_all >= minpc]
        if (length(ok_classes) < 2) stop(sprintf("No-fallback: g=%s has <2 classes with >=%d rows.", g, minpc))
        
        # --- NEW: prune by minimum pooled support + optional top-k cap
        min_support <- .min_train_support_def()
        keep_cands  <- names(cls_counts_all)[cls_counts_all >= min_support]
        ok_classes  <- intersect(ok_classes, keep_cands)
        
        k_cap <- .max_k_cand_def()
        if (!is.na(k_cap) && length(ok_classes) > k_cap) {
            ok_sorted <- names(sort(cls_counts_all[ok_classes], decreasing = TRUE))
            ok_classes <- ok_sorted[seq_len(k_cap)]
        }
        
        if (length(ok_classes) < 2) stop(sprintf("No-fallback: g=%s <2 viable classes after pruning.", g))
        classes_all <- sort(intersect(classes_all, ok_classes))
        
        # --- Build feature matrices on balanced training+test subset
        keep_idx <- sort(unique(c(tr_all, te_idx)))
        X_local  <- build_X_for_rows(keep_idx)
        
        # Drop empty columns on training slice
        nz <- if (ncol(X_local)) which(Matrix::colSums(X_local[match(tr_all, keep_idx), , drop = FALSE]) > 0) else integer(0)
        if (!length(nz)) stop(sprintf("No-fallback: g=%s empty feature matrix after filtering.", g))
        
        # --- NEW: class-balanced training (undersample/oversample toward target)
        target_per_class <- .target_per_class_def()
        tr_by_class <- split(tr_all, dat$ubin[tr_all])
        take <- unlist(lapply(tr_by_class, function(ix) {
            n <- length(ix)
            if (n >= target_per_class) sample(ix, target_per_class)
            else c(ix, sample(ix, target_per_class - n, replace = TRUE))
        }), use.names = FALSE)
        tr_idx2 <- sort(unique(take))
        
        X_tr <- X_local[match(tr_idx2, keep_idx), nz, drop = FALSE]
        X_te <- X_local[match(te_idx,  keep_idx), nz, drop = FALSE]
        rm(X_local); gc(FALSE)
        
        y_tr <- factor(dat$ubin[tr_idx2], levels = sort(unique(ok_classes)))
        
        if (verbose) {
            tbl <- sort(table(y_tr), decreasing = TRUE)
            msg <- paste(head(paste(names(tbl), as.integer(tbl), sep=":"), 6), collapse=", ")
            message(sprintf("[diag] g=%s n_tr=%d n_te=%d k=%d; top classes: %s",
                            g, length(tr_idx2), n_te, length(levels(y_tr)), msg))
        }
        
        # --- NEW: elastic net & lambda choice
        fit <- suppressWarnings(tryCatch(
            glmnet::cv.glmnet(
                X_tr, y_tr,
                family = "multinomial",
                alpha  = .alpha_default(),                 # NEW
                type.multinomial = "grouped",
                nfolds = .cv_folds_default(),
                nlambda = .nlambda_default(),
                lambda.min.ratio = .lmr_default(),
                standardize = FALSE,
                parallel = FALSE
            ),
            error = function(e) NULL
        ))
        if (is.null(fit)) stop(sprintf("No-fallback: glmnet failed for g=%s.", g))
        
        lambda_choice <- .lambda_choice_default()
        M_raw <- .extract_multinomial_probs(fit, X_te, classes_all, s = lambda_choice)
        if (is.null(M_raw)) stop(sprintf("No-fallback: prediction NULL for g=%s.", g))
        
        shp <- .fix_prob_shape(M_raw, n_te = n_te, classes = classes_all,
                               context = paste0("g=", g, ", dims_in=", paste(dim(M_raw), collapse="×")))
        probs <- shp$M
        rs0 <- rowSums(probs); rs0[rs0 == 0] <- 1; probs <- probs / rs0
        
        # --------- age×sex prior from TRAINING rows ----------
        prior_getter <- .make_age_sex_prior_getter(
            classes_all = classes_all,
            tr_idx2     = tr_idx2,
            age_bucket_vec = age_bucket_vec,
            sex_vec        = sex_vec,
            alpha0      = alpha0
        )
        age_te <- age_bucket_vec[te_idx]
        sex_te <- sex_vec[te_idx]
        
        prior_mat <- matrix(NA_real_, nrow = length(te_idx), ncol = length(classes_all))
        for (i in seq_along(te_idx)) prior_mat[i, ] <- prior_getter(age_te[i], sex_te[i])
        colnames(prior_mat) <- classes_all
        
        # ----- Aggregate predicted mass to county×year×bin for DQ_overall -----
        te_meta <- dat[te_idx, c(county_var, "year")]; colnames(te_meta) <- c("cnty", "year")
        cls_summaries <- lapply(seq_along(classes_all), function(j) {
            tibble::tibble(cnty = te_meta$cnty, year = te_meta$year, k = probs[, j]) %>%
                dplyr::group_by(cnty, year) %>% dplyr::summarise(k = sum(k), .groups = "drop") %>%
                dplyr::mutate(bin = classes_all[j])
        })
        prob_counts_list[[length(prob_counts_list)+1L]] <- dplyr::bind_rows(cls_summaries)
        
        # ----- Per-record metrics: KL uses original probs; entropy uses SHARPENED probs -----
        eps <- .Machine$double.eps
        # metrics-only sharpening
        probs_for_metrics <- probs
        tau <- .tau_default()
        if (!is.na(tau) && tau > 0 && tau != 1) {
            probs_for_metrics <- probs_for_metrics ^ tau
            probs_for_metrics <- probs_for_metrics / rowSums(probs_for_metrics)
        }
        
        H_post  <- rowSums(-probs_for_metrics * log(pmax(probs_for_metrics,  eps)))
        k_cand  <- length(classes_all)
        Hnorm_gc_by_kcand <- if (k_cand > 1) H_post / log(k_cand) else rep(0, length(H_post))
        # KL against prior (do NOT sharpen)
        KL_raw_vec  <- rowSums(probs * (log(pmax(probs, eps)) - log(pmax(prior_mat, eps))))
        KL_norm_vec <- if (k_cand > 1) KL_raw_vec / log(k_cand) else rep(0, length(KL_raw_vec))
        
        metr <- tibble::tibble(
            cnty = te_meta$cnty, year = te_meta$year,
            sum_KL_norm           = KL_norm_vec,
            sum_Hnorm_gc_by_kcand = Hnorm_gc_by_kcand,
            one = 1L
        ) %>%
            dplyr::group_by(cnty, year) %>%
            dplyr::summarise(
                sum_KL_norm           = sum(sum_KL_norm,           na.rm = TRUE),
                sum_Hnorm_gc_by_kcand = sum(sum_Hnorm_gc_by_kcand, na.rm = TRUE),
                N_garbage             = sum(one),
                .groups = "drop"
            )
        metrics_list[[length(metrics_list)+1L]] <- metr
    } # end g loop
    
    prob_counts_all <- if (length(prob_counts_list)) dplyr::bind_rows(prob_counts_list) else
        tibble::tibble(cnty = character(), year = numeric(), k = numeric(), bin = character())
    
    # Final county×year×bin distribution = observed valid UCOD + redistributed garbage mass
    cty_year_bin <- dplyr::bind_rows(
        base_counts_all %>% dplyr::select(cnty, year, bin, k),
        prob_counts_all  %>% dplyr::select(cnty, year, bin, k)
    ) %>%
        dplyr::group_by(cnty, year, bin) %>%
        dplyr::summarise(k = sum(k), .groups = "drop")
    
    out_aggregate <- cty_year_bin %>%
        dplyr::group_by(cnty, year) %>%
        dplyr::summarise(
            total      = sum(k),
            DQ_overall = { p <- k / sum(k); .shannon_H(p) / log(K) },
            .groups = "drop"
        ) %>%
        dplyr::mutate(DQ_overall = pmax(0, pmin(1, DQ_overall)))
    
    counts_cy <- dat %>% dplyr::count(.data[[county_var]], year, name = "N_total", .drop = FALSE)
    counts_cy2 <- counts_cy %>% dplyr::rename(cnty = dplyr::all_of(county_var))
    
    metrics_sum <- if (length(metrics_list)) {
        dplyr::bind_rows(metrics_list) %>%
            dplyr::group_by(cnty, year) %>%
            dplyr::summarise(
                sum_KL_norm           = sum(sum_KL_norm,           na.rm = TRUE),
                sum_Hnorm_gc_by_kcand = sum(sum_Hnorm_gc_by_kcand, na.rm = TRUE),
                N_garbage             = sum(N_garbage),
                .groups = "drop"
            )
    } else {
        tibble::tibble(cnty = character(), year = numeric(),
                       sum_KL_norm = numeric(), sum_Hnorm_gc_by_kcand = numeric(),
                       N_garbage = integer())
    }
    
    out_perrecord <- counts_cy2 %>%
        dplyr::left_join(metrics_sum,  by = c("cnty", "year")) %>%
        dplyr::mutate(
            across(c(sum_KL_norm, sum_Hnorm_gc_by_kcand), ~ dplyr::coalesce(.x, 0)),
            N_garbage  = dplyr::coalesce(N_garbage, 0L),
            RI         = dplyr::if_else(N_garbage > 0, sum_KL_norm / N_garbage, NA_real_),
            RI_post_only = dplyr::if_else(N_garbage > 0, 1 - (sum_Hnorm_gc_by_kcand / N_garbage), NA_real_)
        ) %>%
        dplyr::mutate(
            RI           = ifelse(is.na(RI), RI, pmax(0, pmin(1, RI))),
            RI_post_only = ifelse(is.na(RI_post_only), RI_post_only, pmax(0, pmin(1, RI_post_only)))
        ) %>%
        dplyr::select(cnty, year, RI, RI_post_only, N_garbage)
    
    out <- out_aggregate %>%
        dplyr::left_join(out_perrecord, by = c("cnty", "year")) %>%
        dplyr::rename(!!county_var := cnty)
    
    # Optional probability diagnostics hook (leave empty unless you add diag rows)
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
    
    out
}

