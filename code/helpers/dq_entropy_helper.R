# ───────────────── dq_entropy_helper.R — Foreman bins + RI (+ JSD)
#      g3/g8/g9 fully excluded; NO prior-only fallback;
#      CC-presence features, robust quartiles, training backoff, diagnostics
#      Adds pct of certs that are {g1,g2,g4,g5,g6,g7} by county×year
# ─────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(readr); library(stringr)
    library(purrr); library(Matrix); library(tibble); library(jsonlite)
})

# ====== Recommended knobs (edit as needed) ======
# Per-bin strat modes (choices: "age_sex", "age_only", "sex_only", "none")
# No longer reference g3/g8/g9 anywhere, they're excluded
options(DQ_STRAT_MODE_BY_BIN = list())
# Rarity & capacity
options(
    DQ_MIN_PER_CLASS      = 8L,
    DQ_MIN_TRAIN_SUPPORT  = 10L,
    DQ_TARGET_PER_CLASS   = 3000,
    DQ_POOL_YEARS         = 3L
)

# ====== Option helpers ======
.pool_years_default       <- function() as.integer(getOption("DQ_POOL_YEARS", 3L))
.dirichlet_prior_default  <- function() as.numeric(getOption("DQ_DIRICHLET_PRIOR", 0.25))
.min_per_class_default    <- function() as.integer(getOption("DQ_MIN_PER_CLASS", 8L))
.max_train_per_class_def  <- function() as.integer(getOption("DQ_MAX_TRAIN_PER_CLASS", 3000L))
.seed_default             <- function() as.integer(getOption("DQ_SEED", NA_integer_))
.verbose_default          <- function() isTRUE(getOption("DQ_VERBOSE", TRUE))
.diag_dir_default         <- function() { x <- getOption("DQ_DIAG_DIR", NULL); if (is.character(x) && nzchar(x)) x else NULL }
.has_glmnet               <- function() requireNamespace("glmnet", quietly = TRUE)
# +++ NEW: detect doParallel
.has_doParallel           <- function() requireNamespace("doParallel", quietly = TRUE)
.cv_folds_default         <- function() as.integer(getOption("DQ_CV_FOLDS", 3L))
.nlambda_default          <- function() as.integer(getOption("DQ_NLAMBDA", 60L))
.lmr_default              <- function() as.numeric(getOption("DQ_LAMBDA_MIN_RATIO", 0.05))
.alpha_default            <- function() as.numeric(getOption("DQ_ALPHA", 0.7))
.lambda_choice_default    <- function() as.character(getOption("DQ_LAMBDA_CHOICE","lambda.min"))
.tau_default              <- function() as.numeric(getOption("DQ_POST_TAU", 0.75))
.min_train_support_def    <- function() as.integer(getOption("DQ_MIN_TRAIN_SUPPORT", 10L))
.target_per_class_def     <- function() as.integer(getOption("DQ_TARGET_PER_CLASS", 3000L))

# Per-bin stratification mode
.strat_mode_for_bin <- function(g) {
    per <- getOption("DQ_STRAT_MODE_BY_BIN", NULL)
    mode <- if (is.list(per) && !is.null(per[[g]])) as.character(per[[g]]) else "age_sex"
    if (!mode %in% c("age_sex","age_only","sex_only","none")) mode <- "age_sex"
    mode
}

# ====== Utilities ======
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
    base <- if (length(keep)) M[, keep, drop = FALSE] else matrix(0, nrow = nrow(M), ncol = 0)
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
    if (is.matrix(M) && dims[1] == n_te * k && dims[2] == k) {
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
    stop("Bad prob shape: ", paste(dims, collapse="×"),
         " expected ", n_te, "×", k, if (!is.null(context)) paste0(" (", context, ")") else "")
}

# ---------- Sex parsing & stratification guards ----------
.normalize_sex <- function(x) {
    sx <- toupper(trimws(as.character(x)))
    sx[sx %in% c("1","M","MALE")]   <- "M"
    sx[sx %in% c("2","F","FEMALE")] <- "F"
    sx[!(sx %in% c("M","F"))] <- NA
    sx
}
.warn_if_sex_bad <- function(sex2, context = "dataset") {
    n  <- length(sex2); if (!n) return(invisible(NULL))
    n_na <- sum(is.na(sex2)); p_na <- n_na / n
    lev  <- sort(unique(sex2[!is.na(sex2)]))
    if (p_na >= 0.8 || length(lev) == 0) {
        warning(sprintf("[DQ] Sex parsing failed for %s: %d/%d NA (%.1f%%). Stratification by sex will be disabled where needed.",
                        context, n_na, n, 100 * p_na), call. = FALSE)
    } else if (length(lev) == 1) {
        warning(sprintf("[DQ] Sex has a single level ('%s') in %s; treating as effectively missing for stratification.",
                        lev, context), call. = FALSE)
    }
    invisible(NULL)
}
.effective_mode_for_slice <- function(mode, sex_vec, idx) {
    if (!(mode %in% c("age_sex","sex_only"))) return(mode)
    if (length(idx) == 0) return(mode)
    sx <- sex_vec[idx]
    lev <- unique(sx[!is.na(sx)])
    if (length(lev) < 2) {
        if (mode == "age_sex") return("age_only")
        if (mode == "sex_only") return("none")
    }
    mode
}

# Build sparse presence-only CC features for selected rows
.build_X_for_rows <- function(keep_idx, cc_long) {
    if (NROW(cc_long) == 0L || length(keep_idx) == 0L) {
        X_cc <- Matrix(0, nrow = length(keep_idx), ncol = 0, sparse = TRUE)
        colnames(X_cc) <- character(0); return(X_cc)
    }
    cc_sub <- cc_long %>% dplyr::filter(.row %in% keep_idx) %>% dplyr::distinct(.row, cc_bin)
    if (!NROW(cc_sub)) {
        X_cc <- Matrix(0, nrow = length(keep_idx), ncol = 0, sparse = TRUE)
        colnames(X_cc) <- character(0); return(X_cc)
    }
    feat_bins <- sort(unique(cc_sub$cc_bin))
    j_lookup  <- setNames(seq_along(feat_bins), feat_bins)
    i_vec <- match(cc_sub$.row, keep_idx); j_vec <- unname(j_lookup[cc_sub$cc_bin])
    ok <- !is.na(i_vec) & !is.na(j_vec); i_vec <- i_vec[ok]; j_vec <- j_vec[ok]
    if (length(i_vec)) {
        key <- paste0(i_vec, "_", j_vec); keepu <- !duplicated(key)
        i_vec <- i_vec[keepu]; j_vec <- j_vec[keepu]
    }
    Matrix::sparseMatrix(i=i_vec, j=j_vec, x=1,
                         dims=c(length(keep_idx), length(feat_bins)),
                         dimnames=list(NULL, paste0("cc_", feat_bins)))
}

# Robust quartile helpers (tolerant to ties)
.strict_breaks <- function(q) {
    q <- as.numeric(q); if (length(q) < 5) q <- rep_len(q, 5)
    rng <- diff(range(q[is.finite(q)], na.rm = TRUE)); step <- if (is.finite(rng) && rng > 0) rng * 1e-9 else 1e-9
    q[!is.finite(q)] <- 0
    for (i in 2:5) if (!(q[i] > q[i-1])) q[i] <- q[i-1] + step
    q[1] <- q[1] - step; q[5] <- q[5] + step
    q
}
.label_q4 <- function(x, brk) {
    b <- .strict_breaks(brk)
    idx <- findInterval(x, vec = b, rightmost.closed = TRUE, all.inside = TRUE)
    idx <- pmax(1L, pmin(4L, idx))
    factor(c("Q1","Q2","Q3","Q4")[idx], levels = c("Q1","Q2","Q3","Q4"))
}

# Build TEST split per mode
.make_te_split <- function(te_idx_all, mode, age_q_all, sex_vec) {
    if (mode == "age_sex") {
        split_key <- paste(addNA(age_q_all[te_idx_all]),
                           addNA(factor(sex_vec[te_idx_all], levels=c("M","F"))),
                           sep="|")
    } else if (mode == "age_only") {
        split_key <- as.character(addNA(age_q_all[te_idx_all]))
    } else if (mode == "sex_only") {
        split_key <- as.character(addNA(factor(sex_vec[te_idx_all], levels=c("M","F"))))
    } else { # none
        split_key <- rep("ALL", length(te_idx_all))
    }
    split(seq_along(te_idx_all), split_key, drop = TRUE)
}

# Backoff chain for TRAINING per mode
.backoff_levels_for_mode <- function(mode, age_key, sex_key) {
    switch(mode,
           "age_sex" = list(
               list(age=age_key, sex=sex_key, tag="AGE_Q×SEX"),
               list(age=age_key, sex=NA,     tag="AGE_Q"),
               list(age=NA,      sex=sex_key, tag="SEX"),
               list(age=NA,      sex=NA,     tag="ALL")
           ),
           "age_only" = list(
               list(age=age_key, sex=NA, tag="AGE_Q"),
               list(age=NA,      sex=NA, tag="ALL")
           ),
           "sex_only" = list(
               list(age=NA,      sex=sex_key, tag="SEX"),
               list(age=NA,      sex=NA,      tag="ALL")
           ),
           "none" = list(
               list(age=NA,      sex=NA,      tag="ALL")
           )
    )
}

# Prior getter honoring mode (pure: no reference to external `dat`)
.make_prior_getter <- function(classes_all, mode, cls_train, age_tr, sex_tr, alpha0) {
    df <- tibble(
        cls = as.character(cls_train),
        age = as.character(age_tr),
        sex = as.character(sex_tr)
    ) %>% filter(!is.na(cls), cls %in% classes_all)
    
    tabs <- list(
        by_age_sex = df %>% filter(!is.na(age), !is.na(sex)) %>% count(cls, age, sex, name="k"),
        by_age     = df %>% filter(!is.na(age))               %>% count(cls, age,       name="k"),
        by_sex     = df %>% filter(!is.na(sex))               %>% count(cls,     sex,   name="k"),
        by_all     = df %>%                                    count(cls,              name="k")
    )
    
    smoothed_norm <- function(named_counts) {
        v <- setNames(numeric(length(classes_all)), classes_all)
        if (!is.null(named_counts) && length(named_counts)) v[names(named_counts)] <- named_counts
        v <- v + alpha0
        v / sum(v)
    }
    
    function(age_key, sex_key) {
        if (mode == "age_sex") {
            a <- as.character(age_key); s <- as.character(sex_key)
            v <- tabs$by_age_sex %>% filter(age == a, sex == s)
            if (nrow(v)) return(smoothed_norm(setNames(v$k, v$cls)))
            v <- tabs$by_age %>% filter(age == a); if (nrow(v)) return(smoothed_norm(setNames(v$k, v$cls)))
            v <- tabs$by_sex %>% filter(sex == s); if (nrow(v)) return(smoothed_norm(setNames(v$k, v$cls)))
            if (nrow(tabs$by_all)) return(smoothed_norm(setNames(tabs$by_all$k, tabs$by_all$cls)))
        } else if (mode == "age_only") {
            a <- as.character(age_key)
            v <- tabs$by_age %>% filter(age == a)
            if (nrow(v)) return(smoothed_norm(setNames(v$k, v$cls)))
            if (nrow(tabs$by_all)) return(smoothed_norm(setNames(tabs$by_all$k, v$cls)))
        } else if (mode == "sex_only") {
            s <- as.character(sex_key)
            v <- tabs$by_sex %>% filter(sex == s)
            if (nrow(v)) return(smoothed_norm(setNames(v$k, v$cls)))
            if (nrow(tabs$by_all)) return(smoothed_norm(setNames(tabs$by_all$k, v$cls)))
        } else { # none
            if (nrow(tabs$by_all)) return(smoothed_norm(setNames(tabs$by_all$k, tabs$by_all$cls)))
        }
        # fall back uniform
        rep(1/length(classes_all), length(classes_all))
    }
}

# ============================== MAIN FUNCTION ==============================
compute_entropy_county_foreman <- function(
        ds, county_var,
        dict_dir      = NULL,
        icd_map_path  = if (!is.null(dict_dir)) file.path(dict_dir, "foreman-icd10-mapping.csv") else NULL,
        code_map_path = if (!is.null(dict_dir)) file.path(dict_dir, "foreman-table2-map.csv") else NULL
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
    
    # +++ NEW: register parallel backend if available
    if (.has_doParallel()) {
        cores <- parallel::detectCores()
        doParallel::registerDoParallel(cores = cores)
        if (verbose) message(sprintf("[diag] Parallel CV enabled (%d cores).", cores))
    } else {
        if (verbose) message("[diag] doParallel not installed; running CV in serial.")
    }
    
    dropped_log <- list()
    .log_drop <- function(g, key, stage, classes_before, support_vec, kept_classes, note = NA_character_) {
        before  <- sort(classes_before)
        kept    <- sort(kept_classes)
        dropped <- setdiff(before, kept)
        tibble(
            garbage_bin      = g,
            stratum          = key,
            stage            = stage,
            total_candidates = length(before),
            kept_n           = length(kept),
            dropped_n        = length(dropped),
            dropped_classes  = paste(dropped, collapse = "|"),
            kept_classes     = paste(kept, collapse = "|"),
            support_json     = as.character(jsonlite::toJSON(as.list(support_vec[before]), auto_unbox = TRUE)),
            note             = note
        )
    }
    if (verbose) {
        message(sprintf("[settings] cap/class=%d, CV folds=%d, nlambda=%d, pool years=%d, alpha=%.2f, lambda=%s, tau=%.2f",
                        capc, .cv_folds_default(), .nlambda_default(), pool_y,
                        .alpha_default(), .lambda_choice_default(), .tau_default()))
    }
    
    icd_map  <- .read_icd_map(icd_map_path)
    t2_map   <- .read_table2_map(code_map_path)
    icd_lu   <- .make_icd_lookup(icd_map)
    
    # Allowed/garbage sets
    garbage_bins_all <- sort(unique(t2_map$garbage))
    # Remove g3, g8, g9 entirely from redistribution
    garbage_bins     <- setdiff(garbage_bins_all, c("g3","g8","g9"))
    
    all_bins_in_map  <- sort(unique(icd_map$bin))
    valid_bins       <- setdiff(all_bins_in_map, sort(unique(t2_map$garbage)))
    if (length(valid_bins) < 2) stop("Too few valid bins found in dictionaries.")
    
    targets_by_g <- split(t2_map$target, t2_map$garbage, drop = TRUE) |>
        lapply(function(x) sort(unique(x)))
    
    # Basic columns
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
    
    # Normalize basics
    dat <- ds %>%
        dplyr::select(dplyr::all_of(keep)) %>%
        dplyr::mutate(
            uc4  = icd4(ucod),
            uc3  = icd3(ucod),
            ubin = map_icd_to_bin(ucod, icd_lu),
            sex2 = .normalize_sex(sex),
            age_num = dplyr::case_when(
                "age_years" %in% names(.) & !is.na(suppressWarnings(as.numeric(age_years))) ~ suppressWarnings(as.numeric(age_years)),
                "ager27" %in% names(.)    & !is.na(suppressWarnings(as.integer(ager27)))   ~ suppressWarnings(as.integer(ager27)),
                TRUE ~ NA_real_
            )
        )
    
    # warn up-front if sex is problematic
    .warn_if_sex_bad(dat$sex2, context = "full dataset")
    
    # Flags for garbage membership (UCOD-level; used for the new % metric)
    garbage_set_for_pct <- c("g1","g2","g4","g5","g6","g7")
    is_garb_pct <- !is.na(dat$ubin) & dat$ubin %in% garbage_set_for_pct
    
    # CC presence-only (exclude UCOD duplicates)
    cc_long <- if (length(rec_cols) == 0) {
        tibble(.row = integer(), cc_bin = character())
    } else {
        dat %>%
            mutate(.row = dplyr::row_number()) %>%
            dplyr::select(.row, dplyr::all_of(rec_cols), uc4, uc3) %>%
            tidyr::pivot_longer(cols = dplyr::all_of(rec_cols), names_to = "rec", values_to = "cc_raw") %>%
            mutate(cc_icd4 = icd4(cc_raw)) %>%
            filter(!(substr(cc_icd4, 1, 4) == uc4 | substr(cc_icd4, 1, 3) == uc3)) %>%
            transmute(.row, cc_bin = map_icd_to_bin(cc_icd4, icd_lu)) %>%
            filter(!is.na(cc_bin) & cc_bin != "") %>%
            distinct(.row, cc_bin)
    }
    
    row_is_garbage <- !is.na(dat$ubin) & dat$ubin %in% garbage_bins
    
    # Base valid counts (UCOD valid & not one of the active garbage bins)
    base_counts_agg <- dat %>%
        dplyr::filter(!is.na(ubin) & ubin %in% valid_bins & !row_is_garbage) %>%
        dplyr::count(.data[[county_var]], year, bin = ubin, name = "k", .drop = FALSE) %>%
        dplyr::rename(cnty = dplyr::all_of(county_var))
    
    # Per-bin age quartiles (fixed across years)
    .make_q4_breaks_for_bin <- function(g) {
        x <- dat$age_num[dat$ubin == g & is.finite(dat$age_num)]
        if (length(x) < 4 || length(unique(x)) < 4) {
            rng <- range(x, finite = TRUE); if (!all(is.finite(rng))) rng <- c(0, 100)
            q <- seq(rng[1], rng[2], length.out = 5L)
        } else {
            q <- as.numeric(stats::quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE, names = FALSE, type = 7))
        }
        .strict_breaks(q)
    }
    age_q_breaks <- setNames(vector("list", length(garbage_bins)), garbage_bins)
    for (g in garbage_bins) age_q_breaks[[g]] <- .make_q4_breaks_for_bin(g)
    
    prob_counts_list <- list()
    metrics_list     <- list()
    
    # ============================ per garbage bin ============================
    for (g in garbage_bins) {
        targets_raw <- targets_by_g[[g]]; if (is.null(targets_raw) || !length(targets_raw)) next
        classes_all_global <- sort(intersect(targets_raw, valid_bins)); if (!length(classes_all_global)) next
        
        rows_with_g_cc <- if (nrow(cc_long)) unique(cc_long$.row[cc_long$cc_bin == g]) else integer(0)
        te_idx_all <- intersect(which(row_is_garbage), which(dat$ubin == g))
        if (!length(te_idx_all)) next
        
        brk_g <- age_q_breaks[[g]]
        age_q_all <- .label_q4(dat$age_num, brk_g)
        
        # Per-bin strat mode (+ slice-wise effective mode)
        mode_g <- .strat_mode_for_bin(g)
        mode_eff <- .effective_mode_for_slice(mode_g, dat$sex2, te_idx_all)
        if (.verbose_default() && mode_eff != mode_g) {
            warning(sprintf("[DQ] Bin %s: switched stratification %s → %s (sex missing/constant in this slice).",
                            g, mode_g, mode_eff), call. = FALSE)
        }
        
        # TEST split by effective mode
        te_meta_all <- dat[te_idx_all, c(county_var, "year", "sex2")]
        colnames(te_meta_all)[1] <- "cnty"
        te_split <- .make_te_split(te_idx_all, mode_eff, age_q_all, dat$sex2)
        
        for (key in names(te_split)) {
            idx_local <- te_split[[key]]
            te_idx    <- te_idx_all[idx_local]
            n_te      <- length(te_idx); if (!n_te) next
            
            # representative keys for this split
            age_k <- if (mode_eff %in% c("age_sex","age_only")) age_q_all[te_idx][1] else NA
            sex_k <- if (mode_eff %in% c("age_sex","sex_only")) dat$sex2[te_idx][1]   else NA
            
            # TRAINING BACKOFF for this mode
            backoff_levels <- .backoff_levels_for_mode(mode_eff, age_k, sex_k)
            
            tr_idx <- integer(0)
            used_tag <- NA_character_
            ok_classes <- character(0)
            support_vec <- setNames(numeric(length(classes_all_global)), classes_all_global)
            
            for (lvl in backoff_levels) {
                # Vectorized masks
                age_mask <- if (is.na(lvl$age)) rep(TRUE, nrow(dat)) else (!is.na(age_q_all) & (age_q_all == lvl$age))
                sex_mask <- if (is.na(lvl$sex)) rep(TRUE, nrow(dat)) else (!is.na(dat$sex2)  & (dat$sex2  == lvl$sex))
                
                tr_all <- which(!is.na(dat$ubin) & (dat$ubin %in% classes_all_global) & age_mask & sex_mask)
                tr_idx_try <- intersect(rows_with_g_cc, tr_all)
                
                # temporal pooling
                pool_y_use <- .pool_years_default()
                if (pool_y_use > 0 && length(tr_idx_try)) {
                    te_years_probe <- unique(dat$year[te_idx])
                    tr_mask <- dat$year >= min(te_years_probe) - pool_y_use & dat$year <= max(te_years_probe) + pool_y_use
                    tr_idx_try <- intersect(tr_idx_try, which(tr_mask))
                }
                if (!length(tr_idx_try)) next
                
                # support per class in this widened stratum
                cls_counts_all <- table(dat$ubin[tr_idx_try])
                support_vec <- setNames(rep(0, length(classes_all_global)), classes_all_global)
                support_vec[names(cls_counts_all)] <- as.numeric(cls_counts_all)
                
                # Stage 1: MIN_PER_CLASS
                ok_minpc <- names(support_vec)[support_vec >= minpc]
                dropped_log[[length(dropped_log)+1L]] <- .log_drop(g, key, "after_minpc", classes_all_global, support_vec, ok_minpc, note = lvl$tag)
                if (length(ok_minpc) < 2) next
                
                # Stage 2: MIN_TRAIN_SUPPORT
                min_support <- .min_train_support_def()
                ok_support  <- names(support_vec)[support_vec >= max(minpc, min_support)]
                ok_classes_try  <- intersect(ok_minpc, ok_support)
                dropped_log[[length(dropped_log)+1L]] <- .log_drop(g, key, "after_support", classes_all_global, support_vec, ok_classes_try, note = lvl$tag)
                if (length(ok_classes_try) < 2) next
                
                # success at this backoff level (no k-cap)
                tr_idx    <- tr_idx_try
                ok_classes <- sort(ok_classes_try)
                used_tag  <- lvl$tag
                break
            }
            
            # If training still failed, we now simply skip (no fallback)
            if (!length(tr_idx) || length(ok_classes) < 2) {
                dropped_log[[length(dropped_log)+1L]] <- .log_drop(
                    g, key, "skipped_stratum", classes_all_global, support_vec, ok_classes,
                    note = "no_backoff_succeeded"
                )
                next
            }
            
            if (verbose) {
                tbl <- sort(table(dat$ubin[tr_idx]), decreasing = TRUE)
                msg <- paste(head(paste(names(tbl), as.integer(tbl), sep=":"), 6), collapse=", ")
                message(sprintf("[diag] g=%s stratum=%s TRAIN=%s n_tr=%d n_te=%d k=%d; top classes: %s",
                                g, key, used_tag, length(tr_idx), n_te, length(ok_classes), msg))
            }
            
            # Class-balanced training sample
            target_per_class <- .target_per_class_def()
            tr_by_class <- split(tr_idx, dat$ubin[tr_idx])
            take <- unlist(lapply(tr_by_class[names(tr_by_class) %in% ok_classes], function(ix) {
                n <- length(ix)
                if (n >= target_per_class) sample(ix, target_per_class)
                else c(ix, sample(ix, target_per_class - n, replace = TRUE))
            }), use.names = FALSE)
            tr_idx2 <- sort(unique(take))
            if (!length(tr_idx2)) {
                dropped_log[[length(dropped_log)+1L]] <- .log_drop(g, key, "skipped_stratum",
                                                                   classes_all_global, support_vec, ok_classes,
                                                                   note = paste0("empty_after_balancing:", used_tag))
                next
            }
            
            # Build features
            keep_idx <- sort(unique(c(tr_idx2, te_idx)))
            X_local  <- .build_X_for_rows(keep_idx, cc_long)
            nz <- if (ncol(X_local)) which(Matrix::colSums(X_local[match(tr_idx2, keep_idx), , drop = FALSE]) > 0) else integer(0)
            if (!length(nz)) {
                dropped_log[[length(dropped_log)+1L]] <- .log_drop(g, key, "skipped_stratum",
                                                                   classes_all_global, support_vec, ok_classes,
                                                                   note = paste0("no_nz_features:", used_tag))
                next
            }
            X_tr <- X_local[match(tr_idx2, keep_idx), nz, drop = FALSE]
            X_te <- X_local[match(te_idx,  keep_idx), nz, drop = FALSE]
            rm(X_local); gc(FALSE)
            
            y_tr <- factor(dat$ubin[tr_idx2], levels = ok_classes)
            
            fit <- suppressWarnings(tryCatch(
                glmnet::cv.glmnet(
                    X_tr, y_tr,
                    family = "multinomial",
                    alpha  = .alpha_default(),
                    type.multinomial = "grouped",
                    nfolds = .cv_folds_default(),
                    nlambda = .nlambda_default(),
                    lambda.min.ratio = .lmr_default(),
                    standardize = FALSE,
                    parallel = .has_doParallel()   # <── parallel CV when doParallel is available
                ),
                error = function(e) NULL
            ))
            if (is.null(fit)) {
                dropped_log[[length(dropped_log)+1L]] <- .log_drop(g, key, "skipped_stratum",
                                                                   classes_all_global, support_vec, ok_classes,
                                                                   note = paste0("glmnet_fail:", used_tag))
                next
            }
            
            lambda_choice <- .lambda_choice_default()
            M_raw <- .extract_multinomial_probs(fit, X_te, ok_classes, s = lambda_choice)
            if (is.null(M_raw)) {
                dropped_log[[length(dropped_log)+1L]] <- .log_drop(g, key, "skipped_stratum",
                                                                   classes_all_global, support_vec, ok_classes,
                                                                   note = paste0("predict_fail:", used_tag))
                next
            }
            shp <- .fix_prob_shape(M_raw, n_te = n_te, classes = ok_classes,
                                   context = paste0("g=", g, ", stratum=", key, ", dims_in=", paste(dim(M_raw), collapse="×")))
            probs <- shp$M
            rs0 <- rowSums(probs); rs0[rs0 == 0] <- 1; probs <- probs / rs0
            
            # Prior (for KL/JSD; respects effective mode)
            prior_getter <- .make_prior_getter(
                classes_all = ok_classes,
                mode        = mode_eff,
                cls_train   = as.character(dat$ubin[tr_idx2]),
                age_tr      = as.character(age_q_all[tr_idx2]),
                sex_tr      = as.character(dat$sex2[tr_idx2]),
                alpha0      = alpha0
            )
            age_te_q <- age_q_all[te_idx]; sex_te <- dat$sex2[te_idx]
            prior_mat <- matrix(NA_real_, nrow = length(te_idx), ncol = length(ok_classes))
            for (i in seq_along(te_idx)) prior_mat[i, ] <- prior_getter(age_te_q[i], sex_te[i])
            colnames(prior_mat) <- ok_classes
            
            # Aggregate to county×year×bin
            te_meta <- dat[te_idx, c(county_var, "year")]; colnames(te_meta) <- c("cnty", "year")
            cls_summaries <- lapply(seq_along(ok_classes), function(j) {
                tibble(cnty = te_meta$cnty, year = te_meta$year, k = probs[, j]) %>%
                    group_by(cnty, year) %>% summarise(k = sum(k), .groups = "drop") %>%
                    mutate(bin = ok_classes[j])
            })
            prob_counts_list[[length(prob_counts_list)+1L]] <- bind_rows(cls_summaries)
            
            # Metrics
            eps <- .Machine$double.eps
            probs_for_metrics <- probs
            tau <- .tau_default()
            if (!is.na(tau) && tau > 0 && tau != 1) {
                probs_for_metrics <- probs_for_metrics ^ tau
                probs_for_metrics <- probs_for_metrics / rowSums(probs_for_metrics)
            }
            H_post  <- rowSums(-probs_for_metrics * log(pmax(probs_for_metrics,  eps)))
            k_cand  <- length(ok_classes)
            Hnorm_gc_by_kcand <- if (k_cand > 1) H_post / log(k_cand) else rep(0, length(H_post))
            
            P0 <- pmax(probs,     eps); P0 <- P0 / rowSums(P0)
            Q0 <- pmax(prior_mat, eps); Q0 <- Q0 / rowSums(Q0)
            KL_raw_vec  <- rowSums(P0 * (log(P0) - log(Q0)))
            KL_norm_vec <- if (k_cand > 1) KL_raw_vec / log(k_cand) else rep(0, length(KL_raw_vec))
            
            M0 <- 0.5 * (P0 + Q0)
            KL_PM <- rowSums(P0 * (log(P0) - log(pmax(M0, eps))))
            KL_QM <- rowSums(Q0 * (log(Q0) - log(pmax(M0, eps))))
            JSD_vec <- 0.5 * (KL_PM + KL_QM)
            JSD_norm_vec <- JSD_vec / log(2)
            
            metr <- tibble(
                cnty = te_meta$cnty, year = te_meta$year,
                sum_KL_norm           = KL_norm_vec,
                sum_Hnorm_gc_by_kcand = Hnorm_gc_by_kcand,
                sum_JSD_norm          = JSD_norm_vec,
                one = 1L
            ) %>%
                group_by(cnty, year) %>%
                summarise(
                    sum_KL_norm           = sum(sum_KL_norm,           na.rm = TRUE),
                    sum_Hnorm_gc_by_kcand = sum(sum_Hnorm_gc_by_kcand, na.rm = TRUE),
                    sum_JSD_norm          = sum(sum_JSD_norm,          na.rm = TRUE),
                    N_garbage             = sum(one),
                    .groups = "drop"
                )
            metrics_list[[length(metrics_list)+1L]] <- metr
        } # end stratum loop
    } # end g loop
    
    prob_counts_all <- if (length(prob_counts_list)) bind_rows(prob_counts_list) else
        tibble(cnty = character(), year = numeric(), k = numeric(), bin = character())
    
    # Observed valid UCOD + redistributed garbage mass
    cty_year_bin <- bind_rows(
        base_counts_agg %>% select(cnty, year, bin, k),
        prob_counts_all  %>% select(cnty, year, bin, k)
    ) %>%
        group_by(cnty, year, bin) %>%
        summarise(k = sum(k), .groups = "drop")
    
    out_aggregate <- cty_year_bin %>%
        group_by(cnty, year) %>%
        summarise(
            total      = sum(k),
            DQ_overall = { p <- k / sum(k); .shannon_H(p) / log(length(unique(bin))) },
            .groups = "drop"
        ) %>%
        mutate(DQ_overall = pmax(0, pmin(1, DQ_overall)))
    
    # --- NEW: % of certificates whose UCOD is in {g1,g2,g4,g5,g6,g7} ---
    counts_total <- dat %>% count(.data[[county_var]], year, name = "N_total", .drop = FALSE) %>%
        rename(cnty = dplyr::all_of(county_var))
    counts_gsel  <- dat %>% mutate(is_sel = is_garb_pct) %>%
        group_by(.data[[county_var]], year) %>%
        summarise(N_garb_g1g2g4g5g6g7 = sum(is_sel, na.rm = TRUE), .groups = "drop") %>%
        rename(cnty = dplyr::all_of(county_var))
    
    # Per-record/garbage metrics
    metrics_sum <- if (length(metrics_list)) {
        bind_rows(metrics_list) %>%
            group_by(cnty, year) %>%
            summarise(
                sum_KL_norm           = sum(sum_KL_norm,           na.rm = TRUE),
                sum_Hnorm_gc_by_kcand = sum(sum_Hnorm_gc_by_kcand, na.rm = TRUE),
                sum_JSD_norm          = sum(sum_JSD_norm,          na.rm = TRUE),
                N_garbage             = sum(N_garbage),
                .groups = "drop"
            )
    } else {
        tibble(
            cnty = character(), year = numeric(),
            sum_KL_norm = numeric(), sum_Hnorm_gc_by_kcand = numeric(),
            sum_JSD_norm = numeric(), N_garbage = integer()
        )
    }
    
    out_perrecord <- counts_total %>%
        left_join(metrics_sum,  by = c("cnty", "year")) %>%
        left_join(counts_gsel,  by = c("cnty", "year")) %>%
        mutate(
            across(c(sum_KL_norm, sum_Hnorm_gc_by_kcand, sum_JSD_norm), ~ dplyr::coalesce(.x, 0)),
            N_garbage             = dplyr::coalesce(N_garbage, 0L),
            N_garb_g1g2g4g5g6g7   = dplyr::coalesce(N_garb_g1g2g4g5g6g7, 0L),
            pct_garb_g1g2g4g5g6g7 = ifelse(N_total > 0, N_garb_g1g2g4g5g6g7 / N_total, NA_real_),
            RI          = dplyr::if_else(N_garbage > 0, sum_KL_norm / N_garbage, NA_real_),
            RI_post_only= dplyr::if_else(N_garbage > 0, 1 - (sum_Hnorm_gc_by_kcand / N_garbage), NA_real_),
            RI_jsd      = dplyr::if_else(N_garbage > 0, sum_JSD_norm / N_garbage, NA_real_)
        ) %>%
        mutate(
            RI          = ifelse(is.na(RI),          RI,          pmax(0, pmin(1, RI))),
            RI_post_only= ifelse(is.na(RI_post_only),RI_post_only,pmax(0, pmin(1, RI_post_only))),
            RI_jsd      = ifelse(is.na(RI_jsd),      RI_jsd,      pmax(0, pmin(1, RI_jsd)))
        ) %>%
        select(cnty, year, RI, RI_post_only, RI_jsd, N_garbage, N_total,
               N_garb_g1g2g4g5g6g7, pct_garb_g1g2g4g5g6g7)
    
    out <- out_aggregate %>%
        left_join(out_perrecord, by = c("cnty", "year")) %>%
        rename(!!county_var := cnty)
    
    # Diagnostics
    if (length(dropped_log)) {
        drop_df <- bind_rows(dropped_log)
        ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
        f2 <- if (!is.null(diag_dir)) {
            if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
            file.path(diag_dir, paste0("dropped_classes_", ts, ".csv"))
        } else {
            file.path(getwd(), paste0("dropped_classes_", ts, ".csv"))
        }
        readr::write_csv(drop_df, f2)
        if (verbose) message("[diag] wrote dropped-classes report to: ", f2)
        attr(out, "dropped_classes") <- drop_df
    }
    
    out
}
