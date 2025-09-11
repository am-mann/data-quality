# ───────────────── dq_entropy_helper.R — Foreman bins + RI (county×year) ─────────────────
# Inputs:
#   - ds: data frame with 'ucod', county key, optional 'age_years' or 'ager27',
#         and contributing causes in record_1..record_20 (or record1..record20)
#   - county_var: name of the county column in ds (e.g., "county_ihme")
#   - dict_dir: directory containing:
#       * foreman-icd10-mapping.csv  (ICD → USCOD *codes*, NOT names, e.g. A_5, B_3_1, G_1)
#       * foreman-table2-map.csv     (wide: target_cause + G_1..G_9 TRUE/FALSE)
#
# Output columns (joined by county×year):
#   total, DQ_entropy, DQ_K, DQ_overall, DQ_expH_over_K,
#   RI, DQ_rec_expH_over_K_mean, DQ_rec_ig_abs_mean, DQ_rec_ig_frac_mean_garbage,
#   foreman_garbage, foreman_garbage_adj
# ──────────────────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr)
    library(Matrix)
})

.has_glmnet <- function() requireNamespace("glmnet", quietly = TRUE)

# ---------- tiny helpers ----------
canonical_icd <- function(x) stringr::str_remove_all(stringr::str_to_upper(x), "[^A-Z0-9]")
clean_icd     <- function(x) stringr::str_remove_all(stringr::str_to_upper(x), "[^A-Z0-9\\.]")
icd3          <- function(x) substr(canonical_icd(x), 1, 3)
icd4          <- function(x) substr(canonical_icd(x), 1, 4)
# keep letters+digits; drop punctuation/underscores. "B_3_1"→"b31", "G_1"→"g1"
norm_code     <- function(x) tolower(gsub("[^A-Za-z0-9]+", "", as.character(x)))

# Accept record1..record20 or record_1..record_20
.rename_record_cols <- function(df) {
    if (!any(grepl("^record_", names(df))) && any(grepl("^record[0-9]+$", names(df)))) {
        dplyr::rename_with(df, ~ sub("^record([0-9]+)$","record_\\1", .x), .cols = dplyr::matches("^record[0-9]+$"))
    } else df
}

# ---------- robust readers (force character types so USCOD never becomes numeric) ----------
.read_icd_map <- function(path) {
    df <- readr::read_csv(
        path, show_col_types = FALSE,
        col_types = readr::cols(.default = readr::col_character(),
                                ICD10 = readr::col_character(),
                                USCOD = readr::col_character())
    )
    if (!all(c("ICD10","USCOD") %in% names(df))) {
        stop("foreman-icd10-mapping.csv must have columns ICD10 and USCOD (codes like A_5, B_3_1, G_1).")
    }
    out <- df %>%
        transmute(
            icd = canonical_icd(ICD10),
            bin = norm_code(USCOD)   # preserve letters, squash punctuation
        ) %>%
        filter(icd != "", bin != "") %>%
        distinct(icd, .keep_all = TRUE)
    
    # Guard: bins must not be digits-only (would indicate numeric parsing)
    if (any(grepl("^\\d+$", out$bin))) {
        bad <- unique(out$bin[grepl("^\\d+$", out$bin)])
        stop("USCOD normalized to digits-only bins (e.g. ", paste(head(bad, 8), collapse=", "),
             "). Ensure USCOD column is alphanumeric codes like A_5, B_3_1, G_1.")
    }
    out
}

.read_table2_map <- function(path) {
    df <- readr::read_csv(path, show_col_types = FALSE,
                          col_types = readr::cols(.default = readr::col_character()))
    if (!"target_cause" %in% names(df)) {
        stop("foreman-table2-map.csv must contain column 'target_cause'.")
    }
    g_cols <- grep("^G_[0-9]+$", names(df), value = TRUE)
    if (!length(g_cols)) stop("foreman-table2-map.csv must contain columns G_1..G_9.")
    
    # Coerce truthy flags to logical
    df[g_cols] <- lapply(df[g_cols], function(x){
        y <- tolower(as.character(x))
        y %in% c("1","true","t","y","yes","x","✓","check","checked")
    })
    
    t2_long <- df %>%
        mutate(target = norm_code(target_cause)) %>%
        tidyr::pivot_longer(all_of(g_cols), names_to = "garbage_raw", values_to = "flag") %>%
        filter(flag) %>%
        transmute(
            garbage = norm_code(garbage_raw),  # "G_1" -> "g1"
            target  = norm_code(target)
        ) %>%
        filter(garbage != "", target != "") %>%
        distinct(garbage, target)
    
    need <- paste0("g", 1:9)
    miss <- setdiff(need, unique(t2_long$garbage))
    if (length(miss)) stop("Missing garbage bins in Table 2 after normalization: ", paste(miss, collapse=", "))
    t2_long
}

# Lookup tables for 4-char then 3-char ICD matching
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

# Entropy (natural log; base cancels when normalized by log(k))
.shannon_H <- function(p) { p <- p[p > 0]; if (!length(p)) 0 else -sum(p * log(p)) }

# Robust glmnet prediction: fixes dim/labels/row-mismatch
.safe_multinomial_probs <- function(fit, X_te, classes, s = "lambda.1se") {
    arr <- tryCatch(predict(fit, newx = X_te, type = "response", s = s), error = function(e) NULL)
    if (is.null(arr)) return(NULL)
    probs <- if (length(dim(arr)) == 3L) arr[, , 1, drop = FALSE] else arr
    probs <- as.matrix(probs)
    
    if (is.null(colnames(probs)) || anyNA(colnames(probs)) || ncol(probs) == 0L) {
        dn <- dimnames(arr)
        if (!is.null(dn) && length(dn) >= 2L && !is.null(dn[[2]]) && length(dn[[2]]) == ncol(probs)) {
            colnames(probs) <- dn[[2]]
        } else if (!is.null(fit$glmnet.fit$classnames) && length(fit$glmnet.fit$classnames) == ncol(probs)) {
            colnames(probs) <- fit$glmnet.fit$classnames
        } else {
            colnames(probs) <- paste0("cls", seq_len(max(1L, ncol(probs))))
        }
    }
    
    n_te <- nrow(X_te)
    if (nrow(probs) != n_te) {                # align rows if glmnet returns fewer/more
        k <- max(1L, ncol(probs))
        vals <- as.numeric(probs)
        if (length(vals) != n_te * k) vals <- rep_len(vals, n_te * k)
        probs <- matrix(vals, nrow = n_te, ncol = k, byrow = TRUE, dimnames = list(NULL, colnames(probs)))
    }
    
    present <- intersect(classes, colnames(probs))
    base <- if (length(present)) probs[, present, drop = FALSE] else matrix(0, nrow = n_te, ncol = 0)
    missing <- setdiff(classes, colnames(probs))
    if (length(missing)) base <- cbind(base, matrix(0, nrow = n_te, ncol = length(missing), dimnames = list(NULL, missing)))
    base <- base[, classes, drop = FALSE]
    
    rs <- rowSums(base); rs[rs == 0] <- 1
    base / rs
}

.pool_years_default      <- function() as.integer(getOption("DQ_POOL_YEARS", 2L))      # years to include on either side
.dirichlet_prior_default <- function() as.numeric(getOption("DQ_DIRICHLET_PRIOR", 0.5))# symmetric prior for fallback
.min_per_class_default   <- function() as.integer(getOption("DQ_MIN_PER_CLASS", 12L)) # glmnet training guard
.max_train_per_class_def <- function() as.integer(getOption("DQ_MAX_TRAIN_PER_CLASS", 8000L))

# New runtime knobs (override with options(DQ_POOL_YEARS=2, DQ_DIRICHLET_PRIOR=0.5, ...))
.pool_years_default      <- function() as.integer(getOption("DQ_POOL_YEARS", 2L))      # years to include on either side
.dirichlet_prior_default <- function() as.numeric(getOption("DQ_DIRICHLET_PRIOR", 0.5))# symmetric prior for fallback
.min_per_class_default   <- function() as.integer(getOption("DQ_MIN_PER_CLASS", 12L)) # glmnet training guard
.max_train_per_class_def <- function() as.integer(getOption("DQ_MAX_TRAIN_PER_CLASS", 8000L))

compute_entropy_county_foreman<- function(
        ds,
        county_var,
        dict_dir      = NULL,
        icd_map_path  = if (!is.null(dict_dir)) file.path(dict_dir, "foreman-icd10-mapping.csv") else NULL,
        code_map_path = if (!is.null(dict_dir)) file.path(dict_dir, "foreman-table2-map.csv") else NULL,
        age_breaks    = c(-Inf, 1, 5, seq(10, 85, by = 5), Inf)
) {
    suppressPackageStartupMessages({
        library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr); library(Matrix)
    })
    .has_glmnet <- function() requireNamespace("glmnet", quietly = TRUE)
    
    canonical_icd <- function(x) stringr::str_remove_all(stringr::str_to_upper(x), "[^A-Z0-9]")
    icd3          <- function(x) substr(canonical_icd(x), 1, 3)
    icd4          <- function(x) substr(canonical_icd(x), 1, 4)
    norm_code     <- function(x) tolower(gsub("[^A-Za-z0-9]+", "", as.character(x)))
    
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
        if (!all(c("ICD10","USCOD") %in% names(df))) {
            stop("foreman-icd10-mapping.csv must have columns ICD10 and USCOD.")
        }
        out <- df %>%
            transmute(icd = canonical_icd(ICD10), bin = norm_code(USCOD)) %>%
            filter(icd != "", bin != "") %>%
            distinct(icd, .keep_all = TRUE)
        if (any(grepl("^\\d+$", out$bin))) stop("USCOD parsed as digits-only; fix the mapping file.")
        out
    }
    .read_table2_map <- function(path) {
        df <- readr::read_csv(path, show_col_types = FALSE,
                              col_types = readr::cols(.default = readr::col_character()))
        if (!"target_cause" %in% names(df)) stop("foreman-table2-map.csv must contain 'target_cause'.")
        g_cols <- grep("^G_[0-9]+$", names(df), value = TRUE)
        if (!length(g_cols)) stop("foreman-table2-map.csv must contain columns G_1..")
        df[g_cols] <- lapply(df[g_cols], function(x){
            y <- tolower(as.character(x)); y %in% c("1","true","t","y","yes","x","✓","check","checked")
        })
        df %>%
            mutate(target = norm_code(target_cause)) %>%
            pivot_longer(all_of(g_cols), names_to = "garbage_raw", values_to = "flag") %>%
            filter(flag) %>%
            transmute(garbage = norm_code(garbage_raw), target = norm_code(target)) %>%
            distinct(garbage, target)
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
    .shannon_H <- function(p) { p <- p[p > 0]; if (!length(p)) 0 else -sum(p * log(p)) }
    
    if (is.null(icd_map_path) || is.null(code_map_path))
        stop("Provide icd_map_path and code_map_path (or dict_dir).")
    
    icd_map  <- .read_icd_map(icd_map_path)
    t2_map   <- .read_table2_map(code_map_path)
    icd_lu   <- .make_icd_lookup(icd_map)
    
    garbage_bins     <- sort(unique(t2_map$garbage))
    all_bins_in_map  <- sort(unique(icd_map$bin))
    valid_bins       <- setdiff(all_bins_in_map, garbage_bins)
    K <- length(valid_bins)
    if (K < 10) warning("Few valid bins detected (K=", K, "). Check dictionaries.")
    
    targets_by_g <- split(t2_map$target, t2_map$garbage, drop = TRUE) %>%
        lapply(function(x) sort(unique(x)))
    
    # ---- standardize DS ----
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
    
    # ---------- CC long (for filtering only) ----------
    cc_long <- if (length(rec_cols) == 0) {
        tibble::tibble(.row = integer(), cc_bin = character())
    } else {
        dat %>%
            mutate(.row = row_number()) %>%
            select(.row, all_of(rec_cols), uc4, uc3) %>%
            tidyr::pivot_longer(cols = all_of(rec_cols), names_to = "rec", values_to = "cc_raw") %>%
            mutate(cc_icd4 = icd4(cc_raw)) %>%
            filter(!(substr(cc_icd4, 1, 4) == uc4 | substr(cc_icd4, 1, 3) == uc3)) %>%
            transmute(.row, cc_bin = map_icd_to_bin(cc_icd4, icd_lu)) %>%
            filter(!is.na(cc_bin) & cc_bin != "")
    }
    
    row_is_garbage <- !is.na(dat$ubin) & dat$ubin %in% garbage_bins
    
    # ---------- base counts (non-garbage UCODs) ----------
    base_counts_agg <- dat %>%
        filter(!row_is_garbage & !is.na(ubin) & ubin %in% valid_bins) %>%
        count(.data[[county_var]], year, bin = ubin, name = "k", .drop = FALSE)
    
    counts_cy <- dat %>% count(.data[[county_var]], year, name = "N_total", .drop = FALSE)
    
    # ---------- INFORMATIVE PRIORS by g (from all years in `dat`) ----------
    prior_by_g <- lapply(garbage_bins, function(g) {
        targs <- targets_by_g[[g]]
        if (is.null(targs) || !length(targs)) return(NULL)
        if (!nrow(cc_long)) {
            rep(1/length(targs), length(targs)) %>% setNames(targs)
        } else {
            rows_with_g_cc_all <- unique(cc_long$.row[cc_long$cc_bin == g])
            tr_all <- intersect(rows_with_g_cc_all, which(!is.na(dat$ubin) & dat$ubin %in% targs))
            if (!length(tr_all)) {
                rep(1/length(targs), length(targs)) %>% setNames(targs)
            } else {
                tbl <- table(dat$ubin[tr_all])
                vec <- rep(0, length(targs)); names(vec) <- targs
                vec[names(tbl)] <- as.numeric(tbl)
                if (sum(vec) == 0) rep(1/length(targs), length(targs)) else (vec / sum(vec))
            }
        }
    })
    names(prior_by_g) <- garbage_bins
    alpha0 <- getOption("DQ_DIRICHLET_PRIOR", 1.0)
    
    # ---------- holders ----------
    prob_counts_list <- list()
    metrics_list     <- list()
    list_i <- 0L
    list_m <- 0L
    
    # feature builder
    build_X_for_rows <- function(keep_idx) {
        # CC sparse
        if (nrow(cc_long) == 0 || length(keep_idx) == 0) {
            X_cc <- Matrix::Matrix(0, nrow = length(keep_idx), ncol = 0, sparse = TRUE)
            colnames(X_cc) <- character(0)
        } else {
            cc_sub <- cc_long %>%
                filter(.row %in% keep_idx, cc_bin %in% valid_bins) %>%
                distinct(.row, cc_bin)
            if (!nrow(cc_sub)) {
                X_cc <- Matrix::Matrix(0, nrow = length(keep_idx), ncol = 0, sparse = TRUE)
                colnames(X_cc) <- character(0)
            } else {
                feat_bins_local <- sort(unique(cc_sub$cc_bin))
                j_lookup <- setNames(seq_along(feat_bins_local), feat_bins_local)
                i_vec <- match(cc_sub$.row, keep_idx)
                j_vec <- unname(j_lookup[cc_sub$cc_bin])
                keep_ok <- !is.na(i_vec) & !is.na(j_vec)
                i_vec <- i_vec[keep_ok]; j_vec <- j_vec[keep_ok]
                if (length(i_vec)) {
                    key <- paste0(i_vec, "_", j_vec)
                    which_unique <- !duplicated(key)
                    i_vec <- i_vec[which_unique]; j_vec <- j_vec[which_unique]
                }
                X_cc <- Matrix::sparseMatrix(
                    i = i_vec, j = j_vec, x = 1,
                    dims = c(length(keep_idx), length(feat_bins_local)),
                    dimnames = list(NULL, feat_bins_local)
                )
            }
        }
        
        # Age one-hot
        age_vec <- if ("age_years" %in% names(dat)) suppressWarnings(as.numeric(dat$age_years)) else NA_real_
        if (all(is.na(age_vec)) && "ager27" %in% names(dat)) {
            age_vec <- ifelse(is.na(dat$ager27), NA_real_, ifelse(as.integer(dat$ager27) >= 12, 40, 10))
        }
        age_keep <- age_vec[keep_idx]
        age_bucket <- cut(age_keep, breaks = age_breaks, right = FALSE, include.lowest = TRUE)
        levs <- levels(age_bucket); P <- length(levs)
        if (P > 0) {
            ok <- !is.na(age_bucket)
            i  <- which(ok)
            j  <- match(age_bucket[ok], levs)
            age_mm <- Matrix::sparseMatrix(
                i = i, j = j, x = 1,
                dims = c(length(keep_idx), P),
                dimnames = list(NULL, paste0("age_", levs))
            )
        } else {
            age_mm <- Matrix::Matrix(0, nrow = length(keep_idx), ncol = 0, sparse = TRUE)
        }
        
        if (ncol(age_mm) > 0 && ncol(X_cc) > 0) cbind(X_cc, age_mm)
        else if (ncol(age_mm) > 0) age_mm
        else X_cc
    }
    
    # ---------- loop garbage groups ----------
    min_per_class <- as.integer(getOption("DQ_MIN_PER_CLASS", 12L))
    
    for (g in garbage_bins) {
        targets_g <- targets_by_g[[g]]
        if (is.null(targets_g) || !length(targets_g)) next
        
        rows_with_g_cc <- if (nrow(cc_long)) unique(cc_long$.row[cc_long$cc_bin == g]) else integer(0)
        tr_idx <- intersect(rows_with_g_cc, which(!is.na(dat$ubin) & dat$ubin %in% targets_g))
        te_idx <- intersect(which(row_is_garbage), which(dat$ubin == g))
        if (!length(te_idx)) next
        
        classes_all <- sort(unique(targets_g))
        classes_obs <- sort(intersect(classes_all, unique(dat$ubin[tr_idx])))
        
        model_ok <- FALSE; probs <- NULL
        
        # Try glmnet if we have ≥2 classes with enough counts
        if (.has_glmnet() && length(tr_idx) >= 30) {
            cls_counts <- table(dat$ubin[tr_idx])
            ok_classes <- names(cls_counts)[cls_counts >= min_per_class]
            if (length(ok_classes) >= 2) {
                cap <- as.integer(getOption("DQ_MAX_TRAIN_PER_CLASS", 8000L))
                take <- unlist(lapply(split(tr_idx, dat$ubin[tr_idx]),
                                      function(ix) if (length(ix) > cap) sample(ix, cap) else ix),
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
                                          type.multinomial = "ungrouped", nfolds = 3, parallel = FALSE),
                        error = function(e) NULL
                    ))
                    if (!is.null(fit)) {
                        .safe_multinomial_probs <- function(fit, X_te, classes, s = "lambda.1se") {
                            arr <- tryCatch(predict(fit, newx = X_te, type = "response", s = s), error = function(e) NULL)
                            if (is.null(arr)) return(NULL)
                            probs <- if (length(dim(arr)) == 3L) arr[, , 1, drop = FALSE] else arr
                            probs <- as.matrix(probs)
                            if (is.null(colnames(probs)) || anyNA(colnames(probs)) || ncol(probs) == 0L) {
                                dn <- dimnames(arr)
                                if (!is.null(dn) && length(dn) >= 2L && !is.null(dn[[2]]) && length(dn[[2]]) == ncol(probs)) {
                                    colnames(probs) <- dn[[2]]
                                } else if (!is.null(fit$glmnet.fit$classnames) && length(fit$glmnet.fit$classnames) == ncol(probs)) {
                                    colnames(probs) <- fit$glmnet.fit$classnames
                                } else {
                                    colnames(probs) <- paste0("cls", seq_len(max(1L, ncol(probs))))
                                }
                            }
                            n_te <- nrow(X_te)
                            if (nrow(probs) != n_te) {
                                k <- max(1L, ncol(probs))
                                vals <- as.numeric(probs)
                                if (length(vals) != n_te * k) vals <- rep_len(vals, n_te * k)
                                probs <- matrix(vals, nrow = n_te, ncol = k, byrow = TRUE, dimnames = list(NULL, colnames(probs)))
                            }
                            present <- intersect(classes, colnames(probs))
                            base <- if (length(present)) probs[, present, drop = FALSE] else matrix(0, nrow = n_te, ncol = 0)
                            missing <- setdiff(classes, colnames(probs))
                            if (length(missing)) base <- cbind(base, matrix(0, nrow = n_te, ncol = length(missing), dimnames = list(NULL, missing)))
                            base <- base[, classes, drop = FALSE]
                            rs <- rowSums(base); rs[rs == 0] <- 1
                            base / rs
                        }
                        probs <- .safe_multinomial_probs(fit, X_te, classes = classes_all, s = "lambda.1se")
                        model_ok <- !is.null(probs)
                    }
                }
            }
        }
        
        # Dirichlet-smoothed fallback (no NA even if k_obs < 2 or 0)
        if (!model_ok) {
            counts <- as.numeric(table(factor(dat$ubin[tr_idx], levels = classes_all)))
            prior_vec <- prior_by_g[[g]]
            prior_aligned <- prior_vec[match(classes_all, names(prior_vec))]
            if (!length(prior_aligned) || any(is.na(prior_aligned))) prior_aligned <- rep(1/length(classes_all), length(classes_all))
            prior_aligned <- prior_aligned / sum(prior_aligned)
            
            w <- (counts + alpha0 * prior_aligned) / (sum(counts) + alpha0)
            probs <- matrix(rep(w, each = length(te_idx)),
                            nrow = length(te_idx), ncol = length(classes_all),
                            dimnames = list(NULL, classes_all))
        }
        
        # ---- aggregate expected counts for this group's garbage rows ----
        te_meta <- dat[te_idx, c(county_var, "year")]
        colnames(te_meta) <- c("cnty", "year")
        cls_summaries <- lapply(seq_along(classes_all), function(j) {
            tibble::tibble(cnty = te_meta$cnty, year = te_meta$year, k = probs[, j]) %>%
                group_by(cnty, year) %>%
                summarise(k = sum(k), .groups = "drop") %>%
                mutate(bin = classes_all[j])
        })
        list_i <- list_i + 1L
        prob_counts_list[[list_i]] <- dplyr::bind_rows(cls_summaries)
        
        # ---- per-record metric partials (garbage rows only) ----
        eps <- .Machine$double.eps
        H_post <- rowSums(-probs * log(pmax(probs, eps)))
        k_cand <- length(classes_all)
        H0     <- log(k_cand)
        
        Hnorm_post_gc_by_kcand <- if (k_cand > 1) H_post / log(k_cand) else rep(0, length(H_post))
        Hnorm_post_K    <- if (K > 0) H_post / log(K) else rep(NA_real_, length(H_post))
        expH_over_K_vec <- if (K > 0) exp(H_post) / K else rep(NA_real_, length(H_post))
        IG_abs_vec      <- H0 - H_post
        IG_frac_vec     <- if (H0 > 0) (H0 - H_post) / H0 else rep(NA_real_, length(H_post))
        
        metr <- tibble::tibble(
            cnty = te_meta$cnty, year = te_meta$year,
            sum_Hnorm_post_K        = Hnorm_post_K,
            sum_expH_over_K         = expH_over_K_vec,
            sum_IG_abs              = IG_abs_vec,
            sum_IG_frac             = IG_frac_vec,
            sum_Hnorm_gc_by_kcand   = Hnorm_post_gc_by_kcand,
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
        list_m <- list_m + 1L
        metrics_list[[list_m]] <- metr
    } # end for g
    
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
                       sum_Hnorm_gc_by_kcand = numeric(),
                       N_garbage = integer())
    }
    
    counts_cy2 <- counts_cy %>% dplyr::rename(cnty = dplyr::all_of(county_var))
    
    out_perrecord <- counts_cy2 %>%
        dplyr::left_join(metrics_sum,  by = c("cnty", "year")) %>%
        dplyr::mutate(
            sum_Hnorm_post_K      = dplyr::coalesce(sum_Hnorm_post_K, 0),
            sum_expH_over_K       = dplyr::coalesce(sum_expH_over_K, 0),
            sum_IG_abs            = dplyr::coalesce(sum_IG_abs, 0),
            sum_IG_frac           = dplyr::coalesce(sum_IG_frac, 0),
            sum_Hnorm_gc_by_kcand = dplyr::coalesce(sum_Hnorm_gc_by_kcand, 0),
            N_garbage             = dplyr::coalesce(N_garbage, 0L),
            
            # RI uses the per-record normalized entropy over the candidate set
            RI = dplyr::if_else(N_garbage > 0,
                                1 - (sum_Hnorm_gc_by_kcand / N_garbage),
                                NA_real_),
            
            DQ_rec_expH_over_K_mean = dplyr::if_else(N_total > 0, sum_expH_over_K / N_total, NA_real_),
            DQ_rec_ig_abs_mean      = dplyr::if_else(N_total > 0, sum_IG_abs      / N_total, NA_real_),
            DQ_rec_ig_frac_mean_garbage = dplyr::if_else(N_garbage > 0, sum_IG_frac / N_garbage, NA_real_)
        ) %>%
        dplyr::select(cnty, year,
                      RI,
                      DQ_rec_expH_over_K_mean,
                      DQ_rec_ig_abs_mean,
                      DQ_rec_ig_frac_mean_garbage)
    
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
            foreman_garbage_adj = foreman_garbage * (1 - dplyr::coalesce(DQ_rec_ig_frac_mean_garbage, 0))
        ) %>%
        dplyr::rename(!!county_var := cnty)
    
    out
}


