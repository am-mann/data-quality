# ───────────────── dq_entropy_helper.R — Foreman 24/25-bin DQ (county×year) ─────────────────
# Inputs:
#   - ds: data frame with 'ucod', county key, optional 'age_years' or 'ager27',
#         and contributing causes in record_1..record_20 (or record1..record20)
#   - county_var: name of the county column in ds (e.g., "county_ihme")
#   - dict_dir: directory containing:
#       * foreman-icd10-mapping.csv  (ICD → USCOD *codes*, NOT names)
#       * foreman-table2-map.csv     (wide: target_cause + G_1..G_9 TRUE/FALSE)
#
# Output:
#   tibble: county_var, year, total, DQ_entropy, DQ_K, DQ_overall, DQ_expH_over_K,
#           DQ_rec_entropy_mean (slides-consistent), DQ_rec_expH_over_K_mean,
#           DQ_rec_ig_abs_mean, DQ_rec_ig_frac_mean_garbage,
#           foreman_garbage, foreman_garbage_adj
#
# Notes:
#   Robust to: differing header orders, record1 vs record_1, missing dimnames from glmnet
#   Glmnet gating: needs ≥30 training rows, ≥2 classes, and ≥8 per class; else proportional fallback
#   ICD mapping tries 4-char first, then 3-char
#
# Methodology: Foreman KJ, Naghavi M, Ezzati M (2016), Population Health Metrics 14:14.
# ─────────────────────────────────────────────────────────────────────────────────────────────

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
norm_code     <- function(x) { x %>% tolower() %>% stringr::str_replace_all("[^a-z0-9]", "") } # "B_3_1"→"b31"

# ---------- robust readers tailored to your CSVs ----------
.read_icd_map <- function(path) {
    df  <- readr::read_csv(path, show_col_types = FALSE, locale = readr::locale(encoding = "UTF-8"))
    nms <- tolower(names(df))
    
    # Strict preference for exact 'uscod' (code), not 'uscod name'
    col_icd <- NA_character_
    icd_candidates <- c("icd10","code","icd","icd10 code","icd10_name","icd10 name")
    m <- match(TRUE, nms %in% icd_candidates)
    if (!is.na(m)) col_icd <- names(df)[m]
    if (is.na(col_icd) && "icd10" %in% names(df)) col_icd <- "ICD10"
    if (is.na(col_icd)) stop("ICD column not found in ", basename(path), ". Saw: ", paste(names(df), collapse=", "))
    
    if ("uscod" %in% nms) {
        col_uscod <- names(df)[which(nms == "uscod")][1]
    } else {
        cand <- names(df)[nms %in% c("foreman","foreman bin","foreman_bin","uscode")]
        if (length(cand) == 0L) {
            stop("USCOD (code) column not found in ", basename(path),
                 ". Include a column named 'USCOD' (codes like A_1, B_3_1, ...).")
        }
        col_uscod <- cand[1]
    }
    
    df %>%
        transmute(
            icd = canonical_icd(.data[[col_icd]]),
            bin = norm_code(.data[[col_uscod]])
        ) %>%
        filter(icd != "", bin != "") %>%
        distinct(icd, .keep_all = TRUE)
}

.read_table2_map <- function(path) {
    df <- readr::read_csv(path, show_col_types = FALSE, locale = readr::locale(encoding = "UTF-8"))
    if (!"target_cause" %in% names(df)) {
        stop("Expected 'target_cause' column in ", basename(path),
             ". Saw: ", paste(names(df), collapse = ", "))
    }
    g_cols <- grep("^G_[0-9]+$", names(df), value = TRUE, ignore.case = FALSE)
    if (!length(g_cols)) stop("Expected garbage flag columns like G_1..G_9 in ", basename(path))
    
    # Coerce truthy strings to logical
    df[g_cols] <- lapply(df[g_cols], function(x) {
        if (is.logical(x)) return(x)
        y <- tolower(as.character(x))
        truthy <- y %in% c("1","true","t","y","yes","x","✓","check","checked")
        out <- rep(FALSE, length(y))
        out[truthy] <- TRUE
        out
    })
    
    df %>%
        mutate(target = norm_code(target_cause)) %>%
        tidyr::pivot_longer(all_of(g_cols), names_to = "garbage_raw", values_to = "flag") %>%
        filter(flag) %>%
        transmute(
            target  = norm_code(target),
            garbage = norm_code(garbage_raw)
        ) %>%
        filter(target != "", garbage != "") %>%
        distinct(garbage, target)
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

# Sparse one-hot for contributing bins present per row
make_cc_matrix <- function(cc_bins_by_row, feat_bins) {
    if (!length(cc_bins_by_row)) {
        Matrix(0, 0, length(feat_bins), sparse = TRUE, dimnames = list(NULL, feat_bins))
    } else {
        j_lookup <- setNames(seq_along(feat_bins), feat_bins)
        lens <- vapply(cc_bins_by_row, length, 1L)
        if (!length(lens)) return(Matrix(0, 0, length(feat_bins), sparse = TRUE, dimnames = list(NULL, feat_bins)))
        ii <- rep.int(seq_along(cc_bins_by_row), lens)
        jj <- unlist(lapply(cc_bins_by_row, function(v) unname(j_lookup[v])), use.names = FALSE)
        keep <- !is.na(jj)
        Matrix::sparseMatrix(i = ii[keep], j = jj[keep], x = 1,
                             dims = c(length(cc_bins_by_row), length(feat_bins)),
                             dimnames = list(NULL, feat_bins))
    }
}

# Entropy
.shannon_H <- function(p) { p <- p[p > 0]; if (!length(p)) 0 else -sum(p * log(p)) }

# Robust glmnet prediction: fixes dim/labels/row-mismatch
.safe_multinomial_probs <- function(fit, X_te, classes) {
    arr <- tryCatch(predict(fit, newx = X_te, type = "response", s = "lambda.min"),
                    error = function(e) NULL)
    if (is.null(arr)) return(NULL)
    
    # 3D (n × K × L) → take first lambda; else coerce to matrix
    if (length(dim(arr)) == 3L) {
        probs <- arr[, , 1, drop = FALSE]
    } else {
        probs <- arr
    }
    probs <- as.matrix(probs)
    
    # Ensure we have column names (class labels)
    if (is.null(colnames(probs)) || anyNA(colnames(probs)) || ncol(probs) == 0L) {
        dn <- dimnames(arr)
        if (!is.null(dn) && length(dn) >= 2L && !is.null(dn[[2]]) && length(dn[[2]]) == ncol(probs)) {
            colnames(probs) <- dn[[2]]
        } else if (!is.null(fit$glmnet.fit$classnames) &&
                   length(fit$glmnet.fit$classnames) == ncol(probs)) {
            colnames(probs) <- fit$glmnet.fit$classnames
        } else {
            colnames(probs) <- paste0("cls", seq_len(max(1L, ncol(probs))))
        }
    }
    
    n_te <- nrow(X_te)
    # If glmnet gave fewer/more rows than X_te, rebuild a matrix with matching shape
    if (nrow(probs) != n_te) {
        k <- max(1L, ncol(probs))
        vals <- as.numeric(probs)
        if (length(vals) != n_te * k) vals <- rep_len(vals, n_te * k)
        probs <- matrix(vals, nrow = n_te, ncol = k, byrow = TRUE,
                        dimnames = list(NULL, colnames(probs)))
    }
    
    # Align to requested class order, pad any missing classes with zeros
    present <- intersect(classes, colnames(probs))
    base <- if (length(present)) probs[, present, drop = FALSE]
    else matrix(0, nrow = n_te, ncol = 0, dimnames = list(NULL, character()))
    missing <- setdiff(classes, colnames(probs))
    if (length(missing)) {
        base <- cbind(base,
                      matrix(0, nrow = n_te, ncol = length(missing),
                             dimnames = list(NULL, missing)))
    }
    base <- base[, classes, drop = FALSE]
    
    # Row-normalize
    rs <- rowSums(base); rs[rs == 0] <- 1
    base / rs
}

# =============== MAIN (memory-safe, streaming) ===============
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
    
    # ---- read dictionaries ----
    icd_map  <- .read_icd_map(icd_map_path)       # ICD → USCOD codes (normalized to 'bin')
    t2_map   <- .read_table2_map(code_map_path)   # long: (garbage, target)
    icd_lu   <- .make_icd_lookup(icd_map)
    
    garbage_bins     <- sort(unique(t2_map$garbage))           # e.g., g1..g9
    all_bins_in_map  <- sort(unique(icd_map$bin))              # e.g., a1.., b11.., c11.., etc.
    valid_bins       <- setdiff(all_bins_in_map, garbage_bins) # non-garbage only
    K <- length(valid_bins)
    if (K < 22) warning("Foreman map defines a small set of valid bins (K=", K, "). Check dictionaries.")
    
    # Precompute candidate targets per garbage bin
    targets_by_g <- split(unique(t2_map$target), t2_map$garbage)
    
    # ---- standardize DS ----
    if (!county_var %in% names(ds)) {
        if ("county" %in% names(ds)) {
            ds <- dplyr::rename(ds, !!county_var := .data$county)
        } else {
            stop("county_var '", county_var, "' not found in ds (and 'county' not present).")
        }
    }
    if (!"ucod" %in% names(ds)) stop("ds must contain 'ucod'")
    if (!"year" %in% names(ds) || all(is.na(ds$year))) {
        stop("ds must contain 'year' (numeric) for county×year aggregation.")
    }
    
    # Accept record1..record20 or record_1..record_20
    rec_names_underscore <- paste0("record_", 1:20)
    if (!any(rec_names_underscore %in% names(ds)) && any(grepl("^record[0-9]+$", names(ds)))) {
        ds <- ds %>% rename_with(~ sub("^record([0-9]+)$", "record_\\1", .x), .cols = matches("^record[0-9]+$"))
    }
    rec_cols <- intersect(rec_names_underscore, names(ds))
    
    keep <- c("ucod", county_var, "year", "age_years", "ager27", rec_cols)
    keep <- intersect(keep, names(ds))
    dat  <- ds %>% dplyr::select(dplyr::all_of(keep)) %>%
        dplyr::mutate(
            uc4  = icd4(ucod),
            uc3  = icd3(ucod),
            ubin = map_icd_to_bin(ucod, icd_lu)
        )
    
    n_rows <- nrow(dat)
    
    # ---------- Light-weight long CC table (for filtering only) ----------
    cc_long <- if (length(rec_cols) == 0) {
        tibble::tibble(.row = integer(), cc_bin = character())
    } else {
        dat %>%
            dplyr::mutate(.row = dplyr::row_number()) %>%
            dplyr::select(.row, dplyr::all_of(rec_cols), uc4, uc3) %>%
            tidyr::pivot_longer(cols = dplyr::all_of(rec_cols), names_to = "rec", values_to = "cc_raw") %>%
            dplyr::mutate(cc_icd4 = icd4(cc_raw)) %>%
            dplyr::filter(!(substr(cc_icd4, 1, 4) == uc4 | substr(cc_icd4, 1, 3) == uc3)) %>%
            dplyr::transmute(.row, cc_bin = map_icd_to_bin(cc_icd4, icd_lu)) %>%
            dplyr::filter(!is.na(cc_bin) & cc_bin != "")
    }
    
    # Identify garbage UCOD rows
    row_is_garbage <- !is.na(dat$ubin) & dat$ubin %in% garbage_bins
    
    # ---------- count base (non-garbage) UCODs directly ----------
    base_counts_agg <- dat %>%
        dplyr::filter(!row_is_garbage & !is.na(ubin) & ubin %in% valid_bins) %>%
        dplyr::count(.data[[county_var]], year, bin = ubin, name = "k", .drop = FALSE)
    
    # total certificates per county×year (for per-record means)
    counts_cy <- dat %>%
        dplyr::count(.data[[county_var]], year, name = "N_total", .drop = FALSE)
    
    # holders for aggregated expected counts from garbage rows and per-record metric sums
    prob_counts_list <- list()
    metrics_list     <- list()
    list_i <- 0L
    list_m <- 0L
    
    # local helper: build features ONLY for keep_idx rows
    build_X_for_rows <- function(keep_idx) {
        # CC sparse (restrict to keep_idx and valid_bins actually present)
        if (nrow(cc_long) == 0 || length(keep_idx) == 0) {
            X_cc <- Matrix::Matrix(0, nrow = length(keep_idx), ncol = 0, sparse = TRUE)
            colnames(X_cc) <- character(0)
        } else {
            cc_sub <- cc_long %>%
                dplyr::filter(.row %in% keep_idx, cc_bin %in% valid_bins) %>%
                dplyr::distinct(.row, cc_bin)
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
        
        # Age one-hot for keep rows
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
    
    # robust prediction wrapper for glmnet
    .predict_glmnet_probs <- function(fit, X_te, classes, s = "lambda.1se") {
        arr <- tryCatch(predict(fit, newx = X_te, type = "response", s = s),
                        error = function(e) NULL)
        if (is.null(arr)) return(NULL)
        if (length(dim(arr)) == 3L) probs <- arr[, , 1, drop = FALSE] else probs <- arr
        probs <- as.matrix(probs)
        if (is.null(colnames(probs)) || anyNA(colnames(probs)) || ncol(probs) == 0L) {
            if (!is.null(fit$glmnet.fit$classnames) &&
                length(fit$glmnet.fit$classnames) == ncol(probs)) {
                colnames(probs) <- fit$glmnet.fit$classnames
            } else {
                colnames(probs) <- paste0("cls", seq_len(max(1L, ncol(probs))))
            }
        }
        present <- intersect(classes, colnames(probs))
        base <- if (length(present)) probs[, present, drop = FALSE]
        else matrix(0, nrow = nrow(X_te), ncol = 0, dimnames = list(NULL, character()))
        missing <- setdiff(classes, colnames(probs))
        if (length(missing)) {
            base <- cbind(base, matrix(0, nrow = nrow(X_te), ncol = length(missing),
                                       dimnames = list(NULL, missing)))
        }
        base <- base[, classes, drop = FALSE]
        rs <- rowSums(base); rs[rs == 0] <- 1
        base / rs
    }
    
    # ---------- loop garbage groups; fit on-demand; aggregate ----------
    prob_counts_list <- list()
    metrics_list     <- list()
    list_i <- 0L
    list_m <- 0L
    
    for (g in garbage_bins) {
        targets_g <- unique(targets_by_g[[g]])
        if (!length(targets_g)) next
        
        # rows where CC contains g (training filter)
        if (nrow(cc_long)) {
            rows_with_g_cc <- unique(cc_long$.row[cc_long$cc_bin == g])
        } else rows_with_g_cc <- integer(0)
        
        tr_idx <- intersect(rows_with_g_cc, which(!is.na(dat$ubin) & dat$ubin %in% targets_g))
        te_idx <- intersect(which(row_is_garbage), which(dat$ubin == g))
        if (!length(te_idx)) next
        
        model_ok <- FALSE; probs <- NULL
        classes  <- sort(targets_g)
        
        if (.has_glmnet() && length(tr_idx) >= 30) {
            cls_counts <- table(dat$ubin[tr_idx])
            if (length(cls_counts) >= 2 && !any(cls_counts < 8)) {
                # downsample per class (cap)
                cap <- as.integer(getOption("DQ_MAX_TRAIN_PER_CLASS", 8000L))
                take <- unlist(lapply(split(tr_idx, dat$ubin[tr_idx]), function(ix) {
                    if (length(ix) > cap) sample(ix, cap) else ix
                }), use.names = FALSE)
                tr_idx2 <- sort(take)
                
                keep_idx <- sort(unique(c(tr_idx2, te_idx)))
                X_local  <- build_X_for_rows(keep_idx)
                nz       <- if (ncol(X_local)) which(Matrix::colSums(X_local[match(tr_idx2, keep_idx), , drop = FALSE]) > 0) else integer(0)
                
                if (length(nz)) {
                    X_tr <- X_local[match(tr_idx2, keep_idx), nz, drop = FALSE]
                    X_te <- X_local[match(te_idx,  keep_idx), nz, drop = FALSE]
                    y_tr <- factor(dat$ubin[tr_idx2], levels = classes)
                    
                    fit <- tryCatch(glmnet::cv.glmnet(
                        X_tr, y_tr, family = "multinomial", alpha = 0,
                        type.multinomial = "ungrouped", nfolds = 3, parallel = FALSE
                    ), error = function(e) NULL)
                    
                    if (!is.null(fit)) {
                        probs <- .predict_glmnet_probs(fit, X_te, classes, s = "lambda.1se")
                        model_ok <- !is.null(probs)
                    }
                }
            }
        }
        
        if (!model_ok) {
            # proportional or uniform fallback
            if (length(tr_idx)) {
                tbl <- table(dat$ubin[tr_idx])
                w   <- as.numeric(tbl) / sum(tbl)
                probs <- matrix(rep(w, each = length(te_idx)),
                                nrow = length(te_idx), ncol = length(w),
                                dimnames = list(NULL, names(tbl)))
                # align to classes
                tmp <- matrix(0, nrow = nrow(probs), ncol = length(classes),
                              dimnames = list(NULL, classes))
                tmp[, colnames(probs)] <- probs
                probs <- tmp
            } else {
                probs <- matrix(1/length(classes), nrow = length(te_idx), ncol = length(classes),
                                dimnames = list(NULL, classes))
            }
        }
        
        # ---- aggregate expected counts for this group's garbage rows ----
        te_meta <- dat[te_idx, c(county_var, "year")]
        colnames(te_meta) <- c("cnty", "year")
        cls_summaries <- lapply(seq_along(classes), function(j) {
            tibble::tibble(cnty = te_meta$cnty, year = te_meta$year, k = probs[, j]) %>%
                dplyr::group_by(cnty, year) %>%
                dplyr::summarise(k = sum(k), .groups = "drop") %>%
                dplyr::mutate(bin = classes[j])
        })
        list_i <- list_i + 1L
        prob_counts_list[[list_i]] <- dplyr::bind_rows(cls_summaries)
        
        # ---- per-record metric partials (garbage rows only) ----
        eps <- .Machine$double.eps
        H_post <- rowSums(-probs * log(pmax(probs, eps)))   # posterior entropy
        k_cand <- length(classes)
        H0     <- log(k_cand)                               # prior entropy (uniform over candidates)
        
        # Slides-consistent normalization: H_post / log(k_cand) (only for garbage-coded)
        Hnorm_post_gc_by_kcand <- if (k_cand > 1) H_post / log(k_cand) else rep(0, length(H_post))
        
        # Additional diagnostics (kept)
        Hnorm_post_K      <- H_post / log(K)
        expH_over_K_vec   <- exp(H_post) / K
        IG_abs_vec        <- H0 - H_post
        IG_frac_vec       <- if (H0 > 0) (H0 - H_post) / H0 else rep(NA_real_, length(H_post))
        
        metr <- tibble::tibble(
            cnty = te_meta$cnty, year = te_meta$year,
            sum_Hnorm_post_K        = Hnorm_post_K,
            sum_expH_over_K         = expH_over_K_vec,
            sum_IG_abs              = IG_abs_vec,
            sum_IG_frac             = IG_frac_vec,
            sum_Hnorm_gc_by_kcand   = Hnorm_post_gc_by_kcand,  # NEW: slides-consistent term
            one = 1L
        ) %>%
            dplyr::group_by(cnty, year) %>%
            dplyr::summarise(
                sum_Hnorm_post_K      = sum(sum_Hnorm_post_K,      na.rm = TRUE),
                sum_expH_over_K       = sum(sum_expH_over_K,       na.rm = TRUE),
                sum_IG_abs            = sum(sum_IG_abs,            na.rm = TRUE),
                sum_IG_frac           = sum(sum_IG_frac,           na.rm = TRUE),
                sum_Hnorm_gc_by_kcand = sum(sum_Hnorm_gc_by_kcand, na.rm = TRUE),  # NEW
                N_garbage             = sum(one),
                .groups = "drop"
            )
        list_m <- list_m + 1L
        metrics_list[[list_m]] <- metr
    } # end for g
    
    # ---------- combine aggregated pieces ----------
    prob_counts_all <- if (length(prob_counts_list)) dplyr::bind_rows(prob_counts_list) else
        tibble::tibble(cnty = character(), year = numeric(), k = numeric(), bin = character())
    
    base_counts_all <- base_counts_agg %>%
        dplyr::rename(cnty = dplyr::all_of(county_var))
    
    # total expected counts per county×year×bin
    cty_year_bin <- dplyr::bind_rows(
        base_counts_all %>% dplyr::select(cnty, year, bin, k),
        prob_counts_all  %>% dplyr::select(cnty, year, bin, k)
    ) %>%
        dplyr::group_by(cnty, year, bin) %>%
        dplyr::summarise(k = sum(k), .groups = "drop")
    
    # Foreman aggregate DQ
    out_aggregate <- cty_year_bin %>%
        dplyr::group_by(cnty, year) %>%
        dplyr::summarise(
            total          = sum(k),
            DQ_entropy     = { p <- k / sum(k); .shannon_H(p) },
            DQ_K           = K,
            DQ_overall     = ifelse(DQ_K > 0, DQ_entropy / log(DQ_K), NA_real_),
            DQ_expH_over_K = exp(DQ_entropy) / DQ_K,
            .groups = "drop"
        )
    
    # per-record metrics (streamed): sums over garbage rows + N_garbage
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
        dplyr::left_join(metrics_sum, by = c("cnty", "year")) %>%
        dplyr::mutate(
            sum_Hnorm_post_K      = dplyr::coalesce(sum_Hnorm_post_K, 0),
            sum_expH_over_K       = dplyr::coalesce(sum_expH_over_K, 0),
            sum_IG_abs            = dplyr::coalesce(sum_IG_abs, 0),
            sum_IG_frac           = dplyr::coalesce(sum_IG_frac, 0),
            sum_Hnorm_gc_by_kcand = dplyr::coalesce(sum_Hnorm_gc_by_kcand, 0),
            N_garbage             = dplyr::coalesce(N_garbage, 0L),
            
            # === Slides-consistent: mean over garbage-coded of H_r / log(k_cand) ===
            DQ_rec_entropy_mean = dplyr::if_else(N_garbage > 0,
                                                 sum_Hnorm_gc_by_kcand / N_garbage,
                                                 NA_real_),
            
            # Additional diagnostics (global-K based, as before)
            DQ_rec_expH_over_K_mean = dplyr::if_else(N_total > 0, sum_expH_over_K / N_total, NA_real_),
            DQ_rec_ig_abs_mean      = dplyr::if_else(N_total > 0, sum_IG_abs      / N_total, NA_real_),
            DQ_rec_ig_frac_mean_garbage = dplyr::if_else(N_garbage > 0, sum_IG_frac / N_garbage, NA_real_)
        ) %>%
        dplyr::select(cnty, year,
                      DQ_rec_entropy_mean,
                      DQ_rec_expH_over_K_mean,
                      DQ_rec_ig_abs_mean,
                      DQ_rec_ig_frac_mean_garbage)
    
    # ---------- Foreman garbage shares ----------
    fg_share <- dat %>%
        dplyr::mutate(is_fg = !is.na(ubin) & ubin %in% garbage_bins) %>%
        dplyr::count(.data[[county_var]], year, is_fg, name = "k") %>%
        tidyr::pivot_wider(names_from = is_fg, values_from = k, values_fill = 0) %>%
        dplyr::transmute(
            cnty = .data[[county_var]], year,
            N_total = `FALSE` + `TRUE`,
            foreman_garbage = dplyr::if_else(N_total > 0, `TRUE` / N_total, NA_real_)
        )
    
    # ---------- final join + rename county key back ----------
    out <- out_aggregate %>%
        dplyr::left_join(out_perrecord, by = c("cnty", "year")) %>%
        dplyr::left_join(fg_share %>% dplyr::select(cnty, year, foreman_garbage),
                         by = c("cnty", "year")) %>%
        dplyr::mutate(
            # Adjusted garbage share: discount by information recovered from garbage-coded deaths
            foreman_garbage_adj = foreman_garbage * (1 - dplyr::coalesce(DQ_rec_ig_frac_mean_garbage, 0))
        ) %>%
        dplyr::rename(!!county_var := cnty)
    
    out
}

