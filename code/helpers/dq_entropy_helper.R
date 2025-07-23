# ───── dq_entropy_helper.R  ─────────────────────────
#  Fits one multinomial per garbage code per year and returns a county table
# ----------------------------------------------------------------------------

library(Matrix);  library(dplyr);  library(tidyr);  library(nnet);  library(glmnet); library(here)

collapse_illdef_cancer <- function(x) ifelse(x %in% c("C80","C96"), "C80C96", x)

bin_map <- readr::read_csv(here::here("data_raw/cause-codes/foreman-24-bins.csv"),
                           show_col_types = FALSE) |>
    transmute(icd3   = substr(trimws(ICD10), 1, 3),
              bucket = trimws(USCOD)) |>
    distinct()
bucket_lu  <- setNames(bin_map$bucket, bin_map$icd3)
bucket_unk <- "Other"

target_map <- readr::read_csv(
    here::here("data_raw/cause-codes/garbage_target_map.csv"),
    show_col_types = FALSE
) |>
    mutate(across(everything(), ~ trimws(toupper(.x)))) |>
    mutate(garbage = collapse_illdef_cancer(garbage),
           target  = collapse_illdef_cancer(target)) |>
    filter(grepl("^[A-Z][0-9]{2}$", garbage),
           grepl("^[A-Z][0-9]{2}$", target)) |>
    distinct()

valid_garbage <- collapse_illdef_cancer(c("C55","C80","C96","J80","I10","I64","N19","Y12"))

garbage_tbl <- target_map |>
    filter(garbage %in% valid_garbage) |>
    group_by(garbage) |>
    summarise(candidate_codes = list(unique(target)), .groups = "drop")

clean_icd3  <- function(x) substr(gsub("[^A-Z0-9]", "", toupper(x)), 1, 3)
clean_icd4  <- function(x) substr(gsub("[^A-Z0-9]", "", toupper(x)), 1, 4)
entropy_bits <- function(p){p<-p[p>0]; -sum(p*log2(p))}

make_features25 <- function(df){
    long <- df |>
        mutate(.row = row_number()) |>
        pivot_longer(starts_with("record_"), values_to = "code") |>
        mutate(code   = clean_icd3(code),
               bucket = dplyr::coalesce(bucket_lu[code], bucket_unk)) |>
        filter(is.na(uc3) | code != uc3) |>
        distinct(.row, bucket) |>
        mutate(flag = 1L)
    
    wide <- pivot_wider(long, id_cols = .row, names_from = bucket,
                        names_prefix = "bin_", values_from = flag,
                        values_fill = 0L)
    
    meta <- df |>
        mutate(.row = row_number()) |>
        select(.row, uc3, age_grp, sex_male)
    
    left_join(meta, wide, by = ".row") |>
        mutate(across(starts_with("bin_"), ~ tidyr::replace_na(.x, 0L))) |>
        select(-.row)
}

fit_multinom <- function(tr, te, dbg_prefix="") {
    
    # drop factor predictors with <2 levels in either tr or te
    bad_tr <- names(Filter(function(x) is.factor(x) && nlevels(x) < 2, tr))
    bad_te <- names(Filter(function(x) is.factor(x) && nlevels(x) < 2, te))
    drop_cols <- union(bad_tr, bad_te)
    
    if (length(drop_cols)) {
        tr <- tr[ , !names(tr) %in% drop_cols, drop = FALSE]
        te <- te[ , !names(te) %in% drop_cols, drop = FALSE]
    }
    
    # Set up y_tr
    y_tr <- factor(tr$uc3, levels = unique(tr$uc3))
    keep <- names(which(table(y_tr) >= 2))
    idx  <- y_tr %in% keep
    tr   <- tr[idx, , drop = FALSE]
    y_tr <- droplevels(y_tr[idx])
    
    if (nlevels(y_tr) < 2) {
        cat(dbg_prefix, " <2 levels → NA\n", sep = "")
        return(rep(NA_real_, nrow(te)))
    }
    
    tr <- select(tr, -any_of(c("age_grp", "sex_male")))
    te <- select(te, -any_of(c("age_grp", "sex_male")))
    
    # convert every non-numeric predictor to integer
    tr <- mutate(tr, across(!where(is.numeric), ~ as.integer(as.factor(.x))))
    te <- mutate(te, across(!where(is.numeric), ~ as.integer(as.factor(.x))))
    
    X_tr <- sparse.model.matrix(~ . - uc3 - 1, tr)
    X_te <- sparse.model.matrix(~ . - uc3 - 1, te)
    
    
    
    if (ncol(X_tr) == 0) {
        cat(dbg_prefix, " 0 cols → NA\n", sep = "")
        return(rep(NA_real_, nrow(te)))
    }
    
    result <- tryCatch({
        
        if (nlevels(y_tr) <= 10) {
            Y   <- nnet::class.ind(y_tr)               
            mod <- nnet::nnet(x = as.matrix(X_tr), y = Y,
                              size = 0, skip = TRUE, softmax = TRUE,
                              MaxNWts = 1e6, trace = FALSE)
            pr  <- predict(mod, as.matrix(X_te), type = "raw")
        } else {
            lambda <- 1e-6
            path   <- 10^seq(-2, log10(lambda), length.out = 50)
            mod <- suppressWarnings(
                glmnet::glmnet(X_tr, y_tr, family = "multinomial",
                               alpha = 0, lambda = path,
                               standardize = FALSE,
                               type.multinomial = "ungrouped")
            )
            pr <- predict(mod, X_te, s = lambda, type = "response")[,,1]
        }
        
        if (is.null(dim(pr))) pr <- matrix(pr, ncol = 1)
        apply(pr, 1L, entropy_bits)
        
    }, error = function(e){
        rep(NA_real_, nrow(te))
    })
    
    result
}

# ──────────────────────────────────────────────────────────────
#  Weighted-average DQ  +  zero-score for R54/R99
# ──────────────────────────────────────────────────────────────
compute_entropy_county <- function(df, county_var = "countyrs") {
    
    preds <- vector("list", nrow(garbage_tbl))
    
    # predict entropy for each modelled garbage code
    for (g in seq_len(nrow(garbage_tbl))) {
        gc     <- garbage_tbl$garbage[g]
        cand   <- garbage_tbl$candidate_codes[[g]]
        gc_vec <- if (gc == "C80C96") c("C80", "C96") else gc
        
        tr <- dplyr::filter(df, uc3 %in% cand)
        te <- dplyr::filter(df, uc4 %in% gc_vec)
        if (nrow(te) == 0 || nrow(tr) < 50) next
        
        age_to_num <- function(x) {
            xi  <- as.character(x)
            bad <- xi %in% invalid_by_var$age | is.na(xi)
            out <- suppressWarnings(as.numeric(xi))
            out[bad] <- NA_real_
            out
        }
        
        tr$age_num <- age_to_num(tr$age_years)
        te$age_num <- age_to_num(te$age_years)
        
        if (sum(!is.na(tr$age_num)) < 5) {
            tr$age_grp <- te$age_grp <- NA_character_
        } else {
            q <- quantile(tr$age_num, probs = c(0, .25, .50, .75, 1),
                          na.rm = TRUE, type = 7) |> unique()
            if (length(q) < 5) {
                q <- c(min(tr$age_num, na.rm = TRUE) - 1,
                       stats::quantile(tr$age_num, probs = c(.33, .66),
                                       na.rm = TRUE, type = 7),
                       max(tr$age_num, na.rm = TRUE) + 1)
            }
            
            tr$age_grp <- cut(tr$age_num, breaks = q,
                              labels = paste0("Q", seq_len(length(q) - 1)),
                              include.lowest = TRUE, right = TRUE)
            te$age_grp <- cut(te$age_num, breaks = q,
                              labels = paste0("Q", seq_len(length(q) - 1)),
                              include.lowest = TRUE, right = TRUE)
        }
        trf <- make_features25(tr) |>
            mutate(age_grp  = factor(age_grp, levels = paste0("Q", 1:4)),
                   sex_male = factor(sex_male))
        tef <- make_features25(te) |>
            mutate(age_grp  = factor(age_grp, levels = paste0("Q", 1:4)),
                   sex_male = factor(sex_male))
        
        missing_bins <- setdiff(names(trf), names(tef))
        tef[missing_bins[grepl("^bin_", missing_bins)]] <- 0L
        
        H <- fit_multinom(trf, tef, dbg_prefix = paste0("▶ ", gc, " "))
        if (sum(is.na(H)) > 0) {
            cat("⚠ GC", gc, ":",
                sum(is.na(H)), "NA predictions out of", length(H), "\n")
            
            na_ranum <- te$ranum[is.na(H)]
            cat("  Affected ranum (sample):", head(na_ranum, 5), "\n")
            
            # Save to global (optional)
            assign(paste0("NA_H_", gc), te[is.na(H), c("ranum", county_var)], envir = .GlobalEnv)
        }
        
        
        preds[[g]] <- tibble(
            ranum   = te$ranum,
            garbage = gc,
            k       = length(cand),
            H_B     = H
        )
    }
    
    pred_bind <- dplyr::bind_rows(preds) |>
        dplyr::left_join(df[, c("ranum", county_var)], by = "ranum") |>
        mutate(DQ_row = 1 - H_B / log2(k))
    
    # county-level means and counts for modelled GCs 
    by_gc_long <- pred_bind |>
        group_by(.data[[county_var]], garbage) |>
        summarise(
            DQ_gc = mean(DQ_row, na.rm = TRUE),
            n_gc  = n(),
            .groups = "drop"
        )
    
    # add R54 / R99 with DQ = 0
    r_tbl <- df |>
        filter(uc3 %in% c("R54", "R99")) |>
        group_by(.data[[county_var]]) |>
        summarise(
            garbage = "R54R99",
            DQ_gc   = 0,
            n_gc    = n(),
            .groups = "drop"
        )
    by_gc_long <- bind_rows(by_gc_long, r_tbl)
    
    # weights & overall (sum w*DQ_gc)
    by_gc_long <- by_gc_long |>
        group_by(.data[[county_var]]) |>
        mutate(weight = n_gc / sum(n_gc)) 

    overall <- by_gc_long |>
        summarise(
            DQ_overall = sum(DQ_gc * weight, na.rm = TRUE),
            .groups    = "drop"
        )
    
    # wide table with per-GC DQs
    by_gc_wide <- by_gc_long |>
        select(-n_gc, -weight) |>
        pivot_wider(names_from  = garbage,
                    values_from = DQ_gc,
                    names_prefix = "DQ_")     # includes DQ_R54R99
    
    # return
    # ───────────────────────────────────────────────────────────────
    # ▶▶ Diagnostics: entropy prediction coverage
    # ───────────────────────────────────────────────────────────────
    n_gc_total   <- nrow(garbage_tbl)
    n_gc_pred    <- sum(lengths(preds) > 0)
    n_gc_skipped <- n_gc_total - n_gc_pred
    
    if (n_gc_skipped > 0) {
        cat("⚠ compute_entropy_county: Skipped", n_gc_skipped,
            "of", n_gc_total, "garbage codes (no test/train rows or failed)\n")
    }
    
    missing_entropy <- pred_bind |> filter(is.na(H_B))
    if (nrow(missing_entropy) > 0) {
        cat("⚠ compute_entropy_county: Found", nrow(missing_entropy),
            "rows with missing entropy predictions\n")
    }
    
    if (nrow(overall) == 0 || any(is.na(overall$DQ_overall))) {
        cat("🚨 compute_entropy_county: DQ_overall missing or not computed\n")
    }
    
    left_join(overall, by_gc_wide, by = county_var)
}