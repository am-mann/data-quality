# ————————————————————————————————————————————————————————————————
# Date: July 18, 2025
# Code to calculate averages for various data quality metrics by county.
#
# Assumes mortality files in format mortXXXX.parquet and variable names
# "year", "marstat", "placdth", "educ2003", "mandeath", "age", "sex", "race",  "ucod", "record1", ..., "record20"
#
# Requires the following files to run:
# - mortality data files from 1999-2023 (called mortXXX.parquet)
# - overdose-detail-codes-v3.csv containing the required contributing cause codes for overdose deaths
# - light_garbage_codes.csv containing list of contributing cause codes that lack detail
# - accident-detail-codes-v3.csv containing the required contributing cause codes for accident deaths
# - gbd_garbage_codes_without_overdose.csv
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
#   DQ_rec_entropy_mean, DQ_rec_expH_over_K_mean, DQ_rec_ig_abs_mean,
#   DQ_rec_ig_frac_mean_garbage, foreman_garbage, foreman_garbage_adj
# ─────────────────────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(readr); library(stringr)
    library(purrr); library(Matrix); library(here); library(binom); library(arrow)
})

PARALLELIZE <- TRUE
N_WORKERS   <- 8
WRITE_ONE_CSV     <- TRUE 
out_csv           <- here("data", "county_year_quality_metrics.csv.gz")

options(
    DQ_MAX_TRAIN_PER_CLASS = 3000L,  
    DQ_NLAMBDA = 40,               
    DQ_CV_FOLDS = 3                
)


source(here("code/helpers", "dq_entropy_helper.R"))

# ---------- tiny helpers ----------
.has_glmnet <- function() requireNamespace("glmnet", quietly = TRUE)
canonical_icd <- function(x) stringr::str_remove_all(stringr::str_to_upper(x), "[^A-Z0-9]")
clean_icd     <- function(x) stringr::str_remove_all(stringr::str_to_upper(x), "[^A-Z0-9\\.]")
num <- function(x) if (bit64::is.integer64(x)) as.numeric(x) else as.numeric(x)

wilson_lower <- function(k, n) {
    k <- num(k); n <- num(n); out <- rep(NA_real_, length(k))
    ok <- n > 0 & !is.na(n); if (any(ok)) out[ok] <- binom.confint(k[ok], n[ok], methods = "wilson")$lower
    out
}
wilson_upper <- function(k, n) {
    k <- num(k); n <- num(n); out <- rep(NA_real_, length(k))
    ok <- n > 0 & !is.na(n); if (any(ok)) out[ok] <- binom.confint(k[ok], n[ok], methods = "wilson")$upper
    out
}

# Accept record1..record20 or record_1..record_20
rename_record_cols <- function(df) {
    if (!any(grepl("^record_", names(df))) && any(grepl("^record[0-9]+$", names(df)))) {
        dplyr::rename_with(df, ~ sub("^record([0-9]+)$", "record_\\1", .x), .cols = dplyr::matches("^record[0-9]+$"))
    } else df
}

# -------- paths --------
parquet_dir <- if (dir.exists(here("data_private", "mcod"))) {
    here("data_private", "mcod")
} else if (dir.exists(here("data_private", "mcod_sample"))) {
    here("data_private", "mcod_sample")
} else stop("Neither data_private/mcod nor data_private/mcod_sample exists.")

dictionary_dir <- here("data_raw", "cause-codes")
county_var     <- "county_ihme"
years_wanted   <- NULL

key_demo_vars <- c("marstat","placdth","educ","mandeath","age","sex","race","racer40","educ2003")
sec_cols <- paste0("record_", 1:20)

invalid_by_var <- list(
    marstat   = c("", "not_stated", "unknown", NA),
    placdth   = c("", "9", "U", "Unknown", NA),
    educ2003  = c("", "9", "99", "999", "9999", "0", NA),
    educ      = c("", "9", "99", "999", "9999", "0", NA),
    mandeath  = c("", "9", "7", "Unknown", NA),
    age       = c("", "999", "9999", "0000", NA),
    sex       = c("", "U", NA),
    race      = c("", "9", "Unknown", NA),
    racer40   = c("", "9", "Unknown", NA)
)

# ---------- lookups (read once) ----------
lookup_garbage <- read_csv(file.path(dictionary_dir, "gbd_garbage_codes_without_overdose.csv"),
                           show_col_types = FALSE) |>
    transmute(icd10 = canonical_icd(icd10), gbd_severity = as.integer(gbd_severity)) |>
    distinct(icd10, .keep_all = TRUE)

light_garbage <- read_csv(file.path(dictionary_dir, "light_garbage_codes.csv"),
                          show_col_types = FALSE)$code |>
    canonical_icd()

overdose_ucod <- read_csv(file.path(dictionary_dir, "overdose-detail-codes-v3.csv"),
                          show_col_types = FALSE)$ucod |>
    canonical_icd()

accident_detail_lookup <- readr::read_csv(
    file.path(dictionary_dir, "accident-detail-codes-v3.csv"),
    show_col_types = FALSE
) |>
    tidyr::pivot_longer(starts_with("detail_"), values_to = "detail", values_drop_na = TRUE) |>
    transmute(ucod = canonical_icd(stringr::str_trim(ucod)),
              detail = canonical_icd(stringr::str_trim(detail))) |>
    filter(detail != "") |>
    distinct()

accident_ucod_set   <- unique(accident_detail_lookup$ucod)
accident_detail_set <- unique(accident_detail_lookup$detail)
acc_req_pair_keys   <- paste0(accident_detail_lookup$ucod, "|", accident_detail_lookup$detail)

# ─────────────────────────────────────────────────────────────────────────────
#  summarise everything for one year (MEMORY-SAFE)
# ─────────────────────────────────────────────────────────────────────────────
summarise_year_file <- function(file) {
    on.exit({gc()}, add = TRUE)
    
    this_year <- as.integer(stringr::str_extract(basename(file), "\\d{4}"))
    message("Processing ", this_year)
    
    # read only needed columns (avoid pulling everything)
    cols_needed <- c("ucod", county_var, key_demo_vars, sec_cols, "ager27", "sex", "age_years", "year")
    tab0 <- arrow::read_parquet(file, as_data_frame = FALSE)
    available_cols <- names(tab0)              # or: tab0$schema$names
    keep_cols <- intersect(cols_needed, available_cols)
    rm(tab0); gc(FALSE)
    
    ds <- arrow::read_parquet(file, as_data_frame = TRUE, col_select = keep_cols) |>
        rename_record_cols() |>
        mutate(
            ranum    = dplyr::row_number(),
            icd10    = canonical_icd(ucod),
            uc4      = substr(clean_icd(ucod), 1, 4),
            uc3      = substr(clean_icd(ucod), 1, 3),
            sex_male = if_else(sex == "M", 1L, 0L)
        )
    
    # Ensure demo columns & year
    for (v in setdiff(key_demo_vars, names(ds))) ds[[v]] <- NA_character_
    if (!"year" %in% names(ds) || all(is.na(ds$year))) ds$year <- this_year
    if (this_year < 2004) ds$educ2003 <- NA_character_
    demo_vars <- if (this_year < 2004) setdiff(key_demo_vars, "educ2003") else key_demo_vars
    if (this_year %in% 2021:2022) demo_vars <- union(setdiff(demo_vars, "race"), "racer40")
    
    # flags that don't need tall tables
    ds <- ds |>
        left_join(lookup_garbage, by = c("icd10")) |>
        mutate(
            needs_detail = icd10 %in% accident_ucod_set,
            is_overdose  = icd10 %in% overdose_ucod
        )
    
    # ========== SINGLE-PASS over 20 record_* columns (no tall pivots) ==========
    n <- nrow(ds)
    present_sec <- intersect(sec_cols, names(ds))
    total_contrib      <- integer(n)
    light_contrib      <- integer(n)
    has_any_T          <- logical(n)
    has_specific_T     <- logical(n)
    has_required_detail<- logical(n)
    
    # helper closures (fast startsWith for T-codes)
    starts_with <- function(x, pat) substr(x, 1, nchar(pat)) == pat
    
    for (col in present_sec) {
        v <- canonical_icd(ds[[col]])
        not_blank <- !is.na(v) & nzchar(v)
        
        # counts
        if (any(not_blank)) {
            total_contrib[not_blank] <- total_contrib[not_blank] + 1L
            light_contrib[not_blank] <- light_contrib[not_blank] + as.integer(v[not_blank] %in% light_garbage)
            
            # overdose details (T-codes) → "specific" if NOT T50.9 or T40.9
            tmask <- not_blank & startsWith(v, "T")
            if (any(tmask)) {
                has_any_T[tmask] <- TRUE
                spec <- !(starts_with(v[tmask], "T509") | starts_with(v[tmask], "T409"))
                if (any(spec)) has_specific_T[which(tmask)[spec]] <- TRUE
            }
            
            # required accident detail (match ucod/detail pairs)
            if (any(ds$needs_detail)) {
                keys <- paste0(ds$icd10, "|", v)
                hit  <- not_blank & (keys %in% acc_req_pair_keys)
                if (any(hit)) has_required_detail[hit] <- TRUE
            }
        }
        # drop temporaries early
        rm(v, not_blank); gc(FALSE)
    }
    
    unspecific_drug <- has_any_T & !has_specific_T
    
    contrib_counts <- tibble(
        ranum = ds$ranum,
        total_contrib = total_contrib,
        light_contrib = light_contrib
    )
    unspec_tbl <- tibble(ranum = ds$ranum, unspecific_drug = unspecific_drug)
    
    # ========== Merge flags & compute county-year metrics ==========
    cert_tbl <- ds |>
        left_join(contrib_counts, by = "ranum") |>
        left_join(unspec_tbl, by = "ranum") |>
        mutate(
            total_contrib  = tidyr::replace_na(total_contrib, 0L),
            light_contrib  = tidyr::replace_na(light_contrib, 0L),
            unspecific_drug= tidyr::replace_na(unspecific_drug, TRUE),
            has_required_detail = has_required_detail,
            is_garbage = dplyr::case_when(
                is.na(gbd_severity) ~ FALSE,
                !stringr::str_starts(icd10, "X") ~ TRUE,
                stringr::str_starts(icd10, "X") & light_contrib == 0 ~ TRUE,
                TRUE ~ FALSE
            )
        )
    
    demo_vars_exc_race <- setdiff(demo_vars, c("race", "racer40", "mandeath"))
    valid_exc_race <- sapply(demo_vars_exc_race, function(v) {
        vals <- cert_tbl[[v]]; inv <- invalid_by_var[[v]]; !(vals %in% inv)
    })
    cert_tbl$complete_all_exc_race <- rowSums(valid_exc_race) == length(demo_vars_exc_race)
    
    valid_mat <- sapply(demo_vars, function(v) {
        vals <- cert_tbl[[v]]; inv <- invalid_by_var[[v]]; !(vals %in% inv)
    })
    cert_tbl$complete_all <- rowSums(valid_mat) == length(demo_vars)
    
    county_metrics <- cert_tbl |>
        group_by(.data[[county_var]], year) |>
        summarise(
            n_cert = n(),
            garb_k = sum(is_garbage),
            prop_garbage = garb_k / n_cert,
            prop_garbage_low = wilson_lower(garb_k, n_cert),
            prop_garbage_hi  = wilson_upper(garb_k, n_cert),
            
            all_comp_k = sum(complete_all_exc_race),
            prop_all_comp     = all_comp_k / n_cert,
            prop_all_comp_low = wilson_lower(all_comp_k, n_cert),
            prop_all_comp_hi  = wilson_upper(all_comp_k, n_cert),
            
            contrib_n  = sum(total_contrib),
            light_k    = sum(light_contrib),
            prop_light = ifelse(contrib_n > 0, light_k / contrib_n, NA_real_),
            prop_light_low = wilson_lower(light_k, contrib_n),
            prop_light_hi  = wilson_upper(light_k, contrib_n),
            
            acc_n       = sum(needs_detail),
            acc_miss_k  = sum(needs_detail & !has_required_detail),
            pct_acc_miss     = ifelse(acc_n > 0, acc_miss_k / acc_n, NA_real_),
            pct_acc_miss_low = wilson_lower(acc_miss_k, acc_n),
            pct_acc_miss_hi  = wilson_upper(acc_miss_k, acc_n),
            
            overd_n      = sum(is_overdose),
            overd_miss_k = sum(is_overdose & !has_required_detail),
            pct_overd_miss     = ifelse(overd_n > 0, overd_miss_k / overd_n, NA_real_),
            pct_overd_miss_low = wilson_lower(overd_miss_k, overd_n),
            pct_overd_miss_hi  = wilson_upper(overd_miss_k, overd_n),
            
            overdose_unspec_k = sum(is_overdose & unspecific_drug),
            
            across(all_of(demo_vars),
                   ~ { bad <- invalid_by_var[[cur_column()]]
                   good <- !(.x %in% bad | is.na(.x))
                   if (all(is.na(.x))) NA_integer_ else sum(good) },
                   .names = "{.col}_comp_k"
            ),
            .groups = "drop"
        )
    
    # === SIM + Foreman DQ entropy (join by county + year) ===
    entropy_tbl <- compute_entropy_county_foreman(
        ds,
        county_var,
        dict_dir      = dictionary_dir,
        icd_map_path  = file.path(dictionary_dir, "foreman-icd10-mapping.csv"),
        code_map_path = file.path(dictionary_dir, "foreman-table2-map.csv")
    )
    
    result <- county_metrics |>
        left_join(entropy_tbl, by = c(county_var, "year")) |>
        mutate(DQ_prop_garbage = (1 - DQ_overall) * prop_garbage)
    
    
    result
}

# ————————————————————————————————————————————————————————————————
#  run for all years
# ————————————————————————————————————————————————————————————————
parquet_files <- list.files(parquet_dir, "\\.parquet$", full.names = TRUE)
if (!is.null(years_wanted)) {
    parquet_files <- parquet_files[stringr::str_extract(basename(parquet_files), "\\d{4}") %in% years_wanted]
}

if (PARALLELIZE) {
    suppressPackageStartupMessages({ library(future); library(furrr) })
    plan(multisession, workers = max(1, N_WORKERS))
    county_year_all <- furrr::future_map_dfr(
        parquet_files,
        summarise_year_file,
        .options = furrr::furrr_options(
            seed = TRUE,
            packages = c("arrow","tidyr","Matrix","glmnet","stringr","readr","dplyr")
        )
    )
    plan(sequential)
} else {
    county_year_all <- purrr::map_dfr(parquet_files, summarise_year_file)
}

if (WRITE_ONE_CSV) {
    readr::write_csv(county_year_all, out_csv)
    message("County-year quality metrics written: ", out_csv)
} else {
    message("Rows written to dataset dir: ", WRITE_DATASET_DIR)
}

BASE_DIR = Path("/Users/amymann/Documents/Data Quality Project/data/parquet")


# read_parquet(here("data_private/mcod_sample/mcod_1999.parquet"))
read_parquet("/Users/amymann/Documents/Data Quality Project/data/parquet/mort2020.parquet")
