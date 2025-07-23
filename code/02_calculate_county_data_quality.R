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
# - gbd_garbage_codes_without_overdose.csv containing the IHME list of garbage codes
# - dq_entropy_helper.R which returns a table of DQ values
#
# All these files can be found on Github: https://github.com/am-mann/data-quality/
#
# Output: county_year_quality_metrics.csv
#
# ————————————————————————————————————————————————————————————————
library(arrow)
library(dplyr)
library(purrr)
library(stringr)
library(binom)
library(readr)
library(tidyr)
library(future)
library(furrr)
library(here)

## CONSTANTS ----
PARALLELIZE <- TRUE

# ————————————————————————————————————————————————————————————————
#  load things + dictionaries + helpers
# ————————————————————————————————————————————————————————————————
source(here("code/helpers", "dq_entropy_helper.R"))   # loads entropy helper functions
parquet_dir    <- here("data_private", "mcod_sample")
# parquet_dir    <- here("data_private", "mcod_sample") 
dictionary_dir <- here("data_raw", "cause-codes")
out_csv        <- here("data", "county_year_quality_metrics.csv")

county_var     <- "countyrs" 
years_wanted   <- NULL
key_demo_vars  <- c("marstat", "placdth", "educ", "mandeath", "age", "sex", "race", "racer40")

canonical_icd <- function(x) {
    stringr::str_remove_all(stringr::str_to_upper(x), "[^A-Z0-9]")
}

lookup_garbage <- read_csv(file.path(dictionary_dir,
    "gbd_garbage_codes_without_overdose.csv"),
show_col_types = FALSE) %>%
    transmute(icd10 = canonical_icd(icd10),
        gbd_severity = as.integer(gbd_severity)) %>%
    distinct(icd10, .keep_all = TRUE)

light_garbage <- read_csv(
    file.path(dictionary_dir, "light_garbage_codes.csv"),
    show_col_types = FALSE
)$code %>%
    canonical_icd()

overdose_path <- file.path(dictionary_dir, "overdose-detail-codes-v3.csv")
overdose_ucod <- read_csv(overdose_path, show_col_types = FALSE)$ucod %>% canonical_icd()

# Looks up accident details
accident_detail_lookup <- readr::read_csv(
    file.path(dictionary_dir, "accident-detail-codes-v3.csv"),
    show_col_types = FALSE
) |>
    tidyr::pivot_longer(
        starts_with("detail_"),
        values_to      = "detail",
        values_drop_na = TRUE
    ) |>
    transmute(
        ucod   = canonical_icd(stringr::str_trim(ucod)),
        detail = canonical_icd(stringr::str_trim(detail))
    ) |>
    filter(detail != "") |>
    distinct()

accident_ucod_set   <- unique(accident_detail_lookup$ucod)
accident_detail_set <- unique(accident_detail_lookup$detail)

# Helpers
clean_icd <- function(x) str_remove_all(str_to_upper(x), "[^A-Z0-9\\.]")

num <- function(x) if (bit64::is.integer64(x)) as.numeric(x) else as.numeric(x)

wilson_lower <- function(k, n) {
    k <- num(k)
    n <- num(n)
    out <- rep(NA_real_, length(k))
    ok  <- n > 0 & !is.na(n)
    if (any(ok)) out[ok] <- binom.confint(k[ok], n[ok], methods = "wilson")$lower
    out
}
wilson_upper <- function(k, n) {
    k <- num(k)
    n <- num(n)
    out <- rep(NA_real_, length(k))
    ok  <- n > 0 & !is.na(n)
    if (any(ok)) out[ok] <- binom.confint(k[ok], n[ok], methods = "wilson")$upper
    out
}

sec_cols <- paste0("record_", 1:20)
invalid_by_var <- list(
    marstat   = c("", "not_stated", "unknown", NA),
    placdth   = c("", "9", "U", "Unknown", NA),
    educ2003  = c("", "9", "99", "999", "9999", "0", NA),
    educ  = c("", "9", "99", "999", "9999", "0", NA),
    mandeath  = c("", "9", "7", "Unknown", NA),
    age       = c("", "999", "9999", "0000", NA),
    sex       = c("", "U", NA),
    race      = c("", "9", "Unknown", NA),
    racer40      = c("", "9", "Unknown", NA)
)

# ─────────────────────────────────────────────────────────────────────────────
#  summarise everything for one year
# ─────────────────────────────────────────────────────────────────────────────
summarise_year_file <- function(file) {
    this_year <- as.integer(stringr::str_extract(basename(file), "\\d{4}"))
    message("Processing ", this_year)

    cols_needed <- c("ucod", county_var, key_demo_vars, sec_cols, "ager27", "sex", "age_years")
    ds <- arrow::read_parquet(
        file,
        as_data_frame = TRUE,
        col_select = intersect(cols_needed,
            names(arrow::read_parquet(file, as_data_frame = FALSE)$schema)))
    
    ds <- ds %>%
        {                        # rename record1→record_1 etc. if underscore cols are absent
            if (!any(grepl("^record_", names(.))) &&
                any(grepl("^record[0-9]+$", names(.))))
                rename_with(., ~ sub("^record([0-9]+)$", "record_\\1", .x),
                            .cols = matches("^record[0-9]+$"))
            else .
        } %>%
        mutate(
            ranum    = row_number(),
            uc4      = substr(clean_icd(ucod), 1, 4),
            uc3      = substr(clean_icd(ucod), 1, 3),
            across(starts_with("record_"), ~ canonical_icd(.x)),
            sex_male = if_else(sex == "M", 1L, 0L))
    
    entropy_tbl <- compute_entropy_county(ds, county_var)
    
    # from 2003 to 2005 the educ2003 variable is called educ so here I change the name accordingly
    # if (this_year %in% 2003:2005 &&
    #     "educ" %in% names(ds)      &&
    #     !"educ2003" %in% names(ds)) {
    #     ds <- dplyr::rename(ds, educ2003 = educ)
    # }

    for (v in setdiff(key_demo_vars, names(ds))) ds[[v]] <- NA_character_
    if (!"year" %in% names(ds)) ds$year <- this_year
    if (this_year < 2004) ds$educ2003 <- NA_character_
    demo_vars <- if (this_year < 2004) {
        setdiff(key_demo_vars, "educ2003")
    } else {key_demo_vars}
    if (this_year %in% 2021:2022) {
        demo_vars <- union(setdiff(demo_vars, "race"), "racer40")
    }

    # flag accidents and overdoses
    ds <- ds %>% mutate(
        icd10 = canonical_icd(ucod),
        ranum = dplyr::row_number()
    ) %>%
        left_join(lookup_garbage, by = "icd10") %>%
        mutate(
            needs_detail = icd10 %in% accident_ucod_set,
            is_overdose  = icd10 %in% overdose_ucod
        )

    # gets all contributing causes
    present_sec <- intersect(sec_cols, names(ds))

    contrib_all <- if (length(present_sec) == 0) {
        tibble(ranum = integer(), detail_full = character(), detail = character())
    } else {
        present_sec %>%
            purrr::map_dfr(~ ds %>% select(ranum, detail_raw = !!sym(.x))) %>%
            mutate(
                detail_full = canonical_icd(detail_raw),        # full code
                detail      = substr(detail_full, 1, 4)        # 4-char version
            ) %>%
            filter(!is.na(detail_full) & trimws(detail_full) != "")}
    

    # keep relevant details for accident deaths
    contrib_long <- contrib_all %>%
        dplyr::filter(detail %in% accident_detail_set)

    # sum light garbage
    contrib_counts <- contrib_all %>%
        mutate(is_light = detail_full %in% light_garbage) %>%  
        group_by(ranum) %>%
        summarise(
            total_contrib = n(),
            light_contrib = sum(is_light),
            .groups = "drop"
        )

    ## flag unspecified drug for overdoses
    unspec_tbl <- contrib_all %>%
        filter(str_starts(detail_full, "T")) %>%      
        mutate(
            is_specific = !(str_starts(detail_full, "T509") |
                                str_starts(detail_full, "T409"))
        ) %>%
        group_by(ranum) %>%
        summarise(has_specific = any(is_specific), .groups = "drop") %>%
        mutate(unspecific_drug = !has_specific) %>%
        select(ranum, unspecific_drug)
    

    # merge all flags
    cert_tbl <- ds %>%
        dplyr::left_join(contrib_counts,  by = "ranum") %>%
        dplyr::left_join(unspec_tbl,      by = "ranum") %>%
        dplyr::mutate(
            total_contrib   = tidyr::replace_na(total_contrib, 0L),
            light_contrib   = tidyr::replace_na(light_contrib, 0L),
            unspecific_drug = tidyr::replace_na(unspecific_drug, TRUE),
            is_garbage      = dplyr::case_when(
                is.na(gbd_severity)                          ~ FALSE,
                !stringr::str_starts(icd10, "X")             ~ TRUE,
                stringr::str_starts(icd10, "X") & light_contrib == 0 ~ TRUE,
                TRUE                                         ~ FALSE
            )
        )

    # certificates that need detail
    need_tbl <- cert_tbl %>%
        dplyr::filter(needs_detail) %>%
        dplyr::distinct(ranum, icd10)
    
    demo_vars_exc_race <- setdiff(demo_vars, c("race", "racer40", "mandeath"))
    valid_exc_race <- sapply(demo_vars_exc_race, function(v) {
        vals <- cert_tbl[[v]]
        inv  <- invalid_by_var[[v]]
        !(vals %in% inv)
    })
    cert_tbl$complete_all_exc_race <- rowSums(valid_exc_race) == length(demo_vars_exc_race)
    
    # find required detail
    required_hits <- contrib_long %>%
        dplyr::semi_join(accident_detail_lookup, by = "detail") %>%
        dplyr::left_join(ds %>% dplyr::select(ranum, icd10), by = "ranum") %>%
        dplyr::semi_join(accident_detail_lookup,
            by = c("icd10" = "ucod", "detail")) %>%
        dplyr::distinct(ranum) %>%
        dplyr::mutate(has_required_detail = TRUE)

    cert_tbl <- cert_tbl %>%
        dplyr::left_join(required_hits, by = "ranum") %>%
        dplyr::mutate(has_required_detail = tidyr::replace_na(has_required_detail, FALSE))

    # check completeness of demographic information
    valid_mat <- sapply(demo_vars, function(v) {
        vals <- cert_tbl[[v]]
        inv  <- invalid_by_var[[v]]
        !(vals %in% inv)
        
    })
    
    cert_tbl$complete_all <- rowSums(valid_mat) == length(demo_vars)
 
    # make big table
    county_metrics <- cert_tbl %>%
        dplyr::group_by(.data[[county_var]], year) %>%
        dplyr::summarise(
            n_cert               = n(),
            garb_k               = sum(is_garbage),
            prop_garbage         = garb_k / n_cert,
            prop_garbage_low     = wilson_lower(garb_k, n_cert),
            prop_garbage_hi      = wilson_upper(garb_k, n_cert),
            pct_gc_I64        = mean(icd10 == "I64"),
            pct_gc_C_misc     = mean(icd10 %in% c("C80", "C55", "C97")),
            pct_gc_I10        = mean(icd10 == "I10"),
            pct_gc_R_misc     = mean(icd10 %in% c("R99", "R54")),
            pct_gc_N19        = mean(icd10 == "N19"),
            pct_gc_J80        = mean(icd10 == "J80"),
            all_comp_k           = sum(complete_all_exc_race),
            prop_all_comp        = all_comp_k / n_cert,
            prop_all_comp_low    = wilson_lower(all_comp_k, n_cert),
            prop_all_comp_hi     = wilson_upper(all_comp_k, n_cert),
            contrib_n            = sum(total_contrib),
            light_k              = sum(light_contrib),
            prop_light           = ifelse(contrib_n > 0, light_k / contrib_n, NA_real_),
            prop_light_low       = wilson_lower(light_k, contrib_n),
            prop_light_hi        = wilson_upper(light_k, contrib_n),
            acc_n                = sum(needs_detail),
            acc_miss_k           = sum(needs_detail & !has_required_detail),
            pct_acc_miss         = ifelse(acc_n > 0, acc_miss_k / acc_n, NA_real_),
            pct_acc_miss_low     = wilson_lower(acc_miss_k, acc_n),
            pct_acc_miss_hi      = wilson_upper(acc_miss_k, acc_n),
            overd_n              = sum(is_overdose),
            overd_miss_k         = sum(is_overdose & !has_required_detail),
            pct_overd_miss       = ifelse(overd_n > 0, overd_miss_k / overd_n, NA_real_),
            pct_overd_miss_low   = wilson_lower(overd_miss_k, overd_n),
            pct_overd_miss_hi    = wilson_upper(overd_miss_k, overd_n),      overdose_unspec_k    = sum(is_overdose & unspecific_drug),
            across(all_of(demo_vars),
                   ~ {bad_vals <- invalid_by_var[[cur_column()]]
                    good     <- !(.x %in% bad_vals | is.na(.x))
                    if (all(is.na(.x))) NA_integer_ else sum(good) },
                .names = "{.col}_comp_k"
            ),
            .groups = "drop"
        )
    result <- dplyr::left_join(county_metrics, entropy_tbl, by = county_var) %>%
        dplyr::mutate(
            DQ_prop_garbage = (1-DQ_overall) * prop_garbage)
}

# ————————————————————————————————————————————————————————————————
#  run for all years + output
# ————————————————————————————————————————————————————————————————
parquet_files <- list.files(parquet_dir, "\\.parquet$", full.names = TRUE)
if (!is.null(years_wanted)) {
    parquet_files <- parquet_files[
        str_extract(basename(parquet_files), "\\d{4}") %in% years_wanted
    ]
}

if (PARALLELIZE) {
    future::plan(future::multisession(workers = future::availableCores() - 2))
    county_year_all <- furrr::future_map_dfr(
        parquet_files,
        summarise_year_file,
        .options = furrr_options(seed = TRUE, packages = c("arrow", "tidyr", "forcats", "Matrix", "nnet", "glmnet")))
} else {
    county_year_all <- map_dfr(parquet_files, summarise_year_file)
}

write_csv(county_year_all, out_csv)
message("County-year quality metrics written")

# file <- "data_private/mcod_sample/mcod_2006.parquet"
# df <- read_parquet(file)
