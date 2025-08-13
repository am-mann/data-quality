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
library(fs)
library(iNEXT)

## CONSTANTS ----
PARALLELIZE <- TRUE

# ————————————————————————————————————————————————————————————————
#  load things + dictionaries + helpers
# ————————————————————————————————————————————————————————————————
source(here("code/helpers", "dq_entropy_helper.R")) # loads entropy helper functions
if (dir_exists(here("data_private", "mcod"))) {
    parquet_dir <- here("data_private", "mcod")
} else {
    parquet_dir <- here("data_private", "mcod_sample")
}
dictionary_dir <- here("data_raw", "cause-codes")
out_csv <- here("data", "county_year_quality_metrics.csv.gz")

county_var <- "county_ihme"
years_wanted <- NULL
key_demo_vars <- c("marstat", "placdth", "educ", "mandeath", "age", "sex", "race", "racer40")

canonical_icd <- function(x) {
    stringr::str_remove_all(stringr::str_to_upper(x), "[^A-Z0-9]")
}

lookup_garbage <- read_csv(
    file.path(
        dictionary_dir,
        "gbd_garbage_codes_without_overdose.csv"
    ),
    show_col_types = FALSE
) %>%
    transmute(
        icd10 = canonical_icd(icd10),
        gbd_severity = as.integer(gbd_severity)
    ) %>%
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

accident_ucod_set <- unique(accident_detail_lookup$ucod)
accident_detail_set <- unique(accident_detail_lookup$detail)

# Helpers
clean_icd <- function(x) str_remove_all(str_to_upper(x), "[^A-Z0-9\\.]")

num <- function(x) if (bit64::is.integer64(x)) as.numeric(x) else as.numeric(x)

wilson_lower <- function(k, n) {
    k <- num(k)
    n <- num(n)
    out <- rep(NA_real_, length(k))
    ok <- n > 0 & !is.na(n)
    if (any(ok)) out[ok] <- binom.confint(k[ok], n[ok], methods = "wilson")$lower
    out
}
wilson_upper <- function(k, n) {
    k <- num(k)
    n <- num(n)
    out <- rep(NA_real_, length(k))
    ok <- n > 0 & !is.na(n)
    if (any(ok)) out[ok] <- binom.confint(k[ok], n[ok], methods = "wilson")$upper
    out
}

sec_cols <- paste0("record_", 1:20)
invalid_by_var <- list(
    marstat = c("", "not_stated", "unknown", NA),
    placdth = c("", "9", "U", "Unknown", NA),
    educ2003 = c("", "9", "99", "999", "9999", "0", NA),
    educ = c("", "9", "99", "999", "9999", "0", NA),
    mandeath = c("", "9", "7", "Unknown", NA),
    age = c("", "999", "9999", "0000", NA),
    sex = c("", "U", NA),
    race = c("", "9", "Unknown", NA),
    racer40 = c("", "9", "Unknown", NA)
)

# ─────────────────────────────────────────────────────────────────────────────
#  summarise everything for one year
# ─────────────────────────────────────────────────────────────────────────────
summarise_year_file <- function(file) {
    this_year <- as.integer(stringr::str_extract(basename(file), "\\d{4}"))
    
    cat(sprintf("Processing %d ...\n", this_year))
    flush.console()
    
    cols_needed <- c("ucod", county_var, key_demo_vars, sec_cols, "ager27", "sex", "age_years", "ucr358")
    ds <- arrow::read_parquet(
        file,
        as_data_frame = TRUE,
        col_select = intersect(
            cols_needed,
            names(arrow::read_parquet(file, as_data_frame = FALSE)$schema)
        )
    )
    
    ds <- ds %>%
        { # rename record1→record_1 etc. if underscore cols are absent
            if (!any(grepl("^record_", names(.))) &&
                any(grepl("^record[0-9]+$", names(.)))) {
                rename_with(., ~ sub("^record([0-9]+)$", "record_\\1", .x),
                            .cols = matches("^record[0-9]+$")
                )
            } else {
                .
            }
        } %>%
        mutate(
            ranum    = row_number(),
            uc4      = substr(clean_icd(ucod), 1, 4),
            uc3      = substr(clean_icd(ucod), 1, 3),
            across(starts_with("record_"), ~ canonical_icd(.x)),
            sex_male = if_else(sex == "M", 1L, 0L)
        )
    
    entropy_tbl <- compute_entropy_county(ds, county_var)
    
    for (v in setdiff(key_demo_vars, names(ds))) ds[[v]] <- NA_character_
    if (!"year" %in% names(ds)) ds$year <- this_year
    if (this_year < 2004) ds$educ2003 <- NA_character_
    demo_vars <- if (this_year < 2004) setdiff(key_demo_vars, "educ2003") else key_demo_vars
    if (this_year %in% 2021:2022) demo_vars <- union(setdiff(demo_vars, "race"), "racer40")
    
    # ---- Garbage/overdose/accident prep (unchanged) ----
    ds <- ds %>%
        mutate(
            icd10 = canonical_icd(ucod),
            ranum = dplyr::row_number()
        ) %>%
        left_join(lookup_garbage, by = "icd10") %>%
        mutate(
            needs_detail = icd10 %in% accident_ucod_set,
            is_overdose  = icd10 %in% overdose_ucod
        )
    
    present_sec <- intersect(sec_cols, names(ds))
    contrib_all <- if (length(present_sec) == 0) {
        tibble(ranum = integer(), detail_full = character(), detail = character())
    } else {
        present_sec %>%
            purrr::map_dfr(~ ds %>% select(ranum, detail_raw = !!sym(.x))) %>%
            mutate(
                detail_full = canonical_icd(detail_raw),
                detail      = substr(detail_full, 1, 4)
            ) %>%
            filter(!is.na(detail_full) & trimws(detail_full) != "")
    }
    
    # keep only required accident details for diagnostics later
    contrib_long <- contrib_all %>% filter(detail %in% accident_detail_set)
    
    # light garbage counts in MCOD
    contrib_counts <- contrib_all %>%
        mutate(is_light = detail_full %in% light_garbage) %>%
        group_by(ranum) %>%
        summarise(
            total_contrib = n(),
            light_contrib = sum(is_light),
            .groups = "drop"
        )
    
    # unspecific drug flags (optional diagnostic)
    unspec_tbl <- contrib_all %>%
        filter(str_starts(detail_full, "T")) %>%
        mutate(
            is_specific = !(str_starts(detail_full, "T509") | str_starts(detail_full, "T409"))
        ) %>%
        group_by(ranum) %>%
        summarise(has_specific = any(is_specific), .groups = "drop") %>%
        mutate(unspecific_drug = !has_specific) %>%
        select(ranum, unspecific_drug)
    
    # ======================================================================
    # SIM (Self-Injury Mortality) — Rockett et al.
    # SIM = suicides (all mechanisms) + 0.8 * unintentional drug poisonings (X40–X44)
    #                      + 0.9 * undetermined drug poisonings (Y10–Y14)
    # Restricted to ages ≥ 15.
    # ======================================================================
    
    # Age ≥ 15 flag using age_years if available, otherwise ager27 (>= 12 ~ 15–19+)
    adult_15p <- if ("age_years" %in% names(ds)) {
        !is.na(ds$age_years) & as.numeric(ds$age_years) >= 15
    } else if ("ager27" %in% names(ds)) {
        # NCHS ager27: 12 = 15–19, ..., 26 = 85+; treat >=12 as 15+
        !is.na(ds$ager27) & as.integer(ds$ager27) >= 12
    } else {
        # Fallback: assume all ages (conservative—consider adding a warning/log)
        rep(TRUE, nrow(ds))
    }
    
    # ICD-10 helpers
    is_suicide_ucod <- function(uc4) {
        r3 <- toupper(substr(uc4, 1, 3))
        xnum <- suppressWarnings(as.integer(sub("^X", "", r3)))
        in_x60_84 <- substr(r3, 1, 1) == "X" & !is.na(xnum) & xnum >= 60 & xnum <= 84
        y870 <- toupper(substr(uc4, 1, 4)) == "Y870"
        in_x60_84 | y870
    }
    is_unintent_drug_ucod <- function(uc4) {
        r3 <- toupper(substr(uc4, 1, 3))
        xnum <- suppressWarnings(as.integer(sub("^X", "", r3)))
        substr(r3, 1, 1) == "X" & !is.na(xnum) & xnum >= 40 & xnum <= 44
    }
    is_undet_drug_ucod <- function(uc4) {
        r3 <- toupper(substr(uc4, 1, 3))
        ynum <- suppressWarnings(as.integer(sub("^Y", "", r3)))
        substr(r3, 1, 1) == "Y" & !is.na(ynum) & ynum >= 10 & ynum <= 14
    }
    
    # ---- Detect fentanyl in MCOD (T40.4) ----
    fentanyl_ranum <- contrib_all %>%
        filter(toupper(substr(detail_full, 1, 4)) == "T404") %>%
        distinct(ranum) %>%
        pull(ranum)
    
    # ======================================================================
    # SIM (Self-Injury Mortality) — Rockett et al. with fentanyl excluded
    # ======================================================================
    
    sim_flags <- ds %>%
        transmute(
            !!county_var := .data[[county_var]],
            year,
            is_adult15p = adult_15p,
            ranum,
            uc4 = toupper(substr(ucod, 1, 4)),
            suic  = is_suicide_ucod(uc4),
            uninj = is_unintent_drug_ucod(uc4),
            undet = is_undet_drug_ucod(uc4)
        ) %>%
        filter(is_adult15p & !(ranum %in% fentanyl_ranum))
    
    sim_county <- sim_flags %>%
        reframe(
            sim_suic_k_15p_nofent   = sum(suic,  na.rm = TRUE),
            sim_uninj_k_15p_nofent  = sum(uninj, na.rm = TRUE),   # X40–X44
            sim_undet_k_15p_nofent  = sum(undet, na.rm = TRUE),   # Y10–Y14
            sim_added_k_15p_nofent  = 0.8 * sim_uninj_k_15p_nofent +
                0.9 * sim_undet_k_15p_nofent,
            sim_k_15p_nofent        = sim_suic_k_15p_nofent + sim_added_k_15p_nofent,
            sim_added_share_nofent  = ifelse(sim_k_15p_nofent > 0,
                                             sim_added_k_15p_nofent / sim_k_15p_nofent,
                                             NA_real_),
            .by = c(!!sym(county_var), year)
        )
    
    
    # ---- Merge diagnostic flags and completeness (unchanged) ----
    cert_tbl <- ds %>%
        left_join(contrib_counts, by = "ranum") %>%
        left_join(unspec_tbl,    by = "ranum") %>%
        mutate(
            total_contrib   = tidyr::replace_na(total_contrib, 0L),
            light_contrib   = tidyr::replace_na(light_contrib, 0L),
            unspecific_drug = tidyr::replace_na(unspecific_drug, TRUE),
            is_garbage = case_when(
                is.na(gbd_severity) ~ FALSE,
                !stringr::str_starts(icd10, "X") ~ TRUE,
                stringr::str_starts(icd10, "X") & unspecific_drug == TRUE ~ TRUE,
                TRUE ~ FALSE
            )
        )
    
    need_tbl <- cert_tbl %>% filter(needs_detail) %>% distinct(ranum, icd10)
    
    demo_vars_exc_race <- setdiff(demo_vars, c("race", "racer40", "mandeath"))
    valid_exc_race <- sapply(demo_vars_exc_race, function(v) {
        vals <- cert_tbl[[v]]
        inv  <- invalid_by_var[[v]]
        !(vals %in% inv)
    })
    cert_tbl$complete_all_exc_race <- rowSums(valid_exc_race) == length(demo_vars_exc_race)
    
    required_hits <- contrib_long %>%
        semi_join(accident_detail_lookup, by = "detail") %>%
        left_join(ds %>% select(ranum, icd10), by = "ranum") %>%
        semi_join(accident_detail_lookup, by = c("icd10" = "ucod", "detail")) %>%
        distinct(ranum) %>%
        mutate(has_required_detail = TRUE)
    
    cert_tbl <- cert_tbl %>%
        left_join(required_hits, by = "ranum") %>%
        mutate(has_required_detail = tidyr::replace_na(has_required_detail, FALSE))
    
    valid_mat <- sapply(demo_vars, function(v) {
        vals <- cert_tbl[[v]]
        inv  <- invalid_by_var[[v]]
        !(vals %in% inv)
    })
    cert_tbl$complete_all <- rowSums(valid_mat) == length(demo_vars)
    
    county_metrics <- cert_tbl %>%
        group_by(.data[[county_var]], year) %>%
        summarise(
            n_cert = n(),
            garb_k = sum(is_garbage),
            prop_garbage = garb_k / n_cert,
            prop_garbage_low = wilson_lower(garb_k, n_cert),
            prop_garbage_hi  = wilson_upper(garb_k, n_cert),
            pct_gc_I64   = mean(icd10 == "I64"),
            pct_gc_C_misc= mean(icd10 %in% c("C80", "C55", "C97")),
            pct_gc_I10   = mean(icd10 == "I10"),
            pct_gc_R_misc= mean(icd10 %in% c("R99", "R54")),
            pct_gc_N19   = mean(icd10 == "N19"),
            pct_gc_J80   = mean(icd10 == "J80"),
            all_comp_k   = sum(complete_all_exc_race),
            prop_all_comp     = all_comp_k / n_cert,
            prop_all_comp_low = wilson_lower(all_comp_k, n_cert),
            prop_all_comp_hi  = wilson_upper(all_comp_k, n_cert),
            contrib_n   = sum(total_contrib),
            light_k     = sum(light_contrib),
            prop_light      = ifelse(contrib_n > 0, light_k / contrib_n, NA_real_),
            prop_light_low  = wilson_lower(light_k, contrib_n),
            prop_light_hi   = wilson_upper(light_k, contrib_n),
            acc_n       = sum(needs_detail),
            acc_miss_k  = sum(needs_detail & !has_required_detail),
            pct_acc_miss     = ifelse(acc_n > 0, acc_miss_k / acc_n, NA_real_),
            pct_acc_miss_low = wilson_lower(acc_miss_k, acc_n),
            pct_acc_miss_hi  = wilson_upper(acc_miss_k, acc_n),
            overd_n     = sum(is_overdose),
            overd_miss_k= sum(is_overdose & !has_required_detail),
            pct_overd_miss     = ifelse(overd_n > 0, overd_miss_k / overd_n, NA_real_),
            pct_overd_miss_low = wilson_lower(overd_miss_k, overd_n),
            pct_overd_miss_hi  = wilson_upper(overd_miss_k, overd_n),
            overdose_unspec_k  = sum(is_overdose & unspecific_drug),
            across(all_of(demo_vars),
                   ~ {
                       bad_vals <- invalid_by_var[[cur_column()]]
                       good <- !(.x %in% bad_vals | is.na(.x))
                       if (all(is.na(.x))) NA_integer_ else sum(good)
                   },
                   .names = "{.col}_comp_k"
            ),
            .groups = "drop"
        )
    
    county_metrics <- county_metrics %>%
        left_join(sim_county, by = c(county_var, "year"))
    
    result <- left_join(county_metrics, entropy_tbl, by = county_var) %>%
        mutate(DQ_prop_garbage = (1 - DQ_overall) * prop_garbage)
    
    return(result)
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
    future::plan(future::multisession(workers = future::availableCores() - 6))
    county_year_all <- furrr::future_map_dfr(
        parquet_files,
        summarise_year_file,
        .options = furrr_options(seed = TRUE, packages = c("arrow", "tidyr", "forcats", "Matrix", "nnet", "glmnet"))
    )
    future::plan(future::sequential())
} else {
    county_year_all <- map_dfr(parquet_files, summarise_year_file)
}

write_csv(county_year_all, out_csv)
message("County-year quality metrics written")


