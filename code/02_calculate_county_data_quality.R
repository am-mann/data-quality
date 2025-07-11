# ————————————————————————————————————————————————————————————————
# Date: July 7, 2025
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
#
# All these files can be found on Github: https://github.com/am-mann/data-quality/tree/main/cause-codes
#
# Output: county_year_quality_metrics.csv
#
# ***indicates that a line needs to be changed prior to running
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

## CONSTANTS ----
PARALLELIZE <- TRUE

# ————————————————————————————————————————————————————————————————
#  load things + dictionaries + helpers
# ————————————————————————————————————————————————————————————————
parquet_dir    <- here("data_private", "mcod")  #### ***Change this line
dictionary_dir <- here("data_raw", "cause-codes") #### ***Change this line
out_csv        <- here("data", "county_year_quality_metrics.csv")

county_var     <- "countyrs"  #### ***Put name of county variable
years_wanted   <- NULL
key_demo_vars  <- c("marstat", "placdth", "educ", "mandeath", "age", "sex", "race")

lookup_garbage <- read_csv(file.path(dictionary_dir,
    "gbd_garbage_codes_without_overdose.csv"),
show_col_types = FALSE) %>%
    transmute(icd10 = str_to_upper(icd10),
        gbd_severity = as.integer(gbd_severity)) %>%
    distinct(icd10, .keep_all = TRUE)

light_garbage <- read_csv(file.path(dictionary_dir, "light_garbage_codes.csv"),
    show_col_types = FALSE)$icd10 %>%
    str_to_upper()

overdose_path <- file.path(dictionary_dir, "overdose-detail-codes-v3.csv")
overdose_ucod <- read_csv(overdose_path, show_col_types = FALSE)$ucod %>% str_to_upper()

# Looks up accident details
canonical_icd <- function(x) {
    stringr::str_remove_all(stringr::str_to_upper(x), "[^A-Z0-9]")
}

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
invalid  <- c("", "9", "U", "Unknown")

# ─────────────────────────────────────────────────────────────────────────────
#  summarise everything for one year
# ─────────────────────────────────────────────────────────────────────────────
summarise_year_file <- function(file) {
    this_year <- as.integer(stringr::str_extract(basename(file), "\\d{4}"))
    message("Processing ", this_year)

    cols_needed <- c("ranum", "ucod", county_var, key_demo_vars, sec_cols)
    ds <- arrow::read_parquet(
        file,
        as_data_frame = TRUE,
        col_select = intersect(cols_needed,
            names(arrow::read_parquet(file, as_data_frame = FALSE)$schema))
    )

    # from 2003 to 2005 the educ2003 variable is called educ so here I change the name accordingly
    # if (this_year %in% 2003:2005 &&
    #     "educ" %in% names(ds)      &&
    #     !"educ2003" %in% names(ds)) {
    #     ds <- dplyr::rename(ds, educ2003 = educ)
    # }

    for (v in setdiff(key_demo_vars, names(ds))) ds[[v]] <- NA_character_
    # if (!county_var %in% names(ds)) {
    #     set.seed(this_year)
    #     ### *** This line assigns a randomly generated county and should be deleted and replaced with your county variable #####
    #     ds[[county_var]] <- sprintf("%05d", sample(10001:56045, nrow(ds), TRUE))
    # }
    if (!"year" %in% names(ds)) ds$year <- this_year
    if (this_year < 2004) ds$educ2003 <- NA_character_
    demo_vars <- if (this_year < 2004) setdiff(key_demo_vars, "educ2003") else key_demo_vars

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
        tibble(ranum = integer(), detail = character())
    } else {
        present_sec %>%
            purrr::map(~ ds %>% dplyr::select(ranum, detail_raw = !!sym(.x))) %>%
            dplyr::bind_rows() %>%
            dplyr::mutate(detail = canonical_icd(detail_raw)) %>%
            dplyr::filter(!is.na(detail_raw) & trimws(detail_raw) != "")
    }

    # keep relevant details for accident deaths
    contrib_long <- contrib_all %>%
        dplyr::filter(detail %in% accident_detail_set)

    # sum light garbage
    contrib_counts <- contrib_all %>%
        dplyr::mutate(is_light = detail %in% light_garbage) %>%
        dplyr::group_by(ranum) %>%
        dplyr::summarise(
            total_contrib = n(),
            light_contrib = sum(is_light),
            .groups = "drop"
        )

    ## flag unspecified drug for overdoses
    unspec_tbl <- contrib_all %>%
        dplyr::filter(stringr::str_starts(detail, "T")) %>%
        dplyr::mutate(
            is_specific = !(stringr::str_starts(detail, "T509") |
                stringr::str_starts(detail, "T409"))
        ) %>%
        dplyr::group_by(ranum) %>%
        dplyr::summarise(has_specific = any(is_specific), .groups = "drop") %>%
        dplyr::mutate(unspecific_drug = !has_specific) %>%
        dplyr::select(ranum, unspecific_drug)

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
    valid_mat <- !as.matrix(cert_tbl[demo_vars]) %in% invalid &
        !is.na(as.matrix(cert_tbl[demo_vars]))
    cert_tbl$complete_all <- rowSums(valid_mat) == length(demo_vars)

    # make big table
    cert_tbl %>%
        dplyr::group_by(.data[[county_var]], year) %>%
        dplyr::summarise(
            n_cert               = n(),
            garb_k               = sum(is_garbage),
            prop_garbage         = garb_k / n_cert,
            prop_garbage_low     = wilson_lower(garb_k, n_cert),
            prop_garbage_hi      = wilson_upper(garb_k, n_cert),
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
            across(all_of(key_demo_vars),
                ~ {
                    valid <- !(.x %in% invalid | is.na(.x))
                    if (all(is.na(.x))) NA_integer_ else sum(valid)
                },
                .names = "{.col}_comp_k"),
            .groups = "drop"
        )
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
        .options = furrr_options(packages = c("forcats", "arrow")))
} else {
    county_year_all <- map_dfr(parquet_files, summarise_year_file)
}

write_csv(county_year_all, out_csv)
message("County-year quality metrics written")
