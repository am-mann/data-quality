# --------------------------------------------------------------
# demo_quality_timeseries.R
# Memory-efficient demographic time-series + separate plots
# Includes overall composite z-score
# --------------------------------------------------------------

# ── 0 · Packages ──────────────────────────────────────────────
library(dplyr);  library(tidyr);  library(purrr);  library(stringr)
library(arrow);  library(readr);  library(here)
library(ggplot2); library(scales)

# ── 1 · Paths & dictionaries ─────────────────────────────────
parquet_dir <- if (dir.exists(here("data_private", "mcod"))) {
    here("data_private", "mcod")
} else {
    here("data_private", "mcod_sample")
}
# ---- education recode ----------------------------------------------------
recode_educ <- function(x) {
    x <- suppressWarnings(as.numeric(x))   # allow character or numeric input
    
    dplyr::case_when(
        # ── 2003-revision categorical codes (1–8) ─────────────────────────────
        x %in% c(1, 2)                   ~ "<HS",              # ≤HS, no diploma
        x == 3                           ~ "HS",               # HS / GED
        x == 4                           ~ "Some College",     # <1 yr college
        x == 5                           ~ "Associate Degree", # 1–3 yrs / AA
        x == 6                           ~ "Bachelor",
        x %in% c(7, 8)                   ~ "Post-Bachelor",    # master / prof / PhD
        
        # ── Legacy “years of schooling” codes (0–17+) ────────────────────────
        x <= 11                          ~ "<HS",
        x == 12                          ~ "HS",
        x == 13                          ~ "Some College",     # 1 yr college
        x == 14                          ~ "Associate Degree", # 2 yrs college
        x == 15                          ~ "Some College",     # 3 yrs college
        x == 16                          ~ "Bachelor",
        x >= 17                          ~ "Post-Bachelor",
        
        # ── House-keeping ─────────────────────────────────────────────────────
        is.na(x) | x %in% c(9, 99)       ~ "Unknown",
        TRUE                             ~ "Other"
    )
}


dict_dir   <- here("data_raw", "cause-codes")
sec_cols   <- paste0("record_", 1:20)
demo_vars  <- c("educ", "sex", "race", "marstat")   # age_group added later

canonical_icd <- function(x) stringr::str_remove_all(stringr::str_to_upper(x), "[^A-Z0-9]")

recode_race <- function(x) {
    dplyr::case_when(
        x %in% c("1", "01", "WHITE", "W")                      ~ "White",
        x %in% c("2", "02", "BLACK", "B")                      ~ "Black",
        x %in% c("3", "03", "AMERICAN INDIAN", "AIAN")         ~ "Native American",
        x %in% c("4", "05", "HISPANIC", "MEXICAN", "PUERTO RICAN",
                 "CUBAN", "OTHER SPANISH", "H")               ~ "Hispanic",
        TRUE                                                   ~ "Other"
    )
}

lookup_garbage <- read_csv(file.path(dict_dir,
                                     "gbd_garbage_codes_without_overdose.csv"),
                           show_col_types = FALSE) |>
    transmute(icd10 = canonical_icd(icd10),
              gbd_severity = as.integer(gbd_severity))

light_garbage <- read_csv(file.path(dict_dir, "light_garbage_codes.csv"),
                          show_col_types = FALSE)$code |> canonical_icd()

overdose_ucod <- read_csv(file.path(dict_dir,
                                    "overdose-detail-codes-v3.csv"),
                          show_col_types = FALSE)$ucod |> canonical_icd()

accident_detail_lookup <- read_csv(
    file.path(dict_dir, "accident-detail-codes-v3.csv"),
    show_col_types = FALSE) |>
    pivot_longer(starts_with("detail_"), values_to = "detail",
                 values_drop_na = TRUE) |>
    transmute(
        ucod   = canonical_icd(str_trim(ucod)),
        detail = canonical_icd(str_trim(detail))
    ) |> distinct()

accident_ucod_set   <- unique(accident_detail_lookup$ucod)
accident_detail_set <- unique(accident_detail_lookup$detail)

# ── 2 · Helper: record-level flags for ONE file ───────────────
prepare_year_records <- function(file) {
    yr <- as.integer(stringr::str_extract(basename(file), "\\d{4}"))
    df <- read_parquet(file, as_data_frame = TRUE) |>
        select(-any_of("injwork")) |>                    # drop mixed-type col
        mutate(
            year  = yr,
            ranum = row_number(),
            icd10 = canonical_icd(ucod),
            across(starts_with("record"), canonical_icd)
        )
    
    ## contributing causes (long)
    present_sec <- intersect(sec_cols, names(df))
    contrib_all <- if (!length(present_sec)) {
        tibble(ranum = integer(), detail_full = character(), detail = character())
    } else {
        present_sec |>
            map_dfr(~ df |> select(ranum, detail_full = !!sym(.x))) |>
            filter(!is.na(detail_full) & trimws(detail_full) != "") |>
            mutate(detail_full = canonical_icd(detail_full),
                   detail      = substr(detail_full, 1, 4))
    }
    
    contrib_counts <- contrib_all |>
        mutate(is_light = detail_full %in% light_garbage) |>
        group_by(ranum) |>
        summarise(total_contrib = n(),
                  light_contrib = sum(is_light), .groups = "drop")
    
    required_hits <- contrib_all |>
        filter(detail %in% accident_detail_set) |>
        left_join(df |> select(ranum, icd10), by = "ranum") |>
        semi_join(accident_detail_lookup, by = c("icd10" = "ucod", "detail")) |>
        distinct(ranum) |> mutate(has_required_detail = TRUE)
    
    unspec_tbl <- contrib_all |>
        filter(str_starts(detail_full, "T")) |>
        mutate(is_specific = !(str_starts(detail_full, "T509") |
                                   str_starts(detail_full, "T409"))) |>
        group_by(ranum) |>
        summarise(unspecific_drug = !any(is_specific), .groups = "drop")
    
    ## merge back
    df |>
        left_join(lookup_garbage, by = "icd10") |>
        mutate(
            needs_detail = icd10 %in% accident_ucod_set,
            is_overdose  = icd10 %in% overdose_ucod
        ) |>
        left_join(contrib_counts, by = "ranum") |>
        left_join(unspec_tbl,     by = "ranum") |>
        left_join(required_hits,  by = "ranum") |>
        mutate(
            total_contrib       = replace_na(total_contrib, 0L),
            light_contrib       = replace_na(light_contrib, 0L),
            unspecific_drug     = replace_na(unspecific_drug, TRUE),
            has_required_detail = replace_na(has_required_detail, FALSE),
            is_garbage = case_when(
                is.na(gbd_severity)                     ~ FALSE,
                !str_starts(icd10, "X")                 ~ TRUE,
                str_starts(icd10, "X") & unspecific_drug~ TRUE,
                TRUE                                    ~ FALSE
            )
        )
}

# ── 3 · Helper: summarise ONE demographic variable ────────────
summarise_demo <- function(df, var) {
    grp <- sym(var)
    
    df <- df %>%
        mutate(age_group = cut(
            age_years,
            breaks = c(-Inf, 1, 5, 15, 25, 45, 65, Inf),
            right  = FALSE,
            labels = c("<1","1-4","5-14","15-24","25-44","45-64","65+")
        ))
    
    if (var == "race") {
        df <- mutate(df, race = recode_race(as.character(race)))
        
    } else if (var == "educ") {
        # Use modern 'educ2003' when it exists, otherwise the legacy 'educ' field
        df <- mutate(df,
                     educ = recode_educ(as.numeric(educ)))
        
    } else if (var != "age_group") {
        df <- mutate(df, {{ grp }} := as.character({{ grp }}))
    }
    
    df |>
        group_by(year, demo_level = .data[[var]]) |>
        summarise(
            demo_var       = var,
            n_cert         = n(),
            prop_garbage   = sum(is_garbage) /
                n_cert,
            prop_light     = mean(total_contrib > 0 & light_contrib > 0, na.rm = TRUE),
            pct_acc_miss   = mean(needs_detail & !has_required_detail, na.rm = TRUE),
            pct_overd_miss = mean(is_overdose  & !has_required_detail, na.rm = TRUE),
            .groups = "drop"
        )
}

# ── 4 · Stream through files, build summaries ─────────────────
demo_year_metrics <- purrr::map_dfr(
    list.files(parquet_dir, "\\.parquet$", full.names = TRUE),
    function(file) {
        message("Processing ", basename(file))
        df <- prepare_year_records(file)
        
        bind_rows(
            summarise_demo(df, "age_group"),
            purrr::map_dfr(demo_vars, ~ summarise_demo(df, .x))
        )
    }
)

# ── 5 · Within-year z-scores + composite z-score ──────────────
demo_year_metrics <- demo_year_metrics |>
    group_by(year, demo_var) |>
    mutate(
        z_prop_garbage  = scale(prop_garbage)[, 1],
        z_prop_light    = scale(prop_light)[, 1],
        z_pct_acc_miss  = scale(pct_acc_miss)[, 1],
        z_pct_overd_miss= scale(pct_overd_miss)[, 1],
        overall_z_score = rowMeans(across(c(z_prop_garbage,
                                            z_prop_light,
                                            z_pct_acc_miss,
                                            z_pct_overd_miss)), na.rm = TRUE)
    ) |>
    ungroup()

write_csv(demo_year_metrics,
          here("data", "demographic_year_quality_metrics.csv"))

# ── 6 · Separate time-series plots for each metric × demo_var ─
if (!dir.exists(here("figures"))) dir.create(here("figures"), recursive = TRUE)

metrics_to_plot <- c("prop_garbage", "prop_light",
                     "pct_acc_miss", "pct_overd_miss", "overall_z_score")

demo_vars_plot <- unique(demo_year_metrics$demo_var)   # includes age_group

plot_metric_demo <- function(metric, dv) {
    g <- ggplot(
        filter(demo_year_metrics, demo_var == dv, demo_level != "Unknown"),
        aes(x = year, y = .data[[metric]], colour = demo_level)
    ) +
        geom_line(linewidth = 1) +
        labs(title = paste("Time-series of", metric, "by", dv),
             x = NULL, y = NULL, colour = NULL) +
        theme_bw(base_size = 11)
    
    # use % scale for proportion / pct metrics
    if (grepl("^prop_|^pct_", metric)) {
        g <- g + scale_y_continuous(labels = percent_format(accuracy = 0.1))
    }
    
    g
}

purrr::walk(metrics_to_plot, function(met) {
    purrr::walk(demo_vars_plot, function(dv) {
        ggsave(
            filename = here("figures",
                            paste0("ts_", met, "_", dv, ".png")),
            plot     = plot_metric_demo(met, dv),
            width    = 9, height = 6, dpi = 320
        )
    })
})

# ── Top-100 garbage UCODs -----------------------------------------------
library(arrow); library(dplyr); library(purrr); library(readr)

top100_garbage <- list.files(parquet_dir, "\\.parquet$", full.names = TRUE) |>
    map_dfr(function(file) {
        message("Scanning ", basename(file))
        # read **only** the UCOD column, then canonicalise
        read_parquet(file, columns = "ucod") |>
            transmute(icd10 = canonical_icd(ucod))
    }) |>
    filter(icd10 %in% lookup_garbage$icd10) |>      # keep garbage codes only
    count(icd10, sort = TRUE, name = "n_cert") |>
    left_join(lookup_garbage, by = "icd10") |>      # optional: add severity
    slice_head(n = 100)                             # keep top 100

# Inspect
print(top100_garbage, n = 20)

# Or write to disk for later use
write_csv(top100_garbage, here("data", "top100_garbage_codes.csv"))
