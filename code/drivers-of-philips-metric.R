# ──────────────────────────────────────────────────────────────────────────────
# Drivers of ANACONDA "detail" — Age & UCOD mix from PARQUET (hardened joins)
# Date: 2025-09-11
# ──────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(dplyr); library(readr); library(stringr); library(purrr); library(tidyr); library(tibble)
    library(DBI); library(duckdb); library(glue); library(ggplot2); library(scales)
    library(here); library(fs); library(rlang); library(glmnet); library(broom)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b
log <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "│", paste0(..., collapse=""), "\n")

# ── Paths & params ────────────────────────────────────────────────────────────
out_dir   <- here("output")
diag_dir  <- file.path(out_dir, "cstd_driver_analysis")
fig_dir   <- file.path(diag_dir, "figures")
dir_create(diag_dir); dir_create(fig_dir)

parquet_dir <- if (dir_exists(here("data_private","mcod"))) here("data_private","mcod") else here("data_private","mcod_sample")
stopifnot(dir_exists(parquet_dir))

cluster_metrics_path <- file.path(out_dir, "cluster_metrics.csv.gz")
cluster_members_path <- file.path(out_dir, "county_cluster_membership.csv.gz")
stopifnot(file.exists(cluster_metrics_path), file.exists(cluster_members_path))

# Garbage & mapping resources
garbage_csv   <- here("data_raw","cause-codes","gbd_garbage_codes_without_overdose.csv")
GBD_ROOT3_CSV <- here("data_raw","icd10_root3_to_gbdl3.csv")
stopifnot(file.exists(garbage_csv), file.exists(GBD_ROOT3_CSV))

# Periods must match your ANACONDA script
periods <- list("1999_2006" = 1999:2006, "2007_2014" = 2007:2014, "2015_2022" = 2015:2022)

# Age bins
age_breaks <- c(-Inf, 1, 5, 15, 35, 55, 75, Inf)
age_labels <- c("age_u1", "age_1_4", "age_5_14", "age_15_34", "age_35_54", "age_55_74", "age_75p")

# ── Helpers ───────────────────────────────────────────────────────────────────
safe_read <- function(path, ...) if (!file.exists(path)) tibble() else suppressMessages(readr::read_csv(path, show_col_types = FALSE, ...))
scale_01  <- function(x) if (sd(x %||% 0, na.rm=TRUE) > 0) (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE) else x*0

ensure_keys <- function(df, keys = c("period","cluster")) {
    # Ensure df has key columns (character) even if empty; keeps joins from erroring
    for (k in keys) if (!k %in% names(df)) df[[k]] <- character()
    df %>% mutate(across(all_of(keys), as.character))
}

build_xy <- function(df, y_col, x_cols) {
    df2 <- df %>% select(all_of(c("period","cluster", y_col, x_cols))) %>%
        filter(is.finite(.data[[y_col]])) %>% drop_na(all_of(x_cols))
    list(X = as.matrix(df2 %>% select(all_of(x_cols))), y = as.numeric(df2[[y_col]]), df = df2)
}

fit_enet <- function(X, y, family="gaussian", alpha=0.9, nfolds=5, seed=123) {
    set.seed(seed); if (NROW(X) < 10 || NCOL(X) < 2) return(NULL)
    tryCatch(cv.glmnet(X, y, family=family, alpha=alpha, standardize=TRUE, nfolds=min(nfolds, NROW(X))), error=function(e) NULL)
}

coef_table_from_enet <- function(fit, feature_names) {
    if (is.null(fit)) return(tibble(feature=character(), coef=numeric(), lambda=""))
    lambda <- fit$lambda.min %||% fit$lambda.1se
    cfs <- as.numeric(coef(fit, s=lambda))[-1]
    tibble(feature = feature_names, coef = cfs, lambda = ifelse(identical(lambda, fit$lambda.min), "lambda.min", "lambda.1se"))
}

perm_importance <- function(fit, X, y, feature_names, B=100, seed=123) {
    if (is.null(fit)) return(tibble(feature=character(), rel_mse=double()))
    set.seed(seed)
    lambda <- fit$lambda.min %||% fit$lambda.1se
    yhat0 <- drop(predict(fit, newx=X, s=lambda))
    mse0  <- mean((y - yhat0)^2, na.rm=TRUE); if (!is.finite(mse0) || mse0 == 0) mse0 <- 1e-8
    map_dfr(seq_along(feature_names), function(j) {
        mse_j <- replicate(B, { Xp <- X; Xp[,j] <- sample(Xp[,j], FALSE); mean((y - drop(predict(fit, newx=Xp, s=lambda)))^2, na.rm=TRUE) })
        tibble(feature = feature_names[j], rel_mse = mean(mse_j) / mse0)
    }) %>% arrange(desc(rel_mse))
}

icd_root3 <- function(code) {
    x <- toupper(gsub("[^A-Z0-9]", "", as.character(code)))
    m <- stringr::str_match(x, "^([A-Z])([0-9])([0-9A-Z])")
    ifelse(is.na(m[,1]), NA_character_, paste0(m[,2], m[,3], m[,4]))
}

# ── Load cluster metrics & membership ─────────────────────────────────────────
clu_metrics <- safe_read(cluster_metrics_path) %>%
    rename(
        detail_ucod_root3 = dplyr::any_of("detail_ucod_root3"),
        detail_ucod_gbdl3 = dplyr::any_of("detail_ucod_gbdl3"),
        detail_mcod       = dplyr::any_of("detail_mcod")
    ) %>%
    mutate(period = as.character(period))

clu_members <- safe_read(cluster_members_path) %>%
    transmute(period = as.character(period),
              cluster = as.character(cluster),
              fips = stringr::str_pad(as.character(fips), 5, pad="0"),
              cluster_deaths = as.numeric(cluster_deaths))

stopifnot(nrow(clu_metrics) > 0, nrow(clu_members) > 0)

# ── Load garbage roots & GBD map ─────────────────────────────────────────────
garbage_root_set <- readr::read_csv(garbage_csv, show_col_types = FALSE) %>%
    transmute(root3 = icd_root3(icd10)) %>% filter(!is.na(root3)) %>% distinct() %>% pull(root3)

gbd_map_df <- readr::read_csv(GBD_ROOT3_CSV, show_col_types = FALSE) %>%
    transmute(root3 = toupper(substr(root3, 1, 3)), gbd_l3 = as.character(gbd_l3)) %>%
    filter(!is.na(root3), nzchar(root3), !is.na(gbd_l3), nzchar(gbd_l3)) %>%
    distinct(root3, .keep_all = TRUE)

stopifnot(!anyDuplicated(gbd_map_df$root3))

# ── AGE COMPOSITION from parquet ──────────────────────────────────────────────
log("Computing age composition from parquet…")

build_age_comp_for_period <- function(pname, years) {
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir=":memory:"); on.exit(try(DBI::dbDisconnect(con, shutdown=TRUE), silent=TRUE), add=TRUE)
    DBI::dbExecute(con, "PRAGMA threads = 4")
    duckdb::duckdb_register(con, "map",
                            clu_members %>% filter(period == pname) %>% select(fips, cluster) %>% distinct()
    )
    
    sql <- glue("
    WITH cert AS (
      SELECT
        LPAD(CAST(county_ihme AS VARCHAR), 5, '0') AS fips,
        CAST(age_years AS DOUBLE) AS age_years,
        ager27,
        ucod
      FROM parquet_scan('{normalizePath(parquet_dir, winslash = '/') }/*.parquet')
      WHERE year BETWEEN {min(years)} AND {max(years)} AND ucod IS NOT NULL
    )
    SELECT m.cluster, c.fips,
           CASE
             WHEN age_years IS NOT NULL THEN age_years
             WHEN ager27 IS NOT NULL THEN
               CASE
                 WHEN ager27 = 1 THEN 0.0
                 WHEN ager27 = 2 THEN 0.5
                 WHEN ager27 BETWEEN 3 AND 4 THEN 3.0
                 WHEN ager27 BETWEEN 5 AND 6 THEN 10.0
                 WHEN ager27 BETWEEN 7 AND 10 THEN 25.0
                 WHEN ager27 BETWEEN 11 AND 12 THEN 40.0
                 WHEN ager27 BETWEEN 13 AND 15 THEN 60.0
                 WHEN ager27 BETWEEN 16 AND 17 THEN 70.0
                 WHEN ager27 BETWEEN 18 AND 19 THEN 80.0
                 WHEN ager27 BETWEEN 20 AND 22 THEN 85.0
                 WHEN ager27 BETWEEN 23 AND 25 THEN 90.0
                 WHEN ager27 = 26 THEN 95.0
                 ELSE NULL
               END
             ELSE NULL
           END AS age_est
    FROM cert c JOIN map m ON c.fips = m.fips
    WHERE m.cluster IS NOT NULL
  ")
    
    df <- DBI::dbGetQuery(con, sql) %>% as_tibble()
    if (!nrow(df)) return(tibble(period = pname, cluster = character()))
    
    df %>%
        mutate(age_bin = cut(age_est, breaks = age_breaks, labels = age_labels, right = FALSE)) %>%
        filter(!is.na(age_bin)) %>%
        count(cluster, age_bin, name = "n") %>%
        group_by(cluster) %>%
        mutate(share = n / sum(n)) %>%
        ungroup() %>%
        pivot_wider(names_from = age_bin, values_from = share, values_fill = 0) %>%
        mutate(period = pname, .before = 1)
}

age_comp <- map2_dfr(names(periods), periods, build_age_comp_for_period) %>% ensure_keys()
readr::write_csv(age_comp, file.path(diag_dir, "drivers_age_shares.csv"))

# ── UCOD COMPOSITION from parquet (root-3 & GBD L3) ──────────────────────────
log("Computing UCOD composition (root-3 and GBD L3) from parquet…")

build_ucod_comp_for_period <- function(pname, years) {
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir=":memory:"); on.exit(try(DBI::dbDisconnect(con, shutdown=TRUE), silent=TRUE), add=TRUE)
    DBI::dbExecute(con, "PRAGMA threads = 4")
    
    duckdb::duckdb_register(con, "map",
                            clu_members %>% filter(period == pname) %>% select(fips, cluster) %>% distinct()
    )
    duckdb::duckdb_register(con, "garbage", tibble(root3 = garbage_root_set))
    duckdb::duckdb_register(con, "gbd_root3_map", gbd_map_df)
    
    DBI::dbExecute(con, glue("
    CREATE TEMP TABLE ucod_valid AS
    SELECT
      LPAD(CAST(county_ihme AS VARCHAR), 5, '0') AS fips,
      SUBSTR(UPPER(regexp_replace(ucod, '[^A-Za-z0-9]', '')),1,3) AS root3
    FROM parquet_scan('{normalizePath(parquet_dir, winslash = '/') }/*.parquet')
    WHERE year BETWEEN {min(years)} AND {max(years)} AND ucod IS NOT NULL
  "))
    
    DBI::dbExecute(con, "
    CREATE TEMP TABLE ucod_clean AS
    SELECT fips, root3
    FROM ucod_valid
    WHERE root3 BETWEEN 'A00' AND 'Z99'
      AND SUBSTR(root3,2,1) BETWEEN '0' AND '9'
      AND root3 NOT IN (SELECT root3 FROM garbage)
  ")
    
    root3_tbl <- DBI::dbGetQuery(con, "
    SELECT m.cluster, c.root3, COUNT(*) AS n
    FROM ucod_clean c JOIN map m ON c.fips = m.fips
    WHERE m.cluster IS NOT NULL
    GROUP BY 1,2
  ") %>% as_tibble()
    
    gbd_tbl <- DBI::dbGetQuery(con, "
    SELECT m.cluster, g.gbd_l3 AS domain_code, COUNT(*) AS n
    FROM ucod_clean c
    JOIN gbd_root3_map g ON c.root3 = g.root3
    JOIN map m ON c.fips = m.fips
    WHERE m.cluster IS NOT NULL
    GROUP BY 1,2
  ") %>% as_tibble()
    
    root3_shares <- if (nrow(root3_tbl)) {
        root3_tbl %>% group_by(cluster) %>% mutate(share = n / sum(n)) %>% ungroup() %>%
            transmute(cluster, domain_code = root3, share) %>% mutate(period = pname, .before = 1)
    } else tibble(period = pname, cluster = character(), domain_code = character(), share = numeric())
    
    gbd_shares <- if (nrow(gbd_tbl)) {
        gbd_tbl %>% group_by(cluster) %>% mutate(share = n / sum(n)) %>% ungroup() %>%
            transmute(cluster, domain_code, share) %>% mutate(period = pname, .before = 1)
    } else tibble(period = pname, cluster = character(), domain_code = character(), share = numeric())
    
    list(root3 = root3_shares, gbd = gbd_shares)
}

ucod_comp_list   <- map2(names(periods), periods, build_ucod_comp_for_period)
ucod_root3_shares<- map_dfr(ucod_comp_list, "root3") %>% ensure_keys()
ucod_gbdl3_shares<- map_dfr(ucod_comp_list, "gbd")   %>% ensure_keys()

readr::write_csv(ucod_root3_shares, file.path(diag_dir, "ucod_root3_shares.csv"))
readr::write_csv(ucod_gbdl3_shares, file.path(diag_dir, "ucod_gbdl3_shares.csv"))

# ── Build analysis frames ─────────────────────────────────────────────────────
df_base <- safe_read(cluster_metrics_path) %>%
    mutate(period = as.character(period),
           cluster = as.character(cluster)) %>%
    select(period, cluster, detail_ucod_root3, detail_ucod_gbdl3, detail_mcod) %>%
    distinct() %>%
    ensure_keys()

# Age wide
df_age <- age_comp %>% ensure_keys()

# UCOD wide (root-3 & GBD L3)
cause_wide <- function(df, prefix="uc_") {
    if (!nrow(df)) return(tibble(period=character(), cluster=character()))
    df %>%
        mutate(code_col = paste0(prefix, str_replace_all(domain_code, "[^A-Za-z0-9]+", "_"))) %>%
        select(period, cluster, code_col, share) %>%
        distinct() %>%
        pivot_wider(names_from = code_col, values_from = share, values_fill = 0) %>%
        ensure_keys()
}

ucod_r3_wide <- cause_wide(ucod_root3_shares, prefix="uc_r3_")
ucod_gb_wide <- cause_wide(ucod_gbdl3_shares, prefix="uc_gb_")

df_all_r3 <- df_base %>% left_join(df_age, by=c("period","cluster")) %>% left_join(ucod_r3_wide, by=c("period","cluster"))
df_all_gb <- df_base %>% left_join(df_age, by=c("period","cluster")) %>% left_join(ucod_gb_wide, by=c("period","cluster"))

# ── Quick descriptive associations (age ↔ detail) ────────────────────────────
plot_age_vs_detail <- function(df, ycol, fname_prefix) {
    if (!nrow(df) || !ycol %in% names(df)) return(invisible(NULL))
    for (acol in age_labels) {
        if (!acol %in% names(df)) next
        p <- ggplot(df, aes(x=.data[[acol]], y=.data[[ycol]])) +
            geom_point(alpha=0.5) +
            geom_smooth(method="loess", se=TRUE) +
            scale_x_continuous(labels=percent) +
            labs(x=paste0("Share ", acol), y=ycol, title=paste("Detail vs", acol),
                 subtitle="Each point is a cluster×period") +
            theme_minimal(base_size = 12)
        ggsave(file.path(fig_dir, sprintf("%s_%s_vs_%s.png", fname_prefix, ycol, acol)), p, width=7, height=5, dpi=200)
    }
}

invisible(plot_age_vs_detail(df_all_r3, "detail_ucod_root3", "age_detail"))
invisible(plot_age_vs_detail(df_all_r3, "detail_mcod",       "age_detail"))
invisible(plot_age_vs_detail(df_all_gb, "detail_ucod_gbdl3", "age_detail"))

# ── Modeling: ENet + OLS ─────────────────────────────────────────────────────
build_and_run <- function(df, y_col, family_label, ucod_prefix) {
    if (!nrow(df) || !y_col %in% names(df)) return(invisible(TRUE))
    
    age_x <- intersect(age_labels, names(df))
    xy_age <- build_xy(df, y_col, age_x)
    
    coefs_age <- imp_age <- tibble(); fit_age_en <- NULL
    if (length(age_x) >= 2 && NROW(xy_age$X) >= 10) {
        fit_age_en <- fit_enet(xy_age$X, xy_age$y)
        coefs_age  <- coef_table_from_enet(fit_age_en, age_x) %>% mutate(model="ENet_age", target=y_col, family=family_label)
        imp_age    <- perm_importance(fit_age_en, xy_age$X, xy_age$y, age_x, B=50) %>% mutate(model="ENet_age", target=y_col, family=family_label)
    }
    
    ols_age <- tryCatch({
        dd <- xy_age$df
        if (nrow(dd) >= length(age_x) + 3) broom::tidy(lm(reformulate(age_x, y_col), data=dd)) %>%
            mutate(model="OLS_age", target=y_col, family=family_label) else tibble()
    }, error=function(e) tibble())
    
    cause_cols <- names(df)[grepl(paste0("^", ucod_prefix), names(df))]
    x_cols <- unique(c(age_x, cause_cols))
    xy_all <- build_xy(df, y_col, x_cols)
    
    coefs_all <- imp_all <- tibble()
    if (length(x_cols) >= 5 && NROW(xy_all$X) >= 20) {
        fit_all_en <- fit_enet(xy_all$X, xy_all$y)
        coefs_all  <- coef_table_from_enet(fit_all_en, x_cols) %>% mutate(model="ENet_age+cause", target=y_col, family=family_label)
        imp_all    <- perm_importance(fit_all_en, xy_all$X, xy_all$y, x_cols, B=50) %>% mutate(model="ENet_age+cause", target=y_col, family=family_label)
        
        if (nrow(imp_all)) {
            imp_top <- imp_all %>% arrange(desc(rel_mse)) %>% group_by(target, family, model) %>% slice_head(n=20) %>% ungroup()
            readr::write_csv(imp_top, file.path(diag_dir, sprintf("model_importance_%s_%s.csv", family_label, y_col)))
            p <- ggplot(imp_top, aes(x=reorder(feature, rel_mse), y=rel_mse)) +
                geom_col() + coord_flip() +
                labs(x="Feature", y="Relative MSE ↑ (permute feature)", title=sprintf("Permutation importance: %s ~ age + causes", y_col), subtitle=family_label) +
                theme_minimal(base_size = 12)
            ggsave(file.path(fig_dir, sprintf("importance_%s_%s.png", family_label, y_col)), p, width=7.5, height=6, dpi=200)
        }
    }
    
    coefs_out <- bind_rows(coefs_age, coefs_all, ols_age)
    if (nrow(coefs_out)) readr::write_csv(coefs_out, file.path(diag_dir, sprintf("model_coefs_%s_%s.csv", family_label, y_col)))
    
    invisible(TRUE)
}

build_and_run(df_all_r3, "detail_ucod_root3", family_label="ucod_root3", ucod_prefix="uc_r3_")
build_and_run(df_all_r3, "detail_mcod",       family_label="mcod_root3", ucod_prefix="uc_r3_")
if ("detail_ucod_gbdl3" %in% names(df_all_gb)) {
    build_and_run(df_all_gb, "detail_ucod_gbdl3", family_label="ucod_gbdl3", ucod_prefix="uc_gb_")
}

# ── Compact human-readable summaries ─────────────────────────────────────────
log("Writing compact summaries…")

assoc_from_shares <- function(share_df, detail_col, out_csv) {
    df <- share_df %>% ensure_keys() %>% left_join(df_base, by=c("period","cluster"))
    if (!nrow(df) || !detail_col %in% names(df)) return(tibble())
    out <- df %>%
        filter(!is.na(share), is.finite(share), is.finite(.data[[detail_col]])) %>%
        group_by(domain_code) %>%
        summarise(cor_with_detail = {
            z <- suppressWarnings(cor(share, .data[[detail_col]], use="complete.obs"))
            ifelse(is.finite(z), z, NA_real_)
        },
        n = dplyr::n(), .groups="drop") %>%
        filter(is.finite(cor_with_detail)) %>%
        arrange(desc(abs(cor_with_detail)))
    if (nrow(out)) readr::write_csv(out, file.path(diag_dir, out_csv))
    out
}

assoc_r3 <- assoc_from_shares(ucod_root3_shares, "detail_ucod_root3", "ucond_root3_association_with_detail.csv")
if (nrow(assoc_r3)) {
    p <- ggplot(assoc_r3 %>% slice_head(n=15), aes(x=reorder(domain_code, cor_with_detail), y=cor_with_detail)) +
        geom_col() + coord_flip() +
        labs(x="UCOD root-3", y="Correlation with detail (cluster×period)", title="Underlying causes most associated with detail (root-3)") +
        theme_minimal(base_size = 12)
    ggsave(file.path(fig_dir, "ucond_root3_assoc_bar.png"), p, width=7.5, height=6, dpi=200)
}

if ("detail_ucod_gbdl3" %in% names(df_base)) {
    assoc_gb <- assoc_from_shares(ucod_gbdl3_shares, "detail_ucod_gbdl3", "ucond_gbdl3_association_with_detail.csv")
    if (nrow(assoc_gb)) {
        p <- ggplot(assoc_gb %>% slice_head(n=15), aes(x=reorder(domain_code, cor_with_detail), y=cor_with_detail)) +
            geom_col() + coord_flip() +
            labs(x="GBD L3", y="Correlation with detail (cluster×period)", title="Underlying causes most associated with detail (GBD L3)") +
            theme_minimal(base_size = 12)
        ggsave(file.path(fig_dir, "ucond_gbdl3_assoc_bar.png"), p, width=7.5, height=6, dpi=200)
    }
}

log("Done. See: ", diag_dir)

