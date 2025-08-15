# --------------------------------------------------------------
# ANACONDA detail metric
# Date: August 8, 2025
#
# This script is an implementation of "detail" metric from the Phillips ANACONDA paper.
# Because most counties have far too few deaths for this metric, they were put into clusters 
# of at least 2000 by merging proximate counties and averaging overtime larger time periods.
# The code computes Shannon entropy over GBD cause categories and converts them to the 
# effective number of causes of death. It then standardizes to 2000 deaths by interpolating clusters
# with more than 2000 deaths. The script uses greedy clustering which means that there are sometimes counties not put in any cluster.
# This is fixed by adding those counties to the nearest cluster
#
# The only difference between this metric and Phillips' is that they use
# reference size as the sample size at which all populations reach 95% sample completeness. 
# Unfortunately, our clusters are too small for that to be feasible, but plots (see slides) 
# show that large cluster size effects seem to drop-off around 2000 deaths making it a justifiable reference size. 
# Phillips et al note in their paper that size based rarefaction is appropriate for small places. 
# --------------------------------------------------------------

suppressPackageStartupMessages({
    library(fs); library(readr); library(dplyr); library(stringr)
    library(purrr); library(tidyr); library(tibble)
    library(sf); library(tigris); library(igraph)
    library(here); library(iNEXT)
    library(DBI); library(duckdb); library(glue)
})

# -------------------- helpers & params --------------------
log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "│", paste0(..., collapse=""), "\n")
log_err <- function(prefix, e) { cat(format(Sys.time(), "%H:%M:%S"), "│ ERROR │", prefix, "│", conditionMessage(e), "\n") }

# Controls
FIXED_K_GLOBAL             <- TRUE   # <-- Stabilized K on
CLAMP_0_100                <- TRUE
EXCLUDE_UCOD_FROM_MCOD     <- TRUE   # drop UCOD roots from MCOD set
AUTO_DETECT_RECORD1_IS_UCOD<- TRUE
INCLUDE_UCOD_GBDL3         <- TRUE   # compute detail on GBD L3

.clamp100 <- function(x) { if (!CLAMP_0_100) return(x); y <- x; y[x < 0] <- 0; y[x > 100] <- 100; y }

# Paths
out_dir  <- here("output"); dir_create(out_dir)
diag_dir <- file.path(out_dir, "diag"); dir_create(diag_dir)
parquet_dir <- if (dir_exists(here("data_private","mcod"))) here("data_private","mcod") else here("data_private","mcod_sample")

garbage_csv   <- here("data_raw","cause-codes","gbd_garbage_codes_without_overdose.csv")
GBD_ROOT3_CSV <- here("data_raw","lookup_icd10root_to_gbd_lvl3_lvl2.csv")
GBD_ROOT3_COL <- "icd10_root"
GBD_L3_COL    <- "gbd_cause_lvl3"

stopifnot(file.exists(garbage_csv))
stopifnot(file.exists(GBD_ROOT3_CSV))

# Other params
county_var  <- "county_ihme"
min_deaths  <- 2000L
periods <- list("1999_2006" = 1999:2006, "2007_2014" = 2007:2014, "2015_2022" = 2015:2022)
years_all <- range(unlist(periods))
M_REF   <- 2000L  # standardization size (interpolate large clusters down to 2000 only)
crs_proj <- 5070
threads  <- max(1L, parallel::detectCores() - 1L)

icd_root3 <- function(code) {
    x <- toupper(gsub("[^A-Z0-9]", "", as.character(code)))
    m <- stringr::str_match(x, "^([A-Z])([0-9])([0-9A-Z])")
    ifelse(is.na(m[,1]), NA_character_, paste0(m[,2], m[,3], m[,4]))
}

extract_q1 <- function(est) {
    if (inherits(est, "try-error") || is.null(est)) return(NA_real_)
    oq  <- suppressWarnings(as.numeric(est$Order.q))
    col <- intersect(c("Estimator","qD"), names(est))
    val <- suppressWarnings(as.numeric(est[[col[1]]]))
    out <- val[which(oq == 1)]
    ifelse(length(out) && is.finite(out) && out > 0, out, NA_real_)
}

# Phillips detail at stabilized maxH; interpolate only (m_use = min(m_ref, n))
compute_detail_row <- function(tbl_counts, max_H_global, m_ref = M_REF) {
    ftab <- as.numeric(tbl_counts$n); names(ftab) <- tbl_counts$domain_code
    ftab <- ftab[is.finite(ftab) & ftab > 0]
    n <- sum(ftab); S <- length(ftab); f1 <- sum(ftab == 1); f2 <- sum(ftab == 2)
    if (!is.finite(n) || n < 2 || S < 2 || !is.finite(max_H_global) || max_H_global <= 0) {
        return(tibble(
            deaths_no_garbage = n, S = S, f1 = f1, f2 = f2,
            H_raw = NA_real_, detail_phillips_raw = NA_real_, detail_phillips_refsize = NA_real_,
            D1_mref = NA_real_, D1_eff_used = NA_real_,
            m_ref = m_ref, eff_m_used = NA_integer_, used_extrapolation = FALSE, fallback_reason = "degenerate"
        ))
    }
    p <- ftab / n
    H_raw <- -sum(p * log(p))
    detail_raw <- .clamp100(100 * H_raw / max_H_global)
    m_use <- as.integer(min(m_ref, n))  # interpolate only; no extrapolation
    est   <- try(iNEXT::estimateD(ftab, datatype = "abundance", base = "size", level = m_use, conf = FALSE), silent = TRUE)
    D1_use <- extract_q1(est)
    if (is.na(D1_use) && is.finite(H_raw)) D1_use <- exp(H_raw)
    tibble(
        deaths_no_garbage = n, S = S, f1 = f1, f2 = f2, H_raw = H_raw,
        detail_phillips_raw     = detail_raw,
        detail_phillips_refsize = if (is.finite(D1_use) && D1_use > 0) .clamp100(100 * log(D1_use) / max_H_global) else NA_real_,
        D1_mref = if (m_ref <= n) D1_use else NA_real_,
        D1_eff_used = D1_use,
        m_ref = as.integer(m_ref), eff_m_used = m_use, used_extrapolation = FALSE, fallback_reason = "ok"
    )
}

# Entropy contribution diagnostics (per cluster × domain)
contrib_from_counts <- function(df_counts, value_col = c("n","incidence")) {
    value_col <- match.arg(value_col)
    df_counts %>%
        group_by(cluster) %>%
        mutate(
            total = sum(.data[[value_col]], na.rm = TRUE),
            p     = ifelse(total > 0, .data[[value_col]] / total, NA_real_),
            contr_nats = ifelse(is.finite(p) & p > 0, -p * log(p), 0)
        ) %>%
        group_by(cluster) %>%
        mutate(H_cluster = sum(contr_nats, na.rm = TRUE),
               contr_frac = ifelse(H_cluster > 0, contr_nats / H_cluster, NA_real_)) %>%
        ungroup() %>%
        arrange(cluster, desc(contr_frac)) %>%
        group_by(cluster) %>%
        mutate(rank = row_number()) %>%
        ungroup() %>%
        select(cluster, domain_code, !!value_col, p, contr_nats, contr_frac, H_cluster, rank)
}

# -------------------- static inputs: maps, neighbors --------------------
garbage_root_set <- readr::read_csv(garbage_csv, show_col_types = FALSE) %>%
    transmute(root3 = icd_root3(icd10)) %>% filter(!is.na(root3)) %>% distinct() %>% pull(root3)

gbd_map_df <- readr::read_csv(GBD_ROOT3_CSV, show_col_types = FALSE) |>
    transmute(
        root3  = toupper(substr(.data[[GBD_ROOT3_COL]], 1, 3)),
        gbd_l3 = as.character(.data[[GBD_L3_COL]])
    ) |>
    filter(!is.na(root3), nzchar(root3), !is.na(gbd_l3), nzchar(gbd_l3)) |>
    distinct() |>
    add_count(root3, name = "n_l3") |>
    filter(n_l3 == 1L) |>
    select(root3, gbd_l3)

if (!nrow(gbd_map_df)) stop("GBD root3→L3 map is empty after cleaning.")

counties_sf <- tigris::counties(cb = TRUE, year = 2020) %>% select(GEOID, geometry)
adj <- sf::st_touches(counties_sf)
edge_df <- tibble(from = rep(counties_sf$GEOID, lengths(adj)),
                  to   = counties_sf$GEOID[unlist(adj)]) %>% filter(from < to)
g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = counties_sf$GEOID)
nbrs_list <- setNames(lapply(V(g)$name, function(v) names(neighbors(g, v, mode = "all"))), V(g)$name)

sf::sf_use_s2(FALSE)
counties_proj <- st_transform(counties_sf, crs_proj)
county_centroids <- st_centroid(counties_proj) %>% mutate(GEOID = counties_sf$GEOID) %>% select(GEOID, geometry)
coords <- st_coordinates(county_centroids)
centroid_mat <- coords[, c("X","Y"), drop = FALSE]; rownames(centroid_mat) <- county_centroids$GEOID

# -------------------- Stabilized K (global) pass --------------------
log_msg("Computing stabilized K (global) across ", years_all[1], "–", years_all[2])
con0 <- dbConnect(duckdb::duckdb(), dbdir=":memory:")
on.exit(try(dbDisconnect(con0, shutdown=TRUE), silent=TRUE), add=TRUE)
dbExecute(con0, paste0("PRAGMA threads = ", threads))
duckdb::duckdb_register(con0, "garbage", tibble(root3 = garbage_root_set))
duckdb::duckdb_register(con0, "gbd_root3_map", gbd_map_df)

# Build minimal cert_all for global pass (id, fips, ucod, record_1..20)
MCOD_COLS_ALL <- paste0("record_", 1:20)
dbExecute(con0, glue("
  CREATE TEMP TABLE cert_all AS
  SELECT
    row_number() OVER () AS id,
    LPAD(CAST({`county_var`} AS VARCHAR), 5, '0') AS fips,
    ucod,
    {glue_collapse(MCOD_COLS_ALL, sep = ', ')}
  FROM parquet_scan('{normalizePath(parquet_dir, winslash = '/') }/*.parquet')
  WHERE year BETWEEN {years_all[1]} AND {years_all[2]}
    AND ucod IS NOT NULL
"))

# UCOD (non-garbage)
dbExecute(con0, "
  CREATE TEMP TABLE ucod_valid_all AS
  SELECT id, fips,
         SUBSTR(UPPER(regexp_replace(ucod, '[^A-Za-z0-9]', '')),1,3)  AS root3
  FROM cert_all
  WHERE SUBSTR(UPPER(regexp_replace(ucod,'[^A-Za-z0-9]','')),1,3) NOT IN (SELECT root3 FROM garbage)
")

# K_root3_global
K_root3_global <- dbGetQuery(con0, "SELECT COUNT(DISTINCT root3) AS K FROM ucod_valid_all")$K
maxH_root3_global <- if (K_root3_global > 1) log(K_root3_global) else NA_real_

# K_gbd_global via mapping
dbExecute(con0, "
  CREATE TEMP TABLE ucod_l3_all AS
  SELECT DISTINCT m.gbd_l3
  FROM ucod_valid_all u
  JOIN gbd_root3_map m ON u.root3 = m.root3
")
K_gbd_global <- dbGetQuery(con0, "SELECT COUNT(*) AS K FROM ucod_l3_all")$K
maxH_gbd_global <- if (K_gbd_global > 1) log(K_gbd_global) else NA_real_

# MCOD global non-underlying (drop record_1 if it's UCOD duplicate >90% in your periods; we'll detect below in each period anyway)
dbExecute(con0, "CREATE TEMP TABLE mcod_codes_raw_all (id BIGINT, root3 VARCHAR)")
append_mcod_all <- function(cols_vec) {
    sql_cols <- glue_collapse(cols_vec, sep = ", ")
    dbExecute(con0, glue("
    INSERT INTO mcod_codes_raw_all
    SELECT t.id, t.r3 AS root3
    FROM (
      SELECT id,
             SUBSTR(UPPER(regexp_replace(code, '[^A-Za-z0-9]', '')),1,3) AS r3
      FROM (
        SELECT id, UNNEST(LIST_VALUE({sql_cols})) AS code
        FROM cert_all
      ) u
      WHERE code IS NOT NULL AND code <> '' AND LENGTH(code) >= 3
    ) t
    WHERE t.r3 BETWEEN 'A00' AND 'Z99'
      AND SUBSTR(t.r3, 2, 1) BETWEEN '0' AND '9'
      AND t.r3 NOT IN (SELECT root3 FROM garbage)
  "))
}
chunks_all <- split(MCOD_COLS_ALL, ceiling(seq_along(MCOD_COLS_ALL)/10))
invisible(lapply(chunks_all, append_mcod_all))
dbExecute(con0, "CREATE TEMP TABLE mcod_clean_all AS SELECT DISTINCT id, root3 FROM mcod_codes_raw_all")
# Remove underlying roots from MCOD globally
dbExecute(con0, "
  CREATE TEMP TABLE mcod_nonunderlying_all AS
  SELECT mc.id, mc.root3
  FROM mcod_clean_all mc
  LEFT JOIN ucod_valid_all u
    ON mc.id = u.id AND mc.root3 = u.root3
  WHERE u.root3 IS NULL
")
K_mcod_global <- dbGetQuery(con0, "SELECT COUNT(DISTINCT root3) AS K FROM mcod_nonunderlying_all")$K
maxH_mcod_global <- if (K_mcod_global > 1) log(K_mcod_global) else NA_real_

log_msg("K (global): UCOD root3 = ", K_root3_global,
        " | UCOD GBD L3 = ", K_gbd_global,
        " | MCOD root3 = ", K_mcod_global)

dbDisconnect(con0, shutdown=TRUE)

# -------------------- containers for results --------------------
cluster_membership       <- list()
ucod_detail_root3_list   <- list()
ucod_detail_gbdl3_list   <- list()
mcod_results             <- list()

ucod_contrib_list_root3  <- list()
ucod_contrib_list_gbd    <- list()
mcod_contrib_list_root3  <- list()

# -------------------- main loop over periods --------------------
for (pname in names(periods)) {
    log_msg("Processing ", pname)
    yrs <- range(periods[[pname]])
    
    if (!dir_exists(parquet_dir)) stop("parquet_dir does not exist: ", parquet_dir)
    parq_files <- dir(parquet_dir, pattern = "\\.parquet$", full.names = TRUE)
    if (!length(parq_files)) stop("No .parquet files found in: ", parquet_dir)
    
    con <- dbConnect(duckdb::duckdb(), dbdir=":memory:")
    on.exit(try(dbDisconnect(con, shutdown=TRUE), silent=TRUE), add=TRUE)
    dbExecute(con, paste0("PRAGMA threads = ", threads))
    dbExecute(con, "PRAGMA enable_progress_bar = true")
    dbExecute(con, "PRAGMA progress_bar_time = 1")
    duckdb_register(con, "garbage", tibble(root3 = garbage_root_set))
    duckdb_register(con, "gbd_root3_map", gbd_map_df)
    
    # Minimal certificate table for the period
    MCOD_COLS <- paste0("record_", 1:20)
    dbExecute(con, glue("
    CREATE TEMP TABLE cert AS
    SELECT
      row_number() OVER () AS id,
      LPAD(CAST({`county_var`} AS VARCHAR), 5, '0') AS fips,
      ucod,
      {glue_collapse(MCOD_COLS, sep = ', ')}
    FROM parquet_scan('{normalizePath(parquet_dir, winslash = '/') }/*.parquet')
    WHERE year BETWEEN {yrs[1]} AND {yrs[2]}
      AND ucod IS NOT NULL
  "))
    n_cert <- dbGetQuery(con, "SELECT COUNT(*) n FROM cert")$n
    log_msg("[", pname, "] rows in cert: ", format(n_cert, big.mark=","))
    
    # Detect UCOD duplication in record_1 for this period
    if (AUTO_DETECT_RECORD1_IS_UCOD) {
        chk <- dbGetQuery(con, "
      SELECT COUNT(*) AS n,
             SUM(CASE WHEN SUBSTR(UPPER(regexp_replace(record_1,'[^A-Za-z0-9]','')),1,3)
                           = SUBSTR(UPPER(regexp_replace(ucod,'[^A-Za-z0-9]','')),1,3)
                      THEN 1 ELSE 0 END) AS n_match
      FROM cert
      WHERE record_1 IS NOT NULL AND ucod IS NOT NULL
    ")
        rate <- if (!is.null(chk$n) && chk$n > 0) chk$n_match / chk$n else NA_real_
        if (is.finite(rate) && rate > 0.90) {
            MCOD_COLS <- paste0("record_", 2:20)
            log_msg("Dropping record_1 from MCOD (", sprintf("%.1f%%", 100*rate), " matches UCOD)")
        } else {
            log_msg("Keeping record_1 in MCOD (duplication ~ ", if (is.finite(rate)) sprintf("%.1f%%", 100*rate) else "NA", ")")
        }
    }
    
    # UCOD (non-garbage)
    dbExecute(con, "
    CREATE TEMP TABLE ucod_valid AS
    SELECT id, fips,
           SUBSTR(UPPER(regexp_replace(ucod, '[^A-Za-z0-9]', '')),1,3)  AS root3,
           SUBSTR(UPPER(regexp_replace(ucod, '[^A-Za-z0-9]', '')),1,10) AS icd10_clean
    FROM cert
    WHERE SUBSTR(UPPER(regexp_replace(ucod,'[^A-Za-z0-9]','')),1,3) NOT IN (SELECT root3 FROM garbage)
  ")
    
    # Death totals by FIPS for clustering
    death_tbl <- dbGetQuery(con, "SELECT fips, COUNT(*) AS deaths FROM ucod_valid GROUP BY fips") %>%
        as_tibble() %>% mutate(fips = stringr::str_pad(as.character(fips), 5, pad = "0")) %>%
        filter(fips %in% counties_sf$GEOID)
    
    # ---- Clustering ----
    make_greedy_clusters <- function(death_tbl, nbrs_list, min_deaths) {
        deaths <- setNames(death_tbl$deaths, death_tbl$fips)
        to_assign <- names(deaths)
        clusters <- setNames(rep(NA_character_, length(deaths)), to_assign)
        cid <- 1L
        while (length(to_assign) > 0) {
            this <- to_assign[which.max(deaths[to_assign])]
            if (deaths[this] >= min_deaths) { clusters[this] <- paste0("C", cid); to_assign <- setdiff(to_assign, this); cid <- cp <- cid + 1L; next }
            cluster <- this; total <- deaths[this]; avail <- setdiff(to_assign, this)
            repeat {
                nbrs <- unique(unlist(nbrs_list[cluster], use.names = FALSE)); nbrs <- setdiff(intersect(nbrs, avail), cluster)
                if (!length(nbrs)) break
                best <- nbrs[which.min(abs((total + deaths[nbrs]) - min_deaths))]
                new_tot <- total + deaths[best]
                if (is.na(new_tot)) break
                if (new_tot <= min_deaths * 1.5 || total < min_deaths) { cluster <- c(cluster, best); total <- new_tot; avail <- setdiff(avail, best) } else break
            }
            clusters[cluster] <- paste0("C", cid); to_assign <- setdiff(to_assign, cluster); cid <- cid + 1L
        }
        clusters
    }
    merge_small_clusters <- function(clu, death_tbl, nbrs_list, min_deaths, max_iter = 20L) {
        deaths_vec <- setNames(death_tbl$deaths, death_tbl$fips)
        unlabeled <- names(deaths_vec)[is.na(clu[names(deaths_vec)])]
        if (length(unlabeled)) {
            next_id <- suppressWarnings(max(as.integer(sub("^C","", na.omit(unique(clu)))), na.rm = TRUE))
            if (!is.finite(next_id)) next_id <- 0L
            clu[unlabeled] <- paste0("C", seq(next_id + 1L, next_id + length(unlabeled)))
        }
        iter <- 0L
        repeat {
            iter <- iter + 1L
            sums <- tapply(deaths_vec, clu[names(deaths_vec)], sum, na.rm = TRUE); sums <- sums[is.finite(sums)]
            small <- names(sums)[sums < min_deaths]
            if (!length(small) || iter > max_iter) break
            for (sc in small) {
                members <- names(clu)[clu == sc]; if (!length(members)) next
                neigh_vertices <- unique(unlist(nbrs_list[members], use.names = FALSE))
                neigh_clusters <- setdiff(unique(clu[neigh_vertices]), sc)
                neigh_clusters <- neigh_clusters[neigh_clusters %in% names(sums)]
                target <- if (length(neigh_clusters)) neigh_clusters[which.max(sums[neigh_clusters])] else names(sums)[which.max(sums)]
                if (!is.na(target) && nzchar(target)) clu[members] <- target
            }
        }
        clu
    }
    fix_na_by_nearest <- function(clu, county_centroids) {
        cty <- county_centroids %>% dplyr::select(GEOID, geometry); cty$cluster <- unname(clu[cty$GEOID])
        na_idx <- which(is.na(cty$cluster)); ok_idx <- which(!is.na(cty$cluster)); if (!length(na_idx)) return(clu)
        nn <- st_nearest_feature(cty[na_idx, ], cty[ok_idx, ]); clu[cty$GEOID[na_idx]] <- cty$cluster[ok_idx[nn]]; clu
    }
    merge_small_clusters_by_distance <- function(clu, death_tbl, centroid_mat, min_deaths) {
        deaths_vec <- setNames(death_tbl$deaths, death_tbl$fips); keep <- intersect(names(clu), rownames(centroid_mat))
        sums <- tapply(deaths_vec[keep], clu[keep], sum, na.rm = TRUE); sums <- sums[is.finite(sums)]
        small <- names(sums)[sums < min_deaths]; big <- names(sums)[sums >= min_deaths]
        if (!length(small) || !length(big)) return(clu)
        get_cluster_xy <- function(cid) { memb <- names(clu)[clu == cid]; pts <- centroid_mat[memb, , drop = FALSE]; colMeans(pts) }
        small_xy <- t(vapply(small, get_cluster_xy, numeric(2L))); big_xy <- t(vapply(big, get_cluster_xy, numeric(2L)))
        nearest_idx <- vapply(seq_len(nrow(small_xy)), function(i) { dif <- t(big_xy) - small_xy[i, ]; which.min(colSums(dif * dif)) }, integer(1L))
        targets <- big[nearest_idx]; for (i in seq_along(small)) { sc <- small[i]; tg <- targets[i]; members <- names(clu)[clu == sc]; if (!length(members)) next; clu[members] <- tg }
        clu
    }
    
    clu <- make_greedy_clusters(death_tbl, nbrs_list, min_deaths)
    clu <- merge_small_clusters(clu, death_tbl, nbrs_list, min_deaths)
    clu <- fix_na_by_nearest(clu, county_centroids)
    deaths_vec <- setNames(death_tbl$deaths, death_tbl$fips)
    sums_post <- tapply(deaths_vec, clu[names(deaths_vec)], sum, na.rm = TRUE)
    if (any(is.finite(sums_post) & (sums_post < min_deaths))) {
        clu <- merge_small_clusters_by_distance(clu, death_tbl, centroid_mat, min_deaths)
        clu <- fix_na_by_nearest(clu, county_centroids)
    }
    
    clu_df <- tibble(fips = names(clu), cluster = unname(clu)) %>% left_join(death_tbl, by = "fips")
    cluster_membership[[pname]] <- clu_df %>% mutate(period = pname)
    
    # temp mapping table
    dbExecute(con, "DROP TABLE IF EXISTS map_fips_cluster")
    dbExecute(con, "CREATE TEMP TABLE map_fips_cluster (fips VARCHAR, cluster VARCHAR)")
    DBI::dbAppendTable(con, "map_fips_cluster", clu_df %>% select(fips, cluster) %>% as.data.frame())
    
    # ---- UCOD by root3 ----
    dbExecute(con, "DROP TABLE IF EXISTS ucod_by_fips_root3")
    dbExecute(con, "
    CREATE TEMP TABLE ucod_by_fips_root3 AS
    SELECT fips, root3 AS domain_code, COUNT(*) AS n
    FROM ucod_valid
    GROUP BY 1,2
  ")
    counts_root3 <- dbGetQuery(con, "
    SELECT m.cluster, u.domain_code, SUM(u.n) AS n
    FROM ucod_by_fips_root3 u
    JOIN map_fips_cluster m ON u.fips = m.fips
    WHERE m.cluster IS NOT NULL
    GROUP BY 1,2
  ") %>% as_tibble()
    
    # Stabilized maxH
    maxH_root3 <- maxH_root3_global
    
    res_ucod_root3 <- bind_rows(lapply(split(counts_root3, counts_root3$cluster), function(df) {
        out <- compute_detail_row(df, maxH_root3, m_ref = M_REF); out$cluster <- unique(df$cluster)[1]; out
    })) %>% mutate(period = pname, unit = "cluster", unit_id = cluster)
    ucod_detail_root3_list[[pname]] <- res_ucod_root3 %>% transmute(period, cluster=unit_id, detail_ucod_root3 = detail_phillips_refsize)
    
    # contributions (diagnostic)
    uc_contrib_r3 <- contrib_from_counts(counts_root3, value_col = "n") %>%
        left_join(gbd_map_df %>% rename(domain_code = root3, domain_label = gbd_l3), by = "domain_code")
    uc_contrib_r3$period <- pname
    uc_contrib_r3$family <- "ucod_root3"
    ucod_contrib_list_root3[[pname]] <- uc_contrib_r3
    
    # ---- UCOD by GBD L3 ----
    if (INCLUDE_UCOD_GBDL3) {
        dbExecute(con, "DROP TABLE IF EXISTS ucod_by_fips_gbdl3")
        dbExecute(con, "
      CREATE TEMP TABLE ucod_by_fips_gbdl3 AS
      SELECT u.fips, m.gbd_l3 AS domain_code, COUNT(*) AS n
      FROM ucod_valid u
      JOIN gbd_root3_map m ON u.root3 = m.root3
      GROUP BY 1,2
    ")
        counts_gbd <- dbGetQuery(con, "
      SELECT m.cluster, u.domain_code, SUM(u.n) AS n
      FROM ucod_by_fips_gbdl3 u
      JOIN map_fips_cluster m ON u.fips = m.fips
      WHERE m.cluster IS NOT NULL
      GROUP BY 1,2
    ") %>% as_tibble()
        
        maxH_gbd <- maxH_gbd_global
        res_ucod_gbd <- bind_rows(lapply(split(counts_gbd, counts_gbd$cluster), function(df) {
            out <- compute_detail_row(df, maxH_gbd, m_ref = M_REF); out$cluster <- unique(df$cluster)[1]; out
        })) %>% mutate(period = pname, unit = "cluster", unit_id = cluster)
        ucod_detail_gbdl3_list[[pname]] <- res_ucod_gbd %>% transmute(period, cluster=unit_id, detail_ucod_gbdl3 = detail_phillips_refsize)
        
        uc_contrib_gbd <- contrib_from_counts(counts_gbd, value_col = "n")
        uc_contrib_gbd$domain_label <- uc_contrib_gbd$domain_code
        uc_contrib_gbd$period <- pname
        uc_contrib_gbd$family <- "ucod_gbdl3"
        ucod_contrib_list_gbd[[pname]] <- uc_contrib_gbd
    }
    
    # ---- MCOD (non-underlying) root-3 incidence ----
    dbExecute(con, "DROP TABLE IF EXISTS mcod_codes_raw")
    dbExecute(con, "CREATE TEMP TABLE mcod_codes_raw (id BIGINT, fips VARCHAR, root3 VARCHAR)")
    
    append_mcod_chunk <- function(cols_vec) {
        sql_cols <- glue_collapse(cols_vec, sep = ", ")
        dbExecute(con, glue("
      INSERT INTO mcod_codes_raw
      SELECT t.id, t.fips, t.r3 AS root3
      FROM (
        SELECT id, fips, SUBSTR(UPPER(regexp_replace(code, '[^A-Za-z0-9]', '')),1,3) AS r3
        FROM (
          SELECT id, fips, UNNEST(LIST_VALUE({sql_cols})) AS code
          FROM cert
        ) u
        WHERE code IS NOT NULL AND code <> '' AND LENGTH(code) >= 3
      ) t
      WHERE t.r3 BETWEEN 'A00' AND 'Z99'
        AND SUBSTR(t.r3, 2, 1) BETWEEN '0' AND '9'
        AND t.r3 NOT IN (SELECT root3 FROM garbage)
    "))
    }
    chunks <- split(MCOD_COLS, ceiling(seq_along(MCOD_COLS)/10))
    invisible(lapply(chunks, append_mcod_chunk))
    
    dbExecute(con, "DROP TABLE IF EXISTS mcod_clean")
    dbExecute(con, "CREATE TEMP TABLE mcod_clean AS SELECT DISTINCT id, root3 FROM mcod_codes_raw")
    
    if (isTRUE(EXCLUDE_UCOD_FROM_MCOD)) {
        dbExecute(con, "DROP TABLE IF EXISTS mcod_nonunderlying")
        dbExecute(con, "
      CREATE TEMP TABLE mcod_nonunderlying AS
      SELECT mc.id, mc.root3
      FROM mcod_clean mc
      LEFT JOIN ucod_valid u
        ON mc.id = u.id AND mc.root3 = u.root3
      WHERE u.root3 IS NULL
    ")
    } else {
        dbExecute(con, "CREATE TEMP TABLE mcod_nonunderlying AS SELECT id, root3 FROM mcod_clean")
    }
    
    dbExecute(con, "DROP TABLE IF EXISTS mcod_by_fips")
    dbExecute(con, "
    CREATE TEMP TABLE mcod_by_fips AS
    SELECT r.fips, r.root3, COUNT(DISTINCT r.id) AS incidence
    FROM mcod_codes_raw r
    JOIN mcod_nonunderlying nu ON r.id = nu.id AND r.root3 = nu.root3
    GROUP BY 1,2
  ")
    mcod_incidence <- dbGetQuery(con, "
    SELECT m.cluster, b.root3 AS domain_code, SUM(b.incidence) AS incidence
    FROM mcod_by_fips b
    JOIN map_fips_cluster m ON b.fips = m.fips
    WHERE m.cluster IS NOT NULL
    GROUP BY 1,2
  ") %>% as_tibble()
    
    # Stabilized maxH for MCOD
    maxH_mcod <- maxH_mcod_global
    
    mcod_detail <- mcod_incidence %>%
        group_by(cluster) %>%
        summarise(H = { p <- incidence / sum(incidence); p <- p[p > 0]; -sum(p * log(p)) }, .groups = "drop") %>%
        mutate(detail_mcod_refsize = .clamp100(100 * H / maxH_mcod),
               unit = "cluster", unit_id = cluster, period = pname) %>%
        transmute(unit, unit_id, period, detail_mcod_refsize)
    mcod_results[[pname]] <- mcod_detail
    
    mc_contrib_r3 <- contrib_from_counts(mcod_incidence, value_col = "incidence") %>%
        left_join(gbd_map_df %>% rename(domain_code = root3, domain_label = gbd_l3), by = "domain_code")
    mc_contrib_r3$period <- pname
    mc_contrib_r3$family <- "mcod_root3"
    mcod_contrib_list_root3[[pname]] <- mc_contrib_r3
    
    dbDisconnect(con, shutdown = TRUE); gc()
}

# -------------------- outputs --------------------
make_empty_ucod <- function(which = c("root3","gbd")) {
    which <- match.arg(which)
    if (which == "root3") tibble(period=character(), cluster=character(), detail_ucod_root3=double())
    else tibble(period=character(), cluster=character(), detail_ucod_gbdl3=double())
}
safe_bind_ucod <- function(lst, which) {
    if (length(lst) && any(vapply(lst, nrow, integer(1)) > 0, na.rm = TRUE)) bind_rows(lst) else make_empty_ucod(which)
}

ucod_root3 <- safe_bind_ucod(ucod_detail_root3_list, "root3")
ucod_gbd   <- if (INCLUDE_UCOD_GBDL3) safe_bind_ucod(ucod_detail_gbdl3_list, "gbd") else make_empty_ucod("gbd")
mcod_out   <- bind_rows(mcod_results) %>% transmute(period, cluster = unit_id, detail_mcod = detail_mcod_refsize)

cluster_metrics <- ucod_root3 %>%
    full_join(ucod_gbd, by = c("period","cluster")) %>%
    full_join(mcod_out, by = c("period","cluster")) %>%
    arrange(period, cluster)

cluster_members <- bind_rows(cluster_membership) %>%
    transmute(period, cluster, fips, cluster_deaths = deaths)

# Contribution diagnostics (full)
ucod_contrib_root3 <- bind_rows(ucod_contrib_list_root3)
ucod_contrib_gbd   <- bind_rows(ucod_contrib_list_gbd)
mcod_contrib_root3 <- bind_rows(mcod_contrib_list_root3)

# Also write "top 10" per cluster (handy for QA)
topN <- function(df, n = 10) df %>% group_by(period, cluster, family) %>% slice_min(order_by = rank, n = n, with_ties = FALSE) %>% ungroup()

write_csv(cluster_metrics,             file.path(out_dir, "cluster_metrics.csv"))
write_csv(cluster_members,             file.path(out_dir, "county_cluster_membership.csv"))
write_csv(ucod_contrib_root3,          file.path(out_dir, "detail_contrib_ucod_root3_full.csv"))
write_csv(ucod_contrib_gbd,            file.path(out_dir, "detail_contrib_ucod_gbdl3_full.csv"))
write_csv(mcod_contrib_root3,          file.path(out_dir, "detail_contrib_mcod_root3_full.csv"))
write_csv(topN(ucod_contrib_root3),    file.path(out_dir, "detail_contrib_ucod_root3_top10.csv"))
write_csv(topN(ucod_contrib_gbd),      file.path(out_dir, "detail_contrib_ucod_gbdl3_top10.csv"))
write_csv(topN(mcod_contrib_root3),    file.path(out_dir, "detail_contrib_mcod_root3_top10.csv"))

# Summary
if (nrow(cluster_metrics)) {
    cat("\n— Summary counts by period —\n")
    cluster_metrics %>%
        group_by(period) %>%
        summarise(
            n_clusters = n(),
            nonNA_ucod_root3 = sum(!is.na(detail_ucod_root3)),
            nonNA_ucod_gbd   = sum(!is.na(detail_ucod_gbdl3)),
            nonNA_mcod       = sum(!is.na(detail_mcod)),
            .groups = "drop"
        ) %>% print(n=Inf)
    cat("\nK (global) used: UCOD root3 =", K_root3_global,
        "| UCOD GBD L3 =", K_gbd_global,
        "| MCOD root3 =", K_mcod_global, "\n")
}
cat("\nDone.\n")
