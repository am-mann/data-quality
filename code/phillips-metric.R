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
#
# Updates (Aug 12): added logging, GBD-L3, MCODs
# --------------------------------------------------------------

suppressPackageStartupMessages({
    library(fs); library(readr); library(dplyr); library(stringr)
    library(purrr); library(tidyr); library(tibble)
    library(sf); library(tigris); library(igraph)
    library(here); library(iNEXT)
    library(DBI); library(duckdb); library(glue); library(readxl)
})

# helpers
log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "│", paste0(..., collapse=""), "\n")
log_err <- function(prefix, e) { cat(format(Sys.time(), "%H:%M:%S"), "│ ERROR │", prefix, "│", conditionMessage(e), "\n") }

FIXED_K_GLOBAL <- FALSE        
CLAMP_0_100    <- TRUE      
EXCLUDE_UCOD_FROM_MCOD <- TRUE  
AUTO_DETECT_RECORD1_IS_UCOD <- TRUE

INCLUDE_UCOD_GBDL3        <- TRUE

out_dir     <- here("output"); dir_create(out_dir)
diag_dir    <- file.path(out_dir, "diag"); dir_create(diag_dir)

parquet_dir <- if (dir_exists(here("data_private","mcod"))) here("data_private","mcod") else here("data_private","mcod_sample")

garbage_csv <- here("data_raw","cause-codes","gbd_garbage_codes_without_overdose.csv")

gbd_root3_csv_candidates <- c(
    here("data_raw","cause-codes","lookup_icd10root_to_gbd_lvl3_lvl2.csv"),
    here("data_raw","lookup_icd10root_to_gbd_lvl3_lvl2.csv")
)

gbd_icd_map_xlsx  <- here("data_raw","IHME_GBD_2021_COD_CAUSE_ICD_CODE_MAP.XLSX")
gbd_hier_xlsx     <- here("data_raw","GBD_2021_CAUSE_HIERARCHY.XLSX")

# other params
county_var  <- "county_ihme"
min_deaths  <- 2000L
periods <- list(
    "1999_2006" = 1999:2006,
    "2007_2014" = 2007:2014,
    "2015_2022" = 2015:2022
)
M_REF        <- 2000L
crs_proj     <- 5070
threads      <- max(1L, parallel::detectCores() - 1L)
FAST_MCOD    <- TRUE
FAST_UCOD_AGE<- TRUE

first_existing <- function(paths) {
    p <- paths[file_exists(paths)]
    if (length(p)) p[[1]] else NA_character_
}
icd_root3 <- function(code) {
    x <- toupper(gsub("[^A-Z0-9]", "", as.character(code)))
    m <- stringr::str_match(x, "^([A-Z])([0-9])([0-9A-Z])")
    ifelse(is.na(m[,1]), NA_character_, paste0(m[,2], m[,3], m[,4]))
}
clean_icd <- function(code) toupper(gsub("[^A-Z0-9]", "", as.character(code)))
.clamp100 <- function(x) { if (!CLAMP_0_100) return(x); y <- x; y[x<0] <- 0; y[x>100] <- 100; y }
extract_q1 <- function(est) {
    if (inherits(est, "try-error") || is.null(est)) return(NA_real_)
    oq  <- suppressWarnings(as.numeric(est$Order.q))
    col <- intersect(c("Estimator","qD"), names(est))
    val <- suppressWarnings(as.numeric(est[[col[1]]]))
    out <- val[which(oq == 1)]
    ifelse(length(out) && is.finite(out) && out > 0, out, NA_real_)
}

compute_detail_row <- function(tbl_counts, max_H_period, m_ref = M_REF) {
    ftab <- as.numeric(tbl_counts$n); names(ftab) <- tbl_counts$domain_code
    ftab <- ftab[is.finite(ftab) & ftab > 0]
    n <- sum(ftab); S <- length(ftab); f1 <- sum(ftab == 1); f2 <- sum(ftab == 2)
    if (!is.finite(n) || n < 2 || S < 2 || !is.finite(max_H_period) || max_H_period <= 0) {
        return(tibble(deaths_no_garbage = n, S = S, f1 = f1, f2 = f2,
                      H_raw = NA_real_, detail_phillips_raw = NA_real_, detail_phillips_refsize = NA_real_,
                      D1_mref = NA_real_, D1_eff_used = NA_real_,
                      m_ref = m_ref, eff_m_used = NA_integer_, used_extrapolation = FALSE, fallback_reason = "degenerate"))
    }
    p <- ftab / n
    H_raw <- -sum(p * log(p))
    detail_raw <- .clamp100(100 * H_raw / max_H_period)
    m_use <- as.integer(min(m_ref, n))
    est   <- try(iNEXT::estimateD(ftab, datatype = "abundance", base = "size", level = m_use, conf = FALSE), silent = TRUE)
    D1_use <- extract_q1(est)
    if (is.na(D1_use) && is.finite(H_raw)) D1_use <- exp(H_raw)
    tibble(
        deaths_no_garbage = n, S = S, f1 = f1, f2 = f2, H_raw = H_raw,
        detail_phillips_raw = detail_raw,
        detail_phillips_refsize = if (is.finite(D1_use) && D1_use > 0) .clamp100(100 * log(D1_use) / max_H_period) else NA_real_,
        D1_mref = if (m_ref <= n) D1_use else NA_real_, D1_eff_used = D1_use,
        m_ref = as.integer(m_ref), eff_m_used = m_use, used_extrapolation = FALSE, fallback_reason = "ok"
    )
}

build_gbd_l3_map_excel <- function(path_icd_map, path_hierarchy, verbose = TRUE) {
    if (!file.exists(path_icd_map) || !file.exists(path_hierarchy)) return(NULL)
    if (!requireNamespace("readxl", quietly = TRUE)) return(NULL)
    
    pick_col <- function(nms, patterns) {
        for (pat in patterns) {
            hit <- nms[grepl(pat, nms, ignore.case = TRUE)]
            if (length(hit)) return(hit[1])
        }
        NA_character_
    }
    
    h <- try(readxl::read_excel(path_hierarchy), silent = TRUE)
    if (inherits(h, "try-error") || !nrow(h)) return(NULL)
    nms_h <- names(h)
    col_id     <- pick_col(nms_h, c("^cause_?id$","^gbd_?cause_?id$","^id$"))
    col_parent <- pick_col(nms_h, c("^parent_?id$","^parent$"))
    col_level  <- pick_col(nms_h, c("^level$","^gbd_?level$"))
    col_name   <- pick_col(nms_h, c("cause_?name","gbd_?cause_?name"))
    if (any(is.na(c(col_id, col_parent, col_level, col_name)))) return(NULL)
    
    h <- h |>
        transmute(
            cause_id   = suppressWarnings(as.integer(.data[[col_id]])),
            parent_id  = suppressWarnings(as.integer(.data[[col_parent]])),
            level      = suppressWarnings(as.integer(.data[[col_level]])),
            cause_name = as.character(.data[[col_name]])
        ) |>
        filter(is.finite(cause_id), is.finite(level)) |>
        distinct()
    
    sheets <- try(readxl::excel_sheets(path_icd_map), silent = TRUE)
    if (inherits(sheets, "try-error")) return(NULL)
    
    expand_icd_tokens <- function(x) {
        if (is.na(x) || !nzchar(x)) return(character(0))
        x <- gsub("\u2013|\u2014", "-", x)
        toks <- unlist(strsplit(x, "[,;/\\s]+")); toks <- toks[nzchar(toks)]
        out <- character(0)
        for (t in toks) {
            if (grepl("-", t)) {
                parts <- strsplit(t, "-", fixed = TRUE)[[1]]
                if (length(parts) != 2) next
                a <- toupper(gsub("[^A-Z0-9.]", "", parts[1]))
                b <- toupper(gsub("[^A-Z0-9.]", "", parts[2]))
                if (!nzchar(a) || !nzchar(b)) next
                la <- substr(a, 1, 1); lb <- substr(b, 1, 1)
                if (la != lb) next
                aN <- gsub("\\.", "", substr(a, 2, nchar(a)))
                bN <- gsub("\\.", "", substr(b, 2, nchar(b)))
                if (!grepl("^[0-9]+$", aN) || !grepl("^[0-9]+$", bN)) next
                width <- max(nchar(aN), nchar(bN))
                ai <- suppressWarnings(as.integer(aN)); bi <- suppressWarnings(as.integer(bN))
                if (!is.finite(ai) || !is.finite(bi) || ai > bi || (bi - ai) > 2000) next
                out <- c(out, paste0(la, stringr::str_pad(seq(ai, bi), width = width, pad = "0")))
            } else {
                out <- c(out, toupper(gsub("[^A-Z0-9]", "", t)))
            }
        }
        unique(out[nzchar(out)])
    }
    
    pieces <- lapply(sheets, function(snm) {
        x <- try(readxl::read_excel(path_icd_map, sheet = snm), silent = TRUE)
        if (inherits(x, "try-error") || !nrow(x)) return(NULL)
        nmx <- names(x)
        col_icd <- pick_col(nmx, c("^icd.*10", "^icd10", "icd.?codes?", "icd.*code"))
        col_id2 <- pick_col(nmx, c("^cause_?id$", "^gbd_?cause_?id$", "^id$"))
        if (is.na(col_icd) || is.na(col_id2)) return(NULL)
        tibble(
            cause_id = suppressWarnings(as.integer(x[[col_id2]])),
            tokens   = as.character(x[[col_icd]])
        ) |>
            filter(is.finite(cause_id), !is.na(tokens), nzchar(tokens)) |>
            mutate(icd10_clean = lapply(tokens, expand_icd_tokens)) |>
            unnest(icd10_clean) |>
            mutate(icd10_clean = toupper(gsub("[^A-Z0-9]", "", icd10_clean))) |>
            transmute(cause_id, icd10_clean) |>
            distinct()
    })
    icdmap <- bind_rows(pieces)
    if (!nrow(icdmap)) return(NULL)
    
    parent_map <- h |> select(cause_id, parent_id, level)
    ascend_to_l3 <- function(ids) {
        ids <- as.integer(ids); res <- ids
        for (k in 1:8) {
            lev <- parent_map$level[match(res, parent_map$cause_id)]
            need <- which(is.finite(lev) & lev > 3)
            if (!length(need)) break
            res[need] <- parent_map$parent_id[match(res[need], parent_map$cause_id)]
        }
        lev <- parent_map$level[match(res, parent_map$cause_id)]
        low <- which(is.finite(lev) & lev < 3)
        for (i in low) {
            cur <- res[i]
            for (k in 1:8) {
                cur <- parent_map$parent_id[match(cur, parent_map$cause_id)]
                if (is.na(cur)) break
                if (parent_map$level[match(cur, parent_map$cause_id)] == 3) { res[i] <- cur; break }
            }
        }
        res
    }
    
    icdmap$l3_id <- ascend_to_l3(icdmap$cause_id)
    l3_names <- h |> filter(level == 3) |> select(cause_id, gbd_l3 = cause_name)
    
    out <- icdmap |>
        inner_join(l3_names, by = c("l3_id" = "cause_id")) |>
        transmute(
            icd10_clean_10 = substr(icd10_clean, 1, 10),
            root3          = substr(icd10_clean, 1, 3),
            gbd_l3         = as.character(gbd_l3)
        ) |>
        distinct()
    
    out
}

load_gbd_l3_map <- function(csv_candidates, xlsx_map, xlsx_hier) {
    # Try CSV(s)
    csv_path <- first_existing(csv_candidates)
    if (is.character(csv_path) && !is.na(csv_path)) {
        df <- suppressMessages(readr::read_csv(csv_path, show_col_types = FALSE)) %>%
            transmute(
                root3  = toupper(stringr::str_sub(coalesce(icd10_root, icd10, root, ROOT, `icd10-root`, `ICD10_ROOT`), 1, 3)),
                gbd_l3 = as.character(coalesce(gbd_cause_lvl3, gbd_l3, L3, lvl3, cause, Cause))
            ) %>%
            filter(!is.na(root3), nzchar(root3), !is.na(gbd_l3), nzchar(gbd_l3)) %>%
            distinct() %>%
            add_count(root3, name = "n_l3") %>%
            filter(n_l3 == 1L) %>%
            select(root3, gbd_l3) %>%
            distinct()
        if (nrow(df)) {
            attr(df, "mode") <- "root3"

            return(df)
        }
    }
    ex <- build_gbd_l3_map_excel(xlsx_map, xlsx_hier, verbose = TRUE)
    if (!is.null(ex) && nrow(ex)) {
        attr(ex, "mode") <- "excel"
        return(ex)
    }
    NULL
}

paths_to_check <- c(
    garbage_csv = garbage_csv,
    parquet_dir = parquet_dir,
    gbd_root3_csv_candidates[1],
    gbd_root3_csv_candidates[2],
    gbd_icd_map_xlsx = gbd_icd_map_xlsx,
    gbd_hier_xlsx = gbd_hier_xlsx)


stopifnot(file.exists(garbage_csv))
garbage_root_set <- readr::read_csv(garbage_csv, show_col_types = FALSE) %>%
    transmute(root3 = icd_root3(icd10)) %>% filter(!is.na(root3)) %>% distinct() %>% pull(root3)

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

GBD_ROOT3_CSV <- here("data_raw","lookup_icd10root_to_gbd_lvl3_lvl2.csv")


GBD_ROOT3_COL <- "icd10_root"       # e.g., "icd10_root" or "root3"
GBD_L3_COL    <- "gbd_cause_lvl3"   # e.g., "gbd_cause_lvl3" or "gbd_l3"

stopifnot(file.exists(GBD_ROOT3_CSV))

gbd_map_df <- readr::read_csv(GBD_ROOT3_CSV, show_col_types = FALSE) |>
    dplyr::transmute(
        root3  = toupper(substr(.data[[GBD_ROOT3_COL]], 1, 3)),
        gbd_l3 = as.character(.data[[GBD_L3_COL]])
    ) |>
    dplyr::filter(!is.na(root3), nzchar(root3), !is.na(gbd_l3), nzchar(gbd_l3)) |>
    dplyr::distinct() |>
    dplyr::add_count(root3, name = "n_l3") |>
    dplyr::filter(n_l3 == 1L) |>
    dplyr::select(root3, gbd_l3)

if (!nrow(gbd_map_df)) {
    stop(sprintf(
        "GBD root3 map is empty after cleaning. Check columns '%s' and '%s' in %s.",
        GBD_ROOT3_COL, GBD_L3_COL, GBD_ROOT3_CSV
    ))
}

attr(gbd_map_df, "mode") <- "root3"
GBD_JOIN_MODE <- "root3"



age_case_sql <- "
CASE
  WHEN TRY_CAST(age AS DOUBLE) IS NULL THEN NULL
  WHEN CAST(age AS DOUBLE) < 45 THEN '0-44'
  WHEN CAST(age AS DOUBLE) < 65 THEN '45-64'
  WHEN CAST(age AS DOUBLE) < 80 THEN '65-79'
  ELSE '80+'
END
"
years_all <- range(unlist(periods))
conw <- DBI::dbConnect(duckdb::duckdb(), dbdir=":memory:")
DBI::dbExecute(conw, paste0("PRAGMA threads = ", threads))
duckdb::duckdb_register(conw, "garbage", tibble::tibble(root3 = garbage_root_set))
age_w_tbl <- DBI::dbGetQuery(conw, glue("
  SELECT age_bucket, COUNT(*) AS deaths
  FROM (SELECT {age_case_sql} AS age_bucket
        FROM parquet_scan('{normalizePath(parquet_dir, winslash = '/') }/*.parquet')
        WHERE year BETWEEN {years_all[1]} AND {years_all[2]}) t
  WHERE age_bucket IS NOT NULL
  GROUP BY 1
  ORDER BY CASE age_bucket WHEN '0-44' THEN 1 WHEN '45-64' THEN 2 WHEN '65-79' THEN 3 ELSE 4 END
")) %>% as_tibble()
age_weights <- age_w_tbl %>% mutate(w = deaths / sum(deaths)) %>% transmute(age_bucket, w)
write_csv(age_weights, file.path(out_dir, "age_weights_skew4_global.csv"))
DBI::dbDisconnect(conw, shutdown=TRUE)


cluster_membership       <- list()
ucod_detail_root3_list   <- list()
ucod_detail_gbdl3_list   <- list()
phillips_age_results     <- list()
mcod_results             <- list()
mcod_age_results         <- list()

# main
for (pname in names(periods)) {
    log_msg("Processing ", pname)
    yrs <- range(periods[[pname]])
    
    if (!dir.exists(parquet_dir)) stop("parquet_dir does not exist: ", parquet_dir)
    parq_files <- list.files(parquet_dir, pattern = "\\.parquet$", full.names = TRUE)
    if (length(parq_files) == 0L) stop("No .parquet files found in: ", parquet_dir)
    
    con <- dbConnect(duckdb::duckdb(), dbdir=":memory:")
    on.exit(try(dbDisconnect(con, shutdown=TRUE), silent=TRUE), add=TRUE)
    dbExecute(con, paste0("PRAGMA threads = ", threads))
    dbExecute(con, "PRAGMA enable_progress_bar = true")
    dbExecute(con, "PRAGMA progress_bar_time = 1")
    duckdb_register(con, "garbage", tibble(root3 = garbage_root_set))
    if (!is.null(gbd_map_df) && nrow(gbd_map_df)) duckdb_register(con, "gbd_map", gbd_map_df)

    sql_cert <- glue("
    CREATE TEMP TABLE cert AS
    SELECT
      row_number() OVER () AS id,
      LPAD(CAST({`county_var`} AS VARCHAR), 5, '0') AS fips,
      ucod, age,
      {glue_collapse(paste0('record_', 1:20), sep = ', ')}
    FROM parquet_scan('{normalizePath(parquet_dir, winslash = '/') }/*.parquet')
    WHERE year BETWEEN {yrs[1]} AND {yrs[2]}
      AND ucod IS NOT NULL
  ")
    dbExecute(con, sql_cert)
    n_cert <- dbGetQuery(con, "SELECT COUNT(*) n FROM cert")$n

    MCOD_COLS <- paste0("record_", 1:20)
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
        msg  <- if (is.finite(rate)) sprintf("%.1f%%", 100*rate) else "NA"
        if (is.finite(rate) && rate > 0.90) {
            MCOD_COLS <- paste0("record_", 2:20)
        } else {
        }
    }
    
    dbExecute(con, "
    CREATE TEMP TABLE ucod_valid AS
    SELECT id, fips,
           SUBSTR(UPPER(regexp_replace(ucod, '[^A-Za-z0-9]', '')),1,3)  AS root3,
           SUBSTR(UPPER(regexp_replace(ucod, '[^A-Za-z0-9]', '')),1,10) AS icd10_clean
    FROM cert
    WHERE ucod IS NOT NULL
      AND SUBSTR(UPPER(regexp_replace(ucod,'[^A-Za-z0-9]','')),1,3) NOT IN (SELECT root3 FROM garbage)
  ")
    n_uv <- dbGetQuery(con, "SELECT COUNT(*) n FROM ucod_valid")$n
    k_uv <- dbGetQuery(con, "SELECT COUNT(DISTINCT root3) k FROM ucod_valid")$k

    death_tbl <- dbGetQuery(con, "SELECT fips, COUNT(*) AS deaths FROM ucod_valid GROUP BY fips") %>%
        as_tibble() %>% mutate(fips = stringr::str_pad(as.character(fips), 5, pad = "0")) %>%
        filter(fips %in% counties_sf$GEOID)
    
    # clustering
    make_greedy_clusters <- function(death_tbl, nbrs_list, min_deaths) {
        deaths <- setNames(death_tbl$deaths, death_tbl$fips); to_assign <- names(deaths)
        clusters <- setNames(rep(NA_character_, length(deaths)), to_assign); cid <- 1L
        while (length(to_assign) > 0) {
            this <- to_assign[which.max(deaths[to_assign])]
            if (deaths[this] >= min_deaths) { clusters[this] <- paste0("C", cid); to_assign <- setdiff(to_assign, this); cid <- cid + 1L; next }
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
            sums <- tapply(deaths_vec, clu[names(deaths_vec)], sum, na.rm = TRUE)
            sums <- sums[is.finite(sums)]
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
    
    dbExecute(con, "DROP TABLE IF EXISTS map_fips_cluster")
    dbExecute(con, "CREATE TEMP TABLE map_fips_cluster (fips VARCHAR, cluster VARCHAR)")
    DBI::dbAppendTable(con, "map_fips_cluster", clu_df %>% select(fips, cluster) %>% as.data.frame())
    
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
    K_root3_period <- dbGetQuery(con, "SELECT COUNT(DISTINCT root3) AS K FROM ucod_valid")$K
    maxH_root3 <- if (FIXED_K_GLOBAL) NA_real_ else if (K_root3_period>1) log(K_root3_period) else NA_real_
    res_ucod_root3 <- bind_rows(lapply(split(counts_root3, counts_root3$cluster), function(df) {
        out <- compute_detail_row(df, maxH_root3, m_ref = M_REF); out$cluster <- unique(df$cluster)[1]; out
    })) %>% mutate(period = pname, unit = "cluster", unit_id = cluster)
    ucod_detail_root3_list[[pname]] <- res_ucod_root3 %>% transmute(period, cluster=unit_id, detail_ucod_root3 = detail_phillips_refsize)
    
    if (!is.null(gbd_map_df) && nrow(gbd_map_df)) {
        dbExecute(con, "DROP TABLE IF EXISTS ucod_by_fips_gbdl3")
        
        if (GBD_JOIN_MODE == "excel") {
            # Try exact code join first
            dbExecute(con, "
        CREATE TEMP TABLE ucod_by_fips_gbdl3 AS
        SELECT u.fips, gm.gbd_l3 AS domain_code, COUNT(*) AS n
        FROM ucod_valid u
        JOIN (SELECT DISTINCT icd10_clean_10, gbd_l3 FROM gbd_map) gm
          ON u.icd10_clean = gm.icd10_clean_10
        GROUP BY 1,2
      ")
            n_exact <- dbGetQuery(con, "SELECT COUNT(*) n FROM ucod_by_fips_gbdl3")$n
            if (is.na(n_exact) || n_exact == 0L) {
                # Build a UNIQUE root3→L3 mapping from the excel map (drop ambiguous)
                gbd_root3_unique <- gbd_map_df %>%
                    distinct(root3, gbd_l3) %>%
                    add_count(root3, name = "n_l3") %>%
                    filter(n_l3 == 1L) %>% select(root3, gbd_l3)
                if (nrow(gbd_root3_unique)) {
                    duckdb::duckdb_register(con, "gbd_root3_unique", gbd_root3_unique)
                    dbExecute(con, "DELETE FROM ucod_by_fips_gbdl3") # clear
                    dbExecute(con, "
            INSERT INTO ucod_by_fips_gbdl3
            SELECT u.fips, g.gbd_l3 AS domain_code, COUNT(*) AS n
            FROM ucod_valid u
            JOIN gbd_root3_unique g ON u.root3 = g.root3
            GROUP BY 1,2
          ")
                }
            }
        } else {
            # root-3 CSV mapping (already unique)
            duckdb::duckdb_register(con, "gbd_root3_map", gbd_map_df)
            dbExecute(con, "
        CREATE TEMP TABLE ucod_by_fips_gbdl3 AS
        SELECT u.fips, m.gbd_l3 AS domain_code, COUNT(*) AS n
        FROM ucod_valid u
        JOIN gbd_root3_map m ON u.root3 = m.root3
        GROUP BY 1,2
      ")
        }
        
        k_rows <- dbGetQuery(con, "SELECT COUNT(*) n FROM ucod_by_fips_gbdl3")$n
        if (is.finite(k_rows) && k_rows > 0) {
            counts_gbd <- dbGetQuery(con, "
        SELECT m.cluster, u.domain_code, SUM(u.n) AS n
        FROM ucod_by_fips_gbdl3 u
        JOIN map_fips_cluster m ON u.fips = m.fips
        WHERE m.cluster IS NOT NULL
        GROUP BY 1,2
      ") %>% as_tibble()
            
            K_gbd_period <- dbGetQuery(con, "SELECT COUNT(DISTINCT domain_code) AS K FROM ucod_by_fips_gbdl3")$K
            maxH_gbd <- if (K_gbd_period>1) log(K_gbd_period) else NA_real_
            res_ucod_gbd <- bind_rows(lapply(split(counts_gbd, counts_gbd$cluster), function(df) {
                out <- compute_detail_row(df, maxH_gbd, m_ref = M_REF); out$cluster <- unique(df$cluster)[1]; out
            })) %>% mutate(period = pname, unit = "cluster", unit_id = cluster)
            ucod_detail_gbdl3_list[[pname]] <- res_ucod_gbd %>% transmute(period, cluster=unit_id, detail_ucod_gbdl3 = detail_phillips_refsize)
            
        } else {
        }
    } else {
    }
    
    counts_by_age <- dbGetQuery(con, glue("
    WITH u AS (
      SELECT c.id, c.fips, {age_case_sql} AS age_bucket, uv.root3
      FROM cert c JOIN ucod_valid uv ON c.id = uv.id
    )
    SELECT m.cluster, u.age_bucket, u.root3 AS domain_code, COUNT(*) AS n
    FROM u JOIN map_fips_cluster m ON u.fips = m.fips
    WHERE m.cluster IS NOT NULL AND u.age_bucket IS NOT NULL
    GROUP BY 1,2,3
  ")) %>% as_tibble()
    maxH_age_map <- counts_by_age %>%
        group_by(age_bucket) %>% summarise(K = n_distinct(domain_code), .groups = "drop") %>%
        mutate(maxH = ifelse(K > 1, log(K), NA_real_)) %>% { setNames(.$maxH, .$age_bucket) }
    if (FAST_UCOD_AGE) {
        res_phillips_age <- counts_by_age %>%
            group_by(cluster, age_bucket) %>%
            summarise(H = { p <- n / sum(n); p <- p[p > 0]; -sum(p * log(p)) }, .groups = "drop") %>%
            mutate(maxH = as.numeric(maxH_age_map[as.character(age_bucket)]),
                   detail_phillips_refsize = ifelse(is.finite(maxH) & maxH > 0, .clamp100(100 * H / maxH), NA_real_),
                   unit = "cluster", unit_id = cluster, period = pname) %>%
            transmute(unit, unit_id, period, age_bucket, detail_phillips_refsize)
        phillips_age_results[[pname]] <- res_phillips_age
    }
    
    dbExecute(con, "DROP TABLE IF EXISTS mcod_codes_raw")
    dbExecute(con, "CREATE TEMP TABLE mcod_codes_raw (id BIGINT, fips VARCHAR, root3 VARCHAR)")
    append_mcod_chunk <- function(cols_vec) {
        sql_cols <- glue_collapse(cols_vec, sep = ", ")
        sql <- glue("
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
    ")
        dbExecute(con, sql)
    }
    MCOD_COLS <- unname(MCOD_COLS)
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
    
    K_mcod_period <- dbGetQuery(con, "SELECT COUNT(DISTINCT root3) AS K FROM mcod_nonunderlying")$K
    maxH_mcod <- if (K_mcod_period>1) log(K_mcod_period) else NA_real_
    if (isTRUE(FAST_MCOD) && is.finite(maxH_mcod) && maxH_mcod > 0) {
        mcod_detail <- mcod_incidence %>%
            group_by(cluster) %>%
            summarise(H = { p <- incidence / sum(incidence); p <- p[p > 0]; -sum(p * log(p)) }, .groups = "drop") %>%
            mutate(detail_mcod_refsize = .clamp100(100 * H / maxH_mcod),
                   unit = "cluster", unit_id = cluster, period = pname) %>%
            transmute(unit, unit_id, period, detail_mcod_refsize)
        mcod_results[[pname]] <- mcod_detail
    }
    
    # MCOD by age
    mcod_incidence_age <- dbGetQuery(con, glue("
    WITH ages AS ( SELECT id, fips, {age_case_sql} AS age_bucket FROM cert )
    SELECT m.cluster, a.age_bucket, nu.root3 AS domain_code, COUNT(DISTINCT a.id) AS incidence
    FROM mcod_nonunderlying nu
    JOIN ages a ON nu.id = a.id
    JOIN map_fips_cluster m ON a.fips = m.fips
    WHERE m.cluster IS NOT NULL AND a.age_bucket IS NOT NULL
    GROUP BY 1,2,3
  ")) %>% as_tibble()
    maxH_mcod_age_map <- mcod_incidence_age %>%
        group_by(age_bucket) %>% summarise(K = n_distinct(domain_code), .groups = "drop") %>%
        mutate(maxH = ifelse(K > 1, log(K), NA_real_)) %>% { setNames(.$maxH, .$age_bucket) }
    mcod_detail_age <- mcod_incidence_age %>%
        group_by(cluster, age_bucket) %>%
        summarise(H = { p <- incidence / sum(incidence); p <- p[p > 0]; -sum(p * log(p)) }, .groups = "drop") %>%
        mutate(maxH = as.numeric(maxH_mcod_age_map[as.character(age_bucket)]),
               detail_mcod_refsize = ifelse(is.finite(maxH) & maxH > 0, .clamp100(100 * H / maxH), NA_real_),
               unit = "cluster", unit_id = cluster, period = pname) %>%
        transmute(unit, unit_id, period, age_bucket, detail_mcod_refsize)
    mcod_age_results[[pname]] <- mcod_detail_age
    
    dbDisconnect(con, shutdown = TRUE); gc()
}

# outputs
make_empty_ucod <- function(which = c("root3","gbd")) {
    which <- match.arg(which)
    if (which == "root3") tibble(period=character(), cluster=character(), detail_ucod_root3=double())
    else tibble(period=character(), cluster=character(), detail_ucod_gbdl3=double())
}
safe_bind_ucod <- function(lst, which) {
    if (length(lst) && any(vapply(lst, nrow, integer(1)) > 0, na.rm = TRUE)) bind_rows(lst) else make_empty_ucod(which)
}

ucod_root3 <- safe_bind_ucod(ucod_detail_root3_list, "root3")
ucod_gbd   <- safe_bind_ucod(ucod_detail_gbdl3_list, "gbd")

ucod_multi <- ucod_root3 %>%
    full_join(ucod_gbd, by = c("period","cluster")) %>%
    arrange(period, cluster)

phillips_age_out <- bind_rows(phillips_age_results)
mcod_age_out     <- bind_rows(mcod_age_results)
safe_wmean <- function(x, w) { i <- is.finite(x) & is.finite(w); if (!any(i)) return(NA_real_); sum(x[i]*w[i]) / sum(w[i]) }

ucod_age_std <- phillips_age_out %>%
    inner_join(age_weights, by = "age_bucket") %>%
    group_by(unit, unit_id, period) %>%
    summarise(detail_ucod_age_std = safe_wmean(detail_phillips_refsize, w),
              n_buckets = sum(is.finite(detail_phillips_refsize)), .groups = "drop") %>%
    transmute(period, cluster = unit_id, detail_ucod_age_std = .clamp100(detail_ucod_age_std))

mcod_age_std <- mcod_age_out %>%
    inner_join(age_weights, by = "age_bucket") %>%
    group_by(unit, unit_id, period) %>%
    summarise(detail_mcod_age_std = safe_wmean(detail_mcod_refsize, w),
              n_buckets = sum(is.finite(detail_mcod_refsize)), .groups = "drop") %>%
    transmute(period, cluster = unit_id, detail_mcod_age_std = .clamp100(detail_mcod_age_std))

mcod_out <- bind_rows(mcod_results) %>% transmute(period, cluster = unit_id, detail_mcod = detail_mcod_refsize)

cluster_metrics <- ucod_multi %>%
    full_join(mcod_out,     by = c("period","cluster")) %>%
    full_join(ucod_age_std, by = c("period","cluster")) %>%
    full_join(mcod_age_std, by = c("period","cluster")) %>%
    mutate(detail_age_std_COMPOSITE = .clamp100((detail_ucod_age_std + detail_mcod_age_std)/2)) %>%
    arrange(period, cluster)

cluster_members <- bind_rows(cluster_membership) %>%
    transmute(period, cluster, fips, cluster_deaths = deaths)

# write outputs
write_csv(cluster_metrics, file.path(out_dir, "cluster_metrics.csv"))
write_csv(cluster_members, file.path(out_dir, "county_cluster_membership.csv"))

# ---- Save age-bucket detail (long) ----
ucod_age_by_bucket <- phillips_age_out %>%
    transmute(period, cluster = unit_id, age_bucket,
              detail_ucod_age = .clamp100(detail_phillips_refsize))
mcod_age_by_bucket <- mcod_age_out %>%
    transmute(period, cluster = unit_id, age_bucket,
              detail_mcod_age = .clamp100(detail_mcod_refsize))

write_csv(ucod_age_by_bucket, file.path(out_dir, "ucod_detail_by_age_bucket.csv"))
write_csv(mcod_age_by_bucket, file.path(out_dir, "mcod_detail_by_age_bucket.csv"))

ucod_age_wide <- ucod_age_by_bucket %>%
    tidyr::pivot_wider(names_from = age_bucket, values_from = detail_ucod_age,
                       names_prefix = "ucod_age_")
mcod_age_wide <- mcod_age_by_bucket %>%
    tidyr::pivot_wider(names_from = age_bucket, values_from = detail_mcod_age,
                       names_prefix = "mcod_age_")

write_csv(ucod_age_wide, file.path(out_dir, "ucod_detail_by_age_bucket_wide.csv"))
write_csv(mcod_age_wide, file.path(out_dir, "mcod_detail_by_age_bucket_wide.csv"))


if (nrow(cluster_metrics)) {
    cat("\n— Summary counts by period —\n")
    cluster_metrics %>%
        group_by(period) %>%
        summarise(
            n_clusters = n(),
            nonNA_ucod_root3 = sum(!is.na(detail_ucod_root3)),
            nonNA_ucod_gbd   = sum(!is.na(detail_ucod_gbdl3)),
            nonNA_mcod       = sum(!is.na(detail_mcod)),
            nonNA_ucod_age   = sum(!is.na(detail_ucod_age_std)),
            nonNA_mcod_age   = sum(!is.na(detail_mcod_age_std)),
            .groups = "drop"
        ) %>% print(n=Inf)
}

cat("\nDone.\n")

