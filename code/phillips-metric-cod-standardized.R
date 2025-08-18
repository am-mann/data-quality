# --------------------------------------------------------------
# ANACONDA detail metric
# Date: August 18, 2025
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
# CHANGES: It is now standardized by ucr39 cause of death categories.
#
# --------------------------------------------------------------

suppressPackageStartupMessages({
    library(fs); library(readr); library(dplyr); library(stringr)
    library(purrr); library(tidyr); library(tibble)
    library(sf); library(tigris); library(igraph)
    library(here); library(DBI); library(duckdb); library(glue)
})

# -------------------- helpers & params --------------------
log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "│", paste0(..., collapse=""), "\n")
CLAMP_0_100 <- TRUE
.clamp100 <- function(x){ if(!CLAMP_0_100) return(x); y <- x; y[x<0]<-0; y[x>100]<-100; y }

# Paths
out_dir   <- here("output"); dir_create(out_dir)
parquet_dir <- if (dir_exists(here("data_private","mcod"))) here("data_private","mcod") else here("data_private","mcod_sample")

garbage_csv <- here("data_raw","cause-codes","gbd_garbage_codes_without_overdose.csv")
stopifnot(file.exists(garbage_csv))

# Data columns
county_var <- "county_ihme"
UCR39_COL  <- "ucr39"   # <-- already in your files

# Periods and clustering
min_deaths <- 2000L
periods <- list(
    "1999_2005" = 1999:2005,
    "2006_2012" = 2006:2012,
    "2013_2019" = 2013:2019,
    "2020_2022" = 2020:2022
)
years_all <- range(unlist(periods))
crs_proj  <- 5070
threads   <- max(1L, parallel::detectCores() - 1L)

# Cleaners
icd_root3 <- function(code){
    x <- toupper(gsub("[^A-Z0-9]","",as.character(code)))
    m <- stringr::str_match(x, "^([A-Z])([0-9])([0-9A-Z])")
    ifelse(is.na(m[,1]), NA_character_, paste0(m[,2],m[,3],m[,4]))
}
icd4_clean <- function(code){
    x <- toupper(gsub("[^A-Z0-9]","",as.character(code)))
    substr(x,1,4)
}

# -------------------- reference data --------------------
garbage_root_set <- readr::read_csv(garbage_csv, show_col_types = FALSE) %>%
    transmute(root3 = icd_root3(icd10)) %>% filter(!is.na(root3)) %>% distinct() %>% pull(root3)

# Geography & adjacency (for clustering)
sf::sf_use_s2(FALSE)
counties_sf <- tigris::counties(cb = TRUE, year = 2020) %>% select(GEOID, geometry)
adj <- sf::st_touches(counties_sf)
edge_df <- tibble(from = rep(counties_sf$GEOID, lengths(adj)),
                  to   = counties_sf$GEOID[unlist(adj)]) %>% filter(from < to)
g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = counties_sf$GEOID)
nbrs_list <- setNames(lapply(V(g)$name, function(v) names(neighbors(g, v, mode="all"))), V(g)$name)

counties_proj <- st_transform(counties_sf, crs_proj)
county_centroids <- st_centroid(counties_proj) %>% mutate(GEOID = counties_sf$GEOID) %>% select(GEOID, geometry)
coords <- st_coordinates(county_centroids)
centroid_mat <- coords[, c("X","Y"), drop=FALSE]; rownames(centroid_mat) <- county_centroids$GEOID

# -------------------- clustering helpers --------------------
make_greedy_clusters <- function(death_tbl, nbrs_list, min_deaths){
    deaths <- setNames(death_tbl$deaths, death_tbl$fips)
    to_assign <- names(deaths)
    clusters <- setNames(rep(NA_character_, length(deaths)), to_assign)
    cid <- 1L
    while(length(to_assign) > 0){
        this <- to_assign[which.max(deaths[to_assign])]
        if (deaths[this] >= min_deaths){
            clusters[this] <- paste0("C", cid); to_assign <- setdiff(to_assign, this); cid <- cid + 1L; next
        }
        cluster <- this; total <- deaths[this]; avail <- setdiff(to_assign, this)
        repeat{
            nbrs <- unique(unlist(nbrs_list[cluster], use.names = FALSE))
            nbrs <- setdiff(intersect(nbrs, avail), cluster)
            if (!length(nbrs)) break
            best <- nbrs[which.min(abs((total + deaths[nbrs]) - min_deaths))]
            new_tot <- total + deaths[best]
            if (is.na(new_tot)) break
            if (new_tot <= min_deaths * 1.5 || total < min_deaths){
                cluster <- c(cluster, best); total <- new_tot; avail <- setdiff(avail, best)
            } else break
        }
        clusters[cluster] <- paste0("C", cid); to_assign <- setdiff(to_assign, cluster); cid <- cid + 1L
    }
    clusters
}
merge_small_clusters <- function(clu, death_tbl, nbrs_list, min_deaths, max_iter=20L){
    deaths_vec <- setNames(death_tbl$deaths, death_tbl$fips)
    unlabeled <- names(deaths_vec)[is.na(clu[names(deaths_vec)])]
    if (length(unlabeled)){
        next_id <- suppressWarnings(max(as.integer(sub("^C","", na.omit(unique(clu)))), na.rm = TRUE))
        if (!is.finite(next_id)) next_id <- 0L
        clu[unlabeled] <- paste0("C", seq(next_id + 1L, next_id + length(unlabeled)))
    }
    iter <- 0L
    repeat{
        iter <- iter + 1L
        sums <- tapply(deaths_vec, clu[names(deaths_vec)], sum, na.rm = TRUE); sums <- sums[is.finite(sums)]
        small <- names(sums)[sums < min_deaths]
        if (!length(small) || iter > max_iter) break
        for (sc in small){
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
fix_na_by_nearest <- function(clu, county_centroids){
    cty <- county_centroids %>% dplyr::select(GEOID, geometry); cty$cluster <- unname(clu[cty$GEOID])
    na_idx <- which(is.na(cty$cluster)); ok_idx <- which(!is.na(cty$cluster)); if (!length(na_idx)) return(clu)
    nn <- st_nearest_feature(cty[na_idx, ], cty[ok_idx, ])
    clu[cty$GEOID[na_idx]] <- cty$cluster[ok_idx[nn]]; clu
}
merge_small_clusters_by_distance <- function(clu, death_tbl, centroid_mat, min_deaths){
    deaths_vec <- setNames(death_tbl$deaths, death_tbl$fips); keep <- intersect(names(clu), rownames(centroid_mat))
    sums <- tapply(deaths_vec[keep], clu[keep], sum, na.rm = TRUE); sums <- sums[is.finite(sums)]
    small <- names(sums)[sums < min_deaths]; big <- names(sums)[sums >= min_deaths]
    if (!length(small) || !length(big)) return(clu)
    get_cluster_xy <- function(cid){ memb <- names(clu)[clu == cid]; pts <- centroid_mat[memb,,drop=FALSE]; colMeans(pts) }
    small_xy <- t(vapply(small, get_cluster_xy, numeric(2L))); big_xy <- t(vapply(big, get_cluster_xy, numeric(2L)))
    nearest_idx <- vapply(seq_len(nrow(small_xy)), function(i){ dif <- t(big_xy) - small_xy[i,]; which.min(colSums(dif*dif)) }, integer(1L))
    targets <- big[nearest_idx]
    for (i in seq_along(small)){ sc <- small[i]; tg <- targets[i]; members <- names(clu)[clu == sc]; if (!length(members)) next; clu[members] <- tg }
    clu
}

# -------------------- global K for normalization --------------------
log_msg("Computing global K for root3 & icd4 across ", years_all[1], "–", years_all[2])
con0 <- dbConnect(duckdb::duckdb(), dbdir=":memory:"); on.exit(try(dbDisconnect(con0, shutdown=TRUE), silent=TRUE), add=TRUE)
dbExecute(con0, paste0("PRAGMA threads = ", threads))
duckdb::duckdb_register(con0, "garbage", tibble(root3 = garbage_root_set))

dbExecute(con0, glue("
  CREATE TEMP TABLE cert_all AS
  SELECT
    row_number() OVER () AS id,
    LPAD(CAST({`county_var`} AS VARCHAR), 5, '0') AS fips,
    ucod,
    {`UCR39_COL`} AS ucr39
  FROM parquet_scan('{normalizePath(parquet_dir, winslash = '/') }/*.parquet')
  WHERE year BETWEEN {years_all[1]} AND {years_all[2]}
    AND ucod IS NOT NULL
    AND {`UCR39_COL`} IS NOT NULL
"))

dbExecute(con0, "
  CREATE TEMP TABLE ucod_valid_all AS
  SELECT id, fips,
         SUBSTR(UPPER(regexp_replace(ucod,'[^A-Za-z0-9]','')),1,3) AS root3,
         SUBSTR(UPPER(regexp_replace(ucod,'[^A-Za-z0-9]','')),1,4) AS icd4,
         ucr39
  FROM cert_all
  WHERE SUBSTR(UPPER(regexp_replace(ucod,'[^A-Za-z0-9]','')),1,3) NOT IN (SELECT root3 FROM garbage)
")

K_root3_global <- dbGetQuery(con0, "SELECT COUNT(DISTINCT root3) AS K FROM ucod_valid_all")$K
K_icd4_global  <- dbGetQuery(con0, "SELECT COUNT(DISTINCT icd4)  AS K FROM ucod_valid_all")$K
maxH_root3_global <- if (K_root3_global > 1) log(K_root3_global) else NA_real_
maxH_icd4_global  <- if (K_icd4_global  > 1) log(K_icd4_global)  else NA_real_
log_msg("K_global: root3=", K_root3_global, " | icd4=", K_icd4_global)
dbDisconnect(con0, shutdown=TRUE)


# ----- MCOD global K (icd4, non-underlying) for normalization -----
conK4 <- DBI::dbConnect(duckdb::duckdb(), dbdir=":memory:")
on.exit(try(DBI::dbDisconnect(conK4, shutdown=TRUE), silent=TRUE), add=TRUE)
DBI::dbExecute(conK4, paste0("PRAGMA threads = ", threads))
duckdb::duckdb_register(conK4, "garbage", tibble(root3 = garbage_root_set))

DBI::dbExecute(conK4, glue("
  CREATE TEMP TABLE cert_all AS
  SELECT row_number() OVER () AS id, ucod, {`UCR39_COL`} AS ucr39,
         {`county_var`} AS fips
  FROM parquet_scan('{normalizePath(parquet_dir, winslash = '/') }/*.parquet')
  WHERE year BETWEEN {years_all[1]} AND {years_all[2]}
    AND ucod IS NOT NULL
    AND {`UCR39_COL`} IS NOT NULL
"))

DBI::dbExecute(conK4, "
  CREATE TEMP TABLE ucod_valid_all AS
  SELECT id,
         SUBSTR(UPPER(regexp_replace(ucod,'[^A-Za-z0-9]','')),1,3) AS ucod_root3,
         ucr39
  FROM cert_all
  WHERE SUBSTR(UPPER(regexp_replace(ucod,'[^A-Za-z0-9]','')),1,3) NOT IN (SELECT root3 FROM garbage)
")

DBI::dbExecute(conK4, "CREATE TEMP TABLE mcod_codes_all (id BIGINT, icd4 VARCHAR, root3 VARCHAR)")
DBI::dbExecute(conK4, "
  INSERT INTO mcod_codes_all
  SELECT id,
         SUBSTR(UPPER(regexp_replace(code,'[^A-Za-z0-9]','')),1,4) AS icd4,
         SUBSTR(UPPER(regexp_replace(code,'[^A-Za-z0-9]','')),1,3) AS root3
  FROM (
    SELECT id, UNNEST(LIST_VALUE(record_1,record_2,record_3,record_4,record_5,
                                 record_6,record_7,record_8,record_9,record_10,
                                 record_11,record_12,record_13,record_14,record_15,
                                 record_16,record_17,record_18,record_19,record_20)) AS code
    FROM (
      SELECT row_number() OVER () AS id,
             record_1,record_2,record_3,record_4,record_5,
             record_6,record_7,record_8,record_9,record_10,
             record_11,record_12,record_13,record_14,record_15,
             record_16,record_17,record_18,record_19,record_20
      FROM parquet_scan('{normalizePath(parquet_dir, winslash = '/') }/*.parquet')
      WHERE year BETWEEN {years_all[1]} AND {years_all[2]}
    )
  )
  WHERE code IS NOT NULL AND code <> '' AND LENGTH(code) >= 3
")
DBI::dbExecute(conK4, "
  CREATE TEMP TABLE mcod_clean_all AS
  SELECT DISTINCT id, icd4, root3
  FROM mcod_codes_all
  WHERE root3 BETWEEN 'A00' AND 'Z99'
    AND SUBSTR(root3,2,1) BETWEEN '0' AND '9'
    AND root3 NOT IN (SELECT root3 FROM garbage)
")
DBI::dbExecute(conK4, "
  CREATE TEMP TABLE mcod_nonunderlying_all AS
  SELECT m.id, m.icd4
  FROM mcod_clean_all m
  LEFT JOIN ucod_valid_all u ON m.id = u.id AND m.root3 = u.ucod_root3
  WHERE u.ucod_root3 IS NULL
")
K_mcod_icd4_global <- DBI::dbGetQuery(conK4, "SELECT COUNT(DISTINCT icd4) AS K FROM mcod_nonunderlying_all")$K
maxH_mcod_icd4_global <- if (K_mcod_icd4_global > 1) log(K_mcod_icd4_global) else NA_real_
message("K_global (MCOD icd4, non-underlying) = ", K_mcod_icd4_global)
DBI::dbDisconnect(conK4, shutdown=TRUE)


# ----- MCOD global K (non-underlying) for normalization -----
conK <- DBI::dbConnect(duckdb::duckdb(), dbdir=":memory:")
on.exit(try(DBI::dbDisconnect(conK, shutdown=TRUE), silent=TRUE), add=TRUE)
DBI::dbExecute(conK, paste0("PRAGMA threads = ", threads))
duckdb::duckdb_register(conK, "garbage", tibble(root3 = garbage_root_set))

# Minimal cert with ucr39
DBI::dbExecute(conK, glue("
  CREATE TEMP TABLE cert_all AS
  SELECT row_number() OVER () AS id,
         ucod,
         {`UCR39_COL`} AS ucr39
  FROM parquet_scan('{normalizePath(parquet_dir, winslash = '/') }/*.parquet')
  WHERE year BETWEEN {years_all[1]} AND {years_all[2]}
    AND ucod IS NOT NULL
    AND {`UCR39_COL`} IS NOT NULL
"))

DBI::dbExecute(conK, "
  CREATE TEMP TABLE ucod_valid_all AS
  SELECT id,
         SUBSTR(UPPER(regexp_replace(ucod,'[^A-Za-z0-9]','')),1,3) AS ucod_root3,
         ucr39
  FROM cert_all
  WHERE SUBSTR(UPPER(regexp_replace(ucod,'[^A-Za-z0-9]','')),1,3) NOT IN (SELECT root3 FROM garbage)
")

# Extract MCOD root-3 per id, drop garbage, drop UCOD root if present
DBI::dbExecute(conK, "CREATE TEMP TABLE mcod_codes_all (id BIGINT, root3 VARCHAR)")
DBI::dbExecute(conK, "
  INSERT INTO mcod_codes_all
  SELECT id,
         SUBSTR(UPPER(regexp_replace(code,'[^A-Za-z0-9]','')),1,3) AS root3
  FROM (
    SELECT id, UNNEST(LIST_VALUE(record_1,record_2,record_3,record_4,record_5,
                                 record_6,record_7,record_8,record_9,record_10,
                                 record_11,record_12,record_13,record_14,record_15,
                                 record_16,record_17,record_18,record_19,record_20)) AS code
    FROM (
      SELECT row_number() OVER () AS id,
             record_1,record_2,record_3,record_4,record_5,
             record_6,record_7,record_8,record_9,record_10,
             record_11,record_12,record_13,record_14,record_15,
             record_16,record_17,record_18,record_19,record_20
      FROM parquet_scan('{normalizePath(parquet_dir, winslash = '/') }/*.parquet')
      WHERE year BETWEEN {years_all[1]} AND {years_all[2]}
    )
  )
  WHERE code IS NOT NULL AND code <> '' AND LENGTH(code) >= 3
")
DBI::dbExecute(conK, "
  CREATE TEMP TABLE mcod_clean_all AS
  SELECT DISTINCT id, root3
  FROM mcod_codes_all
  WHERE root3 BETWEEN 'A00' AND 'Z99'
    AND SUBSTR(root3,2,1) BETWEEN '0' AND '9'
    AND root3 NOT IN (SELECT root3 FROM garbage)
")
DBI::dbExecute(conK, "
  CREATE TEMP TABLE mcod_nonunderlying_all AS
  SELECT m.id, m.root3
  FROM mcod_clean_all m
  LEFT JOIN ucod_valid_all u ON m.id = u.id AND m.root3 = u.ucod_root3
  WHERE u.ucod_root3 IS NULL
")
K_mcod_global <- DBI::dbGetQuery(conK, "SELECT COUNT(DISTINCT root3) AS K FROM mcod_nonunderlying_all")$K
maxH_mcod_global <- if (K_mcod_global > 1) log(K_mcod_global) else NA_real_
message("K_global (MCOD root3, non-underlying) = ", K_mcod_global)
DBI::dbDisconnect(conK, shutdown=TRUE)


# -------------------- containers --------------------
cluster_membership <- list()
metrics_list       <- list()
std_ucr_list       <- list()
us_within_root3    <- list()
us_within_icd4     <- list()

# -------------------- main loop --------------------
for (pname in names(periods)) {
    log_msg("Processing ", pname)
    yrs <- range(periods[[pname]])
    stopifnot(dir_exists(parquet_dir))
    if (!length(dir(parquet_dir, pattern="\\.parquet$", full.names=TRUE)))
        stop("No .parquet files found in: ", parquet_dir)
    
    con <- dbConnect(duckdb::duckdb(), dbdir=":memory:")
    on.exit(try(dbDisconnect(con, shutdown=TRUE), silent=TRUE), add=TRUE)
    dbExecute(con, paste0("PRAGMA threads = ", threads))
    duckdb::duckdb_register(con, "garbage", tibble(root3 = garbage_root_set))
    
    dbExecute(con, glue("
    CREATE TEMP TABLE cert AS
    SELECT
      row_number() OVER () AS id,
      LPAD(CAST({`county_var`} AS VARCHAR), 5, '0') AS fips,
      ucod,
      {`UCR39_COL`} AS ucr39
    FROM parquet_scan('{normalizePath(parquet_dir, winslash = '/') }/*.parquet')
    WHERE year BETWEEN {yrs[1]} AND {yrs[2]}
      AND ucod IS NOT NULL
      AND {`UCR39_COL`} IS NOT NULL
  "))
    n_cert <- dbGetQuery(con, "SELECT COUNT(*) n FROM cert")$n
    log_msg("[", pname, "] rows in cert: ", format(n_cert, big.mark=","))
    
    dbExecute(con, "
    CREATE TEMP TABLE ucod_valid AS
    SELECT id, fips,
           SUBSTR(UPPER(regexp_replace(ucod,'[^A-Za-z0-9]','')),1,3) AS root3,
           SUBSTR(UPPER(regexp_replace(ucod,'[^A-Za-z0-9]','')),1,4) AS icd4,
           ucr39
    FROM cert
    WHERE SUBSTR(UPPER(regexp_replace(ucod,'[^A-Za-z0-9]','')),1,3) NOT IN (SELECT root3 FROM garbage)
  ")
    
    # deaths per county for clustering
    death_tbl <- dbGetQuery(con, "SELECT fips, COUNT(*) AS deaths FROM ucod_valid GROUP BY fips") %>%
        as_tibble() %>% mutate(fips = stringr::str_pad(as.character(fips), 5, pad="0")) %>%
        filter(fips %in% counties_sf$GEOID)
    
    # clusters
    clu <- make_greedy_clusters(death_tbl, nbrs_list, min_deaths)
    clu <- merge_small_clusters(clu, death_tbl, nbrs_list, min_deaths)
    clu <- fix_na_by_nearest(clu, county_centroids)
    deaths_vec <- setNames(death_tbl$deaths, death_tbl$fips)
    sums_post <- tapply(deaths_vec, clu[names(deaths_vec)], sum, na.rm=TRUE)
    if (any(is.finite(sums_post) & (sums_post < min_deaths))) {
        clu <- merge_small_clusters_by_distance(clu, death_tbl, centroid_mat, min_deaths)
        clu <- fix_na_by_nearest(clu, county_centroids)
    }
    clu_df <- tibble(fips = names(clu), cluster = unname(clu)) %>% left_join(death_tbl, by="fips")
    cluster_membership[[pname]] <- clu_df %>% mutate(period = pname)
    
    # map for joins
    dbExecute(con, "DROP TABLE IF EXISTS map_fips_cluster")
    dbExecute(con, "CREATE TEMP TABLE map_fips_cluster (fips VARCHAR, cluster VARCHAR)")
    DBI::dbAppendTable(con, "map_fips_cluster", clu_df %>% select(fips, cluster) %>% as.data.frame())
    
    # ----- US standard shares over UCR-39 (per period) -----
    std_ucr <- dbGetQuery(con, "
    SELECT ucr39, COUNT(*) AS n
    FROM ucod_valid
    GROUP BY 1
  ") %>%
        as_tibble() %>%
        mutate(period = pname, s = n / sum(n))
    std_ucr_list[[pname]] <- std_ucr %>% select(period, ucr39, s)
    
    # ----- National within-cause (fallback) shares for root3 & icd4 -----
    us_r3 <- dbGetQuery(con, "
    SELECT ucr39, root3 AS domain, COUNT(*) AS n
    FROM ucod_valid
    GROUP BY 1,2
  ") %>% as_tibble() %>%
        group_by(ucr39) %>% mutate(w_us = n / sum(n)) %>% ungroup() %>%
        mutate(period = pname) %>% select(period, ucr39, domain, w_us)
    
    us_r4 <- dbGetQuery(con, "
    SELECT ucr39, icd4 AS domain, COUNT(*) AS n
    FROM ucod_valid
    GROUP BY 1,2
  ") %>% as_tibble() %>%
        group_by(ucr39) %>% mutate(w_us = n / sum(n)) %>% ungroup() %>%
        mutate(period = pname) %>% select(period, ucr39, domain, w_us)
    # ---- MCOD (non-underlying) by (cluster, ucr39, root3) ----
    # 1) Build MCOD tables in this period
    DBI::dbExecute(con, "DROP TABLE IF EXISTS mcod_codes_raw")
    DBI::dbExecute(con, "CREATE TEMP TABLE mcod_codes_raw (id BIGINT, fips VARCHAR, root3 VARCHAR)")
    
    # Pull all 20 record fields; drop garbage early
    DBI::dbExecute(con, "
  INSERT INTO mcod_codes_raw
  SELECT id, fips,
         SUBSTR(UPPER(regexp_replace(code,'[^A-Za-z0-9]','')),1,3) AS root3
  FROM (
    SELECT id, fips, UNNEST(LIST_VALUE(record_1,record_2,record_3,record_4,record_5,
                                       record_6,record_7,record_8,record_9,record_10,
                                       record_11,record_12,record_13,record_14,record_15,
                                       record_16,record_17,record_18,record_19,record_20)) AS code
    FROM (
      SELECT row_number() OVER () AS id,
             LPAD(CAST({`county_var`} AS VARCHAR), 5, '0') AS fips,
             record_1,record_2,record_3,record_4,record_5,
             record_6,record_7,record_8,record_9,record_10,
             record_11,record_12,record_13,record_14,record_15,
             record_16,record_17,record_18,record_19,record_20
      FROM parquet_scan('{normalizePath(parquet_dir, winslash = '/') }/*.parquet')
      WHERE year BETWEEN {yrs[1]} AND {yrs[2]}
    )
  )
  WHERE code IS NOT NULL AND code <> '' AND LENGTH(code) >= 3
")
    DBI::dbExecute(con, "
  CREATE TEMP TABLE mcod_clean AS
  SELECT DISTINCT id, fips, root3
  FROM mcod_codes_raw
  WHERE root3 BETWEEN 'A00' AND 'Z99'
    AND SUBSTR(root3,2,1) BETWEEN '0' AND '9'
    AND root3 NOT IN (SELECT root3 FROM garbage)
")
    
    # Drop underlying root if it appears among MCOD for that id
    DBI::dbExecute(con, "
  CREATE TEMP TABLE mcod_nonunderlying AS
  SELECT m.id, m.fips, m.root3
  FROM mcod_clean m
  LEFT JOIN ucod_valid u
    ON m.id = u.id AND m.root3 = u.root3
  WHERE u.root3 IS NULL
")
    
    # 2) Join MCOD to the death's ucr39 via id, then roll up by cluster
    DBI::dbExecute(con, "
  CREATE TEMP TABLE mcod_with_cause AS
  SELECT n.id, n.fips, n.root3, u.ucr39
  FROM mcod_nonunderlying n
  JOIN ucod_valid u ON n.id = u.id
")
    
    # Per (cluster, ucr39, root3): deaths with that MCOD present
    mcod_counts <- DBI::dbGetQuery(con, "
  SELECT m.cluster, w.ucr39, w.root3 AS domain, COUNT(DISTINCT w.id) AS n
  FROM mcod_with_cause w
  JOIN map_fips_cluster m ON w.fips = m.fips
  WHERE m.cluster IS NOT NULL
  GROUP BY 1,2,3
") |> tibble::as_tibble() |> dplyr::mutate(period = pname)
    
    # 3) National within-cause MCOD shares (fallback) at root3
    us_mcod_r3 <- DBI::dbGetQuery(con, "
  SELECT ucr39, root3 AS domain, COUNT(DISTINCT id) AS n
  FROM mcod_with_cause
  GROUP BY 1,2
") |> tibble::as_tibble() |>
        dplyr::group_by(ucr39) |>
        dplyr::mutate(w_us = n / sum(n)) |>
        dplyr::ungroup() |>
        dplyr::mutate(period = pname) |>
        dplyr::select(period, ucr39, domain, w_us)
    
    # 4) Cause-standardized MCOD detail (reuse helper)
    det_mcod_r3 <- compute_cstd_detail(
        counts_df    = mcod_counts |> dplyr::select(period, cluster, ucr39, domain, n),
        std_ucr_df   = std_ucr     |> dplyr::select(period, ucr39, s),
        us_within_df = us_mcod_r3,
        maxH_global  = maxH_mcod_global
    ) |> dplyr::mutate(metric = "detail_mcod_root3_cstd")
    
    # 5) Add to metrics
    metrics_list[[pname]] <- dplyr::bind_rows(metrics_list[[pname]], det_mcod_r3)
    
    # ---- MCOD icd4 (non-underlying) by (cluster, ucr39, icd4) ----
    # Reuse the mcod_nonunderlying already built above, but derive icd4 alongside root3
    DBI::dbExecute(con, "DROP TABLE IF EXISTS mcod_nonunderlying_icd4")
    DBI::dbExecute(con, "
  CREATE TEMP TABLE mcod_nonunderlying_icd4 AS
  SELECT DISTINCT n.id, n.fips,
         SUBSTR(UPPER(regexp_replace(code,'[^A-Za-z0-9]','')),1,4) AS icd4,
         SUBSTR(UPPER(regexp_replace(code,'[^A-Za-z0-9]','')),1,3) AS root3
  FROM (
    SELECT id, fips, UNNEST(LIST_VALUE(record_1,record_2,record_3,record_4,record_5,
                                       record_6,record_7,record_8,record_9,record_10,
                                       record_11,record_12,record_13,record_14,record_15,
                                       record_16,record_17,record_18,record_19,record_20)) AS code
    FROM (
      SELECT row_number() OVER () AS id,
             LPAD(CAST({`county_var`} AS VARCHAR), 5, '0') AS fips,
             record_1,record_2,record_3,record_4,record_5,
             record_6,record_7,record_8,record_9,record_10,
             record_11,record_12,record_13,record_14,record_15,
             record_16,record_17,record_18,record_19,record_20
      FROM parquet_scan('{normalizePath(parquet_dir, winslash = '/') }/*.parquet')
      WHERE year BETWEEN {yrs[1]} AND {yrs[2]}
    )
  ) z
  JOIN mcod_nonunderlying nu ON z.id = nu.id AND SUBSTR(UPPER(regexp_replace(z.code,'[^A-Za-z0-9]','')),1,3) = nu.root3
  WHERE z.icd4 IS NOT NULL AND LENGTH(z.icd4) = 4
    AND z.root3 BETWEEN 'A00' AND 'Z99'
    AND SUBSTR(z.root3,2,1) BETWEEN '0' AND '9'
    AND z.root3 NOT IN (SELECT root3 FROM garbage)
")
    
    DBI::dbExecute(con, "DROP TABLE IF EXISTS mcod_with_cause_icd4")
    DBI::dbExecute(con, "
  CREATE TEMP TABLE mcod_with_cause_icd4 AS
  SELECT n.id, n.fips, n.icd4, u.ucr39
  FROM mcod_nonunderlying_icd4 n
  JOIN ucod_valid u ON n.id = u.id
")
    
    mcod4_counts <- DBI::dbGetQuery(con, "
  SELECT m.cluster, w.ucr39, w.icd4 AS domain, COUNT(DISTINCT w.id) AS n
  FROM mcod_with_cause_icd4 w
  JOIN map_fips_cluster m ON w.fips = m.fips
  WHERE m.cluster IS NOT NULL
  GROUP BY 1,2,3
") |> tibble::as_tibble() |> dplyr::mutate(period = pname)
    
    us_mcod_icd4 <- DBI::dbGetQuery(con, "
  SELECT ucr39, icd4 AS domain, COUNT(DISTINCT id) AS n
  FROM mcod_with_cause_icd4
  GROUP BY 1,2
") |> tibble::as_tibble() |>
        dplyr::group_by(ucr39) |>
        dplyr::mutate(w_us = n / sum(n)) |>
        dplyr::ungroup() |>
        dplyr::mutate(period = pname) |>
        dplyr::select(period, ucr39, domain, w_us)
    
    det_mcod_icd4 <- compute_cstd_detail(
        counts_df    = mcod4_counts |> dplyr::select(period, cluster, ucr39, domain, n),
        std_ucr_df   = std_ucr      |> dplyr::select(period, ucr39, s),
        us_within_df = us_mcod_icd4,
        maxH_global  = maxH_mcod_icd4_global
    ) |> dplyr::mutate(metric = "detail_mcod_icd4_cstd")
    
    metrics_list[[pname]] <- dplyr::bind_rows(metrics_list[[pname]], det_mcod_icd4)
    
    
    us_within_root3[[pname]] <- us_r3
    us_within_icd4[[pname]]  <- us_r4
    
    # ----- Cluster-level counts by (ucr39, domain) -----
    counts_r3 <- dbGetQuery(con, "
    SELECT m.cluster, v.ucr39, v.root3 AS domain, COUNT(*) AS n
    FROM ucod_valid v
    JOIN map_fips_cluster m ON v.fips = m.fips
    WHERE m.cluster IS NOT NULL
    GROUP BY 1,2,3
  ") %>% as_tibble() %>% mutate(period = pname, domain_type = "root3")
    
    counts_r4 <- dbGetQuery(con, "
    SELECT m.cluster, v.ucr39, v.icd4 AS domain, COUNT(*) AS n
    FROM ucod_valid v
    JOIN map_fips_cluster m ON v.fips = m.fips
    WHERE m.cluster IS NOT NULL
    GROUP BY 1,2,3
  ") %>% as_tibble() %>% mutate(period = pname, domain_type = "icd4")
    
    # ----- Cause-standardized detail helper -----
    compute_cstd_detail <- function(counts_df, std_ucr_df, us_within_df, maxH_global){
        if (!nrow(counts_df)) return(tibble(period=character(), cluster=character(), detail=numeric()))
        # cluster within-cause shares
        w_k <- counts_df %>%
            group_by(period, cluster, ucr39) %>%
            mutate(w = n / sum(n)) %>% ungroup() %>%
            select(period, cluster, ucr39, domain, w)
        
        # base: all (ucr39, domain) combos from national within-cause + attach s
        base <- us_within_df %>% left_join(std_ucr_df, by=c("period","ucr39"))
        
        # expand to ensure every cluster has every (ucr39, domain); fallback to national w_us
        clusters_all <- distinct(counts_df, cluster) %>% pull(cluster)
        expanded <- base %>%
            left_join(w_k, by=c("period","ucr39","domain")) %>%
            group_by(period, ucr39, domain) %>%
            tidyr::complete(cluster = clusters_all, fill = list(w = NA_real_)) %>%
            ungroup() %>%
            mutate(w_eff = dplyr::coalesce(w, w_us))
        
        # p*(domain) = sum_c s(c) * w_eff(domain|c)
        p_star <- expanded %>%
            group_by(period, cluster, domain) %>%
            summarise(p = sum(s * w_eff, na.rm=TRUE), .groups = "drop")
        
        # Entropy
        detail <- p_star %>%
            group_by(period, cluster) %>%
            summarise(H = {
                pp <- p / sum(p, na.rm=TRUE)
                pp <- pp[is.finite(pp) & pp > 0]
                if (!length(pp)) NA_real_ else -sum(pp * log(pp))
            }, .groups = "drop") %>%
            mutate(detail = .clamp100(100 * H / maxH_global)) %>%
            select(period, cluster, detail)
        detail
    }
    
    det_r3 <- compute_cstd_detail(
        counts_df    = counts_r3 %>% select(period, cluster, ucr39, domain, n),
        std_ucr_df   = std_ucr %>% select(period, ucr39, s),
        us_within_df = us_r3,
        maxH_global  = maxH_root3_global
    ) %>% mutate(metric = "detail_ucod_root3_cstd")
    
    det_r4 <- compute_cstd_detail(
        counts_df    = counts_r4 %>% select(period, cluster, ucr39, domain, n),
        std_ucr_df   = std_ucr %>% select(period, ucr39, s),
        us_within_df = us_r4,
        maxH_global  = maxH_icd4_global
    ) %>% mutate(metric = "detail_ucod_icd4_cstd")
    
    metrics_list[[pname]] <- bind_rows(det_r3, det_r4)
    
    dbDisconnect(con, shutdown=TRUE); gc()
}

# -------------------- outputs --------------------
cluster_members <- bind_rows(cluster_membership) %>%
    transmute(period, cluster, fips, cluster_deaths = deaths)

metrics_out <- bind_rows(metrics_list) %>%
    tidyr::pivot_wider(names_from = metric, values_from = detail) %>%
    arrange(period, cluster)

std_ucr_all      <- bind_rows(std_ucr_list) %>% arrange(period, ucr39)
us_within_r3_all <- bind_rows(us_within_root3) %>% arrange(period, ucr39, domain)
us_within_r4_all <- bind_rows(us_within_icd4) %>% arrange(period, ucr39, domain)

write_csv(metrics_out,      file.path(out_dir, "cluster_metrics_ucr39_cstd.csv.gz"))
write_csv(cluster_members,  file.path(out_dir, "county_cluster_membership.csv.gz"))
write_csv(std_ucr_all,      file.path(out_dir, "ucr39_standard_by_period.csv.gz"))
write_csv(us_within_r3_all, file.path(out_dir, "national_within_cause_shares_root3.csv.gz"))
write_csv(us_within_r4_all, file.path(out_dir, "national_within_cause_shares_icd4.csv.gz"))

if (nrow(metrics_out)) {
    cat("\n— Summary counts by period —\n")
    metrics_out %>%
        summarise(.by = period,
                  n_clusters = n_distinct(cluster),
                  nonNA_cstd_root3 = sum(!is.na(detail_ucod_root3_cstd)),
                  nonNA_cstd_icd4  = sum(!is.na(detail_ucod_icd4_cstd))
        ) %>% print(n=Inf)
    cat("\nK (global) used: root3 =", K_root3_global, " | icd4 =", K_icd4_global, "\n")
}
cat("\nDone.\n")

