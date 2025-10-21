# --------------------------------------------------------------
# ANACONDA detail metric
# Date: August 18, 2025
#
# This script is an implementation of "detail" metric from the Phillips ANACONDA paper.
# Because most counties have far too few deaths for this metric, they were put into clusters 
# using the US cluster groups. Then for clusters with less than 2000 deaths adjacent ones were merged until reaching 2000 
# across all time periods
# The code computes Shannon entropy over GBD cause categories and converts them to the 
# effective number of causes of death. It then standardizes to 2000 deaths by taking a sample of 2000
# when there is more than 2000 deaths. 
#
# The only difference between this metric and Phillips' is that they use
# reference size as the sample size at which all populations reach 95% sample completeness. 
# Unfortunately, our clusters are too small for that to be feasible, but plots (see slides) 
# show that large cluster size effects seem to drop-off around 2000 deaths making it a justifiable reference size. 
# Phillips et al note in their paper that size based rarefaction is appropriate for small places. 
#
# It is standardized by ucr39 cause of death categories.
#
# --------------------------------------------------------------

suppressPackageStartupMessages({
    library(fs); library(readr); library(dplyr); library(stringr)
    library(purrr); library(tidyr); library(tibble)
    library(sf); library(tigris); library(igraph)
    library(here); library(DBI); library(duckdb); library(glue)
})

min_deaths <- 2000L
periods <- list(
    "1999_2005" = 1999:2005,
    "2006_2012" = 2006:2012,
    "2013_2019" = 2013:2019,
    "2020_2022" = 2020:2022
)

out_dir   <- here("output"); dir_create(out_dir)
diag_out_dir <- here("output","debug_missing_fips_diagnostics"); dir_create(diag_out_dir, recurse = TRUE)
parquet_dir <- if (dir_exists(here("data_private","mcod"))) here("data_private","mcod") else here("data_private","mcod_sample")

garbage_csv <- here("data_raw","cause-codes","gbd_garbage_codes_without_overdose.csv")
stopifnot(file.exists(garbage_csv))

log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "│", paste0(..., collapse=""), "\n")
CLAMP_0_100 <- TRUE
.clamp100 <- function(x){ if(!CLAMP_0_100) return(x); y <- x; y[x<0]<-0; y[x>100]<-100; y }

county_var <- "county_ihme"
UCR39_COL  <- "ucr39"

crs_proj  <- 5070
threads   <- max(1L, parallel::detectCores() - 1L)

# ---------------- FIPS patch ----------------
fips_patch <- c(
    "02201" = "02158", "02232" = "02275", "02261" = "02164",
    "02270" = "02158", "02280" = "02105", "02282" = "02195", "02290" = "02122",
    "46113" = "46102",
    "51515" = "51019", "51560" = "51163"
)
apply_fips_patch <- function(f){
    fch <- as.character(f)
    mapped <- ifelse(!is.na(fch) & fch %in% names(fips_patch), unname(fips_patch[fch]), fch)
    mapped <- stringr::str_trim(mapped)
    mapped <- gsub("[^0-9]", "", mapped)
    mapped[mapped == ""] <- NA_character_
    mapped <- ifelse(is.na(mapped), NA_character_, stringr::str_pad(mapped, 5, pad = "0"))
    mapped
}

normalize_fips <- function(x) {
    x <- as.character(x)
    x <- stringr::str_trim(x)
    x <- gsub("[^0-9]", "", x)
    x[x == ""] <- NA_character_
    x <- ifelse(is.na(x), NA_character_, stringr::str_pad(x, width = 5, side = "left", pad = "0"))
    x
}

local_shapefile_dir <- here::here("data_raw", "tigris_counties_local")
if (!dir.exists(local_shapefile_dir)) local_shapefile_dir <- NULL

try_read_local_shapefile <- function(year) {
    if (is.null(local_shapefile_dir)) return(NULL)
    pattern <- paste0("cb_", year, "_us_county_")
    candidates <- list.files(local_shapefile_dir, pattern = pattern, full.names = TRUE)
    shp <- candidates[grepl("\\.shp$", candidates)]
    if (length(shp) == 0) return(NULL)
    tryCatch(sf::st_read(shp[1], quiet = TRUE), error = function(e) NULL)
}

# ---------------- GEOID helper ----------------
ensure_geoid <- function(sf_obj) {
    nm <- names(sf_obj)
    if ("GEOID" %in% nm) return(sf_obj)
    if ("GEOID10" %in% nm) { sf_obj <- dplyr::rename(sf_obj, GEOID = GEOID10); return(sf_obj) }
    if ("GEOID20" %in% nm) { sf_obj <- dplyr::rename(sf_obj, GEOID = GEOID20); return(sf_obj) }
    if (all(c("STATEFP","COUNTYFP") %in% nm)) {
        sf_obj <- sf_obj %>% mutate(GEOID = paste0(STATEFP, COUNTYFP))
        return(sf_obj)
    }
    if (all(c("STATE","COUNTY") %in% nm)) {
        st <- as.character(sf_obj$STATE); ct <- as.character(sf_obj$COUNTY)
        st <- stringr::str_pad(gsub("\\D","",st), 2, pad="0"); ct <- stringr::str_pad(gsub("\\D","",ct), 3, pad="0")
        sf_obj <- sf_obj %>% mutate(GEOID = paste0(st, ct))
        return(sf_obj)
    }
    stop("Shapefile lacks GEOID / STATEFP+COUNTYFP columns and cannot auto-construct GEOID.")
}

# ---------------- Robust tigris cache & optional downloader ----------------
options(tigris_use_cache = TRUE)
tigris_cache_dir <- here::here("data_raw","tigris_cache")
if (!dir.exists(tigris_cache_dir)) dir.create(tigris_cache_dir, recursive = TRUE)
Sys.setenv(TIGRIS_CACHE_DIR = tigris_cache_dir)

download_census_cb <- function(year, dest_dir = here::here("data_raw","tigris_counties_local")) {
    if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
    destfile <- file.path(dest_dir, paste0("cb_", year, "_us_county_500k.zip"))
    urls <- c(
        sprintf("https://www2.census.gov/geo/tiger/GENZ%d/cb_%d_us_county_500k.zip", year, year),
        sprintf("https://www2.census.gov/geo/tiger/GENZ%d/cb_%d_us_county_20m.zip", year, year),
        sprintf("https://www2.census.gov/geo/tiger/GENZ%d/shp/cb_%d_us_county_500k.zip", year, year)
    )
    for (u in urls) {
        try({
            message("Trying download: ", u)
            utils::download.file(u, destfile, mode = "wb", quiet = TRUE)
            if (file.exists(destfile) && file.info(destfile)$size > 2048) {
                unzip(destfile, exdir = dest_dir)
                message("Downloaded and unzipped to ", dest_dir)
                return(TRUE)
            }
        }, silent = TRUE)
    }
    message("All download attempts failed for year ", year, ". Please place shapefile into ", dest_dir, " manually.")
    FALSE
}

# ---------------- Robust spatial builder (bidirectional fallback + cache + local) ----------------
build_year_spatial_pref_earlier_geoid_fallback <- function(yr, crs_proj = 5070, fallback_window = 10, min_year = 1990) {
    sf::sf_use_s2(FALSE)
    
    safe_normalize_counties <- function(counties_sf) {
        if (!inherits(counties_sf, "sf")) stop("safe_normalize_counties needs an sf object")
        if (nrow(counties_sf) == 0) stop("empty sf object provided")
        if (is.null(sf::st_geometry(counties_sf))) stop("no geometry column")
        if ("st_zm" %in% ls(getNamespace("sf"), all.names = TRUE)) {
            try({
                geom <- sf::st_geometry(counties_sf)
                sf::st_geometry(counties_sf) <- sf::st_zm(geom, drop = TRUE, what = "ZM")
            }, silent = TRUE)
        }
        try({
            is_empty <- sf::st_is_empty(counties_sf)
            if (any(is_empty, na.rm = TRUE)) counties_sf <- counties_sf[!is_empty, , drop = FALSE]
        }, silent = TRUE)
        try({
            counties_sf <- sf::st_cast(counties_sf, "MULTIPOLYGON", warn = FALSE)
        }, silent = TRUE)
        try({
            geoms_valid <- sf::st_make_valid(sf::st_geometry(counties_sf))
            sf::st_geometry(counties_sf) <- geoms_valid
        }, silent = TRUE)
        try({
            keep_idx <- which(!is.na(sf::st_geometry(counties_sf)) & !sf::st_is_empty(sf::st_geometry(counties_sf)))
            if (length(keep_idx) == 0) stop("No valid geometries after normalization")
            if (length(keep_idx) < nrow(counties_sf)) counties_sf <- counties_sf[keep_idx, , drop = FALSE]
        }, silent = TRUE)
        counties_sf
    }
    
    safe_st_transform <- function(x, crs_proj) {
        out <- tryCatch(sf::st_transform(x, crs_proj), error = function(e) e)
        if (!inherits(out, "error")) return(out)
        message("    safe_st_transform: initial st_transform failed: ", out$message)
        try({
            x2 <- sf::st_make_valid(x)
            if ("st_zm" %in% ls(getNamespace("sf"), all.names = TRUE)) {
                sf::st_geometry(x2) <- sf::st_zm(sf::st_geometry(x2), drop = TRUE, what = "ZM")
            }
            out2 <- tryCatch(sf::st_transform(x2, crs_proj), error = function(e) e)
            if (!inherits(out2, "error")) return(out2)
        }, silent = TRUE)
        try({
            x3 <- tryCatch(sf::st_cast(x, "MULTIPOLYGON", warn = FALSE), error = function(e) x)
            out3 <- tryCatch(sf::st_transform(x3, crs_proj), error = function(e) e)
            if (!inherits(out3, "error")) return(out3)
        }, silent = TRUE)
        stop("safe_st_transform: unable to transform geometry to target CRS (", crs_proj, ")")
    }
    
    try_download_tigris <- function(y) {
        tryCatch({
            message("    trying tigris::counties(cb=TRUE, year=", y, ") via HTTPS...")
            ct <- tigris::counties(cb = TRUE, year = y)
            if (is.null(ct) || nrow(ct) == 0) stop("empty result")
            ct
        }, error = function(e) {
            message("      tigris download failed for year ", y, ": ", e$message)
            NULL
        })
    }
    
    try_local <- function(year) {
        try_read_local_shapefile(year)
    }
    
    local_sf <- try_local(yr)
    year_used <- NA_integer_
    if (!is.null(local_sf)) {
        counties_sf <- local_sf
        year_used <- yr
    } else {
        counties_sf <- try_download_tigris(yr)
        if (!is.null(counties_sf)) year_used <- yr
        
        if (is.null(counties_sf)) {
            years_before <- seq(from = yr - 1, to = max(min_year, yr - fallback_window), by = -1)
            years_after  <- seq(from = yr + 1, to = yr + fallback_window, by = 1)
            nmax <- max(length(years_before), length(years_after))
            interleaved <- c()
            for (i in seq_len(nmax)) {
                if (i <= length(years_before)) interleaved <- c(interleaved, years_before[i])
                if (i <= length(years_after))  interleaved <- c(interleaved,  years_after[i])
            }
            years_try_candidates <- unique(interleaved)
            for (y2 in years_try_candidates) {
                counties_sf <- try_download_tigris(y2)
                if (!is.null(counties_sf)) { year_used <- y2; break }
                lf <- try_local(y2)
                if (!is.null(lf)) { counties_sf <- lf; year_used <- y2; break }
            }
        }
        
        if (is.null(counties_sf)) {
            try({
                message("    attempting tigris cache-only read for year ", yr, " ...")
                counties_sf <- tigris::counties(cb = TRUE, year = yr, cache = TRUE)
                if (!is.null(counties_sf) && nrow(counties_sf) > 0) year_used <- yr
            }, silent = TRUE)
            if (is.null(counties_sf)) {
                for (y2 in c(seq(yr-1, max(min_year, yr-fallback_window), -1), seq(yr+1, yr+fallback_window, 1))) {
                    try({
                        counties_sf <- tigris::counties(cb = TRUE, year = y2, cache = TRUE)
                        if (!is.null(counties_sf) && nrow(counties_sf) > 0) { year_used <- y2; break }
                    }, silent = TRUE)
                }
            }
        }
    }
    
    if (is.null(counties_sf) || nrow(counties_sf) == 0) {
        stop("Failed to obtain county shapefile for year ", yr, " or nearby fallback years. Provide a local shapefile.")
    }
    
    counties_sf <- tryCatch({ ensure_geoid(counties_sf) }, error = function(e) { stop(e$message) })
    counties_sf <- counties_sf %>% mutate(GEOID = as.character(GEOID), GEOID = stringr::str_pad(GEOID, 5, pad = "0")) %>% dplyr::select(GEOID, geometry)
    
    counties_sf <- tryCatch({
        safe_normalize_counties(counties_sf)
    }, error = function(e) {
        message("    Warning: geometry normalization failed: ", e$message, " — proceeding with raw shapefile (may fail later).")
        counties_sf
    })
    
    counties_proj <- tryCatch({
        safe_st_transform(counties_sf, crs_proj)
    }, error = function(e) {
        message("    safe_st_transform ultimately failed: ", e$message)
        stop("st_transform failed and recovery attempts exhausted for year ", yr)
    })
    
    adj <- tryCatch({
        sf::st_touches(counties_proj)
    }, error = function(e) {
        message("    st_touches failed: ", e$message, " — trying st_intersects() as fallback")
        tryCatch(sf::st_intersects(counties_proj), error = function(e2) stop("st_intersects also failed: ", e2$message))
    })
    
    edge_df <- tibble::tibble(from = rep(counties_sf$GEOID, lengths(adj)),
                              to   = counties_sf$GEOID[unlist(adj)]) %>% filter(from < to)
    g <- tryCatch(igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = counties_sf$GEOID),
                  error = function(e) stop("igraph construction failed: ", e$message))
    nbrs_list <- setNames(lapply(V(g)$name, function(v) names(neighbors(g, v, mode="all"))), V(g)$name)
    
    county_centroids <- tryCatch({
        sf::st_centroid(counties_proj) %>% mutate(GEOID = counties_sf$GEOID) %>% dplyr::select(GEOID, geometry)
    }, error = function(e) {
        stop("st_centroid failed: ", e$message)
    })
    
    coords <- tryCatch(sf::st_coordinates(county_centroids), error = function(e) stop("st_coordinates failed: ", e$message))
    centroid_mat <- coords[, c("X","Y"), drop = FALSE]; rownames(centroid_mat) <- county_centroids$GEOID
    
    tryCatch({
        readr::write_csv(tibble::tibble(requested_year = yr, year_used = year_used),
                         file.path(diag_out_dir, paste0("year_lookup_", yr, ".csv")))
    }, silent = TRUE)
    
    message("    spatial built for requested year ", yr, "; year used: ", year_used)
    list(requested_year = yr, year_used = year_used, counties_sf = counties_sf, centroid_mat = centroid_mat, nbrs_list = nbrs_list)
}

period_spatial_list <- list()
for (pname in names(periods)) {
    yrs <- periods[[pname]]
    message("Building spatials for period:", pname, "years:", min(yrs), "-", max(yrs))
    year_objs <- lapply(yrs, function(y) {
        if (y == 2022) {
            message("  skipping 2022 for Connecticut — using 2021 shapefile instead (planning region change)")
            y_use <- 2021
        } else {
            y_use <- y
        }
        
        message("  building spatial for year ", y, " (effective shapefile year: ", y_use, ")")
        build_year_spatial_pref_earlier_geoid_fallback(y_use, crs_proj = crs_proj, fallback_window = 10, min_year = 1990)
    })    
    names(year_objs) <- as.character(yrs)
    fips_by_year <- lapply(year_objs, function(x) rownames(x$centroid_mat))
    stable_fips <- Reduce(intersect, fips_by_year)
    readr::write_csv(tibble::tibble(year = as.character(yrs), n_counties = vapply(fips_by_year, length, integer(1))),
                     file.path(diag_out_dir, paste0("period_year_counts_", pname, ".csv")))
    readr::write_csv(tibble::tibble(stable_fips = stable_fips),
                     file.path(diag_out_dir, paste0("temporally_stable_fips_", pname, ".csv")))
    period_spatial_list[[pname]] <- list(years = yrs, year_objs = year_objs, stable_fips = stable_fips)
    message("  temporally-stable counties for ", pname, ": ", length(stable_fips))
}

# ---------------- Build ONE global cluster mapping (≥ min_deaths aggregated across all periods) ----------------
message("Building global cluster mapping with ≥", min_deaths, " deaths aggregated across all periods")

# canonical year for centroids/adjacency: earliest available year among retrieved year_objs
all_year_objs <- unlist(lapply(period_spatial_list, function(p) names(p$year_objs)))
all_years_sorted <- sort(unique(all_year_objs))
if (length(all_years_sorted) == 0) stop("No year objects available in period_spatial_list — spatial build failed.")
canonical_year_for_global <- all_years_sorted[1]
message("Using canonical year for centroids/adjacency: ", canonical_year_for_global)

# obtain year_obj that contains canonical_year_for_global
any_year_obj <- NULL
for (p in period_spatial_list) {
    if (canonical_year_for_global %in% names(p$year_objs)) { any_year_obj <- p$year_objs[[canonical_year_for_global]]; break }
}
if (is.null(any_year_obj)) stop("Could not find canonical year spatial for global mapping: ", canonical_year_for_global)
centroid_mat_global <- any_year_obj$centroid_mat
nbrs_list_global    <- any_year_obj$nbrs_list
valid_fips_by_centroid <- rownames(centroid_mat_global)
if (length(valid_fips_by_centroid) == 0) stop("No centroid rows available in canonical year's centroid matrix")

# temporally-stable intersection across all periods
stable_fips_each_period <- lapply(period_spatial_list, function(x) x$stable_fips)
valid_fips_global <- Reduce(intersect, stable_fips_each_period)
message("Number of temporally-stable FIPS in intersection (used for global mapping): ", length(valid_fips_global))
if (length(valid_fips_global) == 0) stop("No temporally-stable counties common to all periods - cannot build global mapping with intersection. Consider using union instead.")
valid_fips_global <- intersect(valid_fips_global, valid_fips_by_centroid)
message("After intersecting with centroid GEOIDs: ", length(valid_fips_global), " counties available for global mapping")
if (length(valid_fips_global) == 0) stop("After intersecting with centroid GEOIDs, no valid fips remain for global mapping — check shapefiles.")

# aggregated deaths across whole study (all periods) using DuckDB / parquet
con_tmp <- dbConnect(duckdb::duckdb(), dbdir=":memory:")
dbExecute(con_tmp, paste0("PRAGMA threads = ", threads))
yrs_all <- range(unlist(periods))
dbExecute(con_tmp, glue::glue("
  CREATE TEMP TABLE base_all AS
  SELECT LPAD(CAST({`county_var`} AS VARCHAR), 5, '0') AS fips
  FROM parquet_scan('{normalizePath(parquet_dir, winslash = '/') }/*.parquet')
  WHERE year BETWEEN {yrs_all[1]} AND {yrs_all[2]}
    AND {`UCR39_COL`} IS NOT NULL
"))
death_tbl_global <- DBI::dbGetQuery(con_tmp, "SELECT fips, COUNT(*) AS deaths FROM base_all GROUP BY fips") %>%
    as_tibble() %>%
    mutate(fips = stringr::str_pad(as.character(fips), 5, pad = "0"),
           fips = apply_fips_patch(fips)) %>%
    filter(!is.na(fips) & fips %in% valid_fips_global)

# read provided mapping and restrict to valid_fips_global
merged_file <- here::here("data_raw","all_county_groupings.csv")
stopifnot(file.exists(merged_file))
merged_df <- readr::read_csv(merged_file, col_types = readr::cols(.default = "c"), show_col_types = FALSE)
cluster_col_name <- if ("CS Area Code" %in% names(merged_df)) "CS Area Code" else "County Set Name"

mapping_df_all <- merged_df %>%
    mutate(
        FIPS_raw = stringr::str_trim(as.character(FIPS)),
        FIPS_digits = ifelse(FIPS_raw == "" | is.na(FIPS_raw), NA_character_, gsub("[^0-9]", "", FIPS_raw)),
        FIPS5 = ifelse(is.na(FIPS_digits), NA_character_, stringr::str_pad(FIPS_digits, width = 5, side = "left", pad = "0")),
        cluster_id = stringr::str_trim(as.character(.data[[cluster_col_name]]))
    ) %>%
    filter(!is.na(FIPS5) & FIPS5 != "" & !is.na(cluster_id) & cluster_id != "") %>%
    distinct(FIPS5, cluster_id) %>%
    rename(fips = FIPS5, cluster = cluster_id) %>%
    mutate(fips = apply_fips_patch(fips))

mapping_global <- mapping_df_all %>% filter(fips %in% valid_fips_global)

# if any valid_fips_global are missing from mapping_global, add them as singletons
missing_in_mapping <- setdiff(valid_fips_global, mapping_global$fips)
if (length(missing_in_mapping) > 0) {
    message("Adding ", length(missing_in_mapping), " temporally-stable counties not present in all_county_groupings.csv as singleton clusters")
    add_df <- tibble(fips = missing_in_mapping, cluster = paste0("C_singleton_", missing_in_mapping))
    mapping_global <- bind_rows(mapping_global, add_df)
}

# initial table with aggregated deaths
clu_global <- tibble(fips = mapping_global$fips, cluster = mapping_global$cluster) %>%
    left_join(death_tbl_global %>% select(fips, deaths), by = "fips") %>%
    mutate(deaths = as.numeric(replace_na(deaths, 0)))

# greedy merge helpers (global)
compute_cluster_totals_global <- function(mapping_df_local) {
    mapping_df_local %>%
        filter(!is.na(cluster)) %>%
        group_by(cluster) %>%
        summarise(cluster_deaths = sum(as.numeric(deaths), na.rm = TRUE),
                  n_counties = n(), .groups = "drop") %>%
        arrange(cluster_deaths)
}
adjacent_clusters_global <- function(cluster_id, mapping_df_local) {
    fips_in_cluster_raw <- mapping_df_local$fips[mapping_df_local$cluster == cluster_id]
    fips_in_cluster <- intersect(fips_in_cluster_raw, names(nbrs_list_global))
    if (length(fips_in_cluster) == 0) return(character(0))
    neigh_fips <- unique(unlist(nbrs_list_global[fips_in_cluster], use.names = FALSE))
    neigh_fips <- intersect(neigh_fips, mapping_df_local$fips)
    neigh_clusters <- unique(mapping_df_local$cluster[mapping_df_local$fips %in% neigh_fips])
    neigh_clusters <- neigh_clusters[!is.na(neigh_clusters) & neigh_clusters != cluster_id]
    neigh_clusters
}
nearest_cluster_by_centroid_global <- function(cluster_id, mapping_df_local) {
    members_raw <- mapping_df_local$fips[mapping_df_local$cluster == cluster_id]
    members <- intersect(members_raw, rownames(centroid_mat_global))
    if (length(members) == 0) return(NA_character_)
    centroid_this <- tryCatch(colMeans(centroid_mat_global[members, , drop = FALSE]), error = function(e) NA_real_)
    if (any(is.na(centroid_this))) return(NA_character_)
    other_cluster_ids <- unique(mapping_df_local$cluster[!is.na(mapping_df_local$cluster) & mapping_df_local$cluster != cluster_id])
    if (length(other_cluster_ids) == 0) return(NA_character_)
    coords_list <- lapply(other_cluster_ids, function(cid) {
        memb <- intersect(mapping_df_local$fips[mapping_df_local$cluster == cid], rownames(centroid_mat_global))
        if (length(memb) == 0) return(NULL)
        as.numeric(colMeans(centroid_mat_global[memb, , drop = FALSE]))
    })
    valid_idx <- which(!vapply(coords_list, is.null, logical(1)))
    if (length(valid_idx) == 0) return(NA_character_)
    coords_mat <- do.call(rbind, coords_list[valid_idx])
    candidate_ids <- other_cluster_ids[valid_idx]
    dists <- sqrt((coords_mat[,1] - centroid_this[1])^2 + (coords_mat[,2] - centroid_this[2])^2)
    candidate_ids[which.min(dists)]
}

mapping_current_global <- clu_global %>% select(fips, cluster, deaths)
max_iters <- 10000L
for (iter in seq_len(max_iters)) {
    totals <- compute_cluster_totals_global(mapping_current_global %>% filter(!is.na(cluster)))
    small <- totals$cluster[which(totals$cluster_deaths < min_deaths)]
    if (length(small) == 0) break
    c_small <- small[1]
    neigh <- adjacent_clusters_global(c_small, mapping_current_global)
    if (length(neigh) > 0) {
        neigh_totals <- totals %>% filter(cluster %in% neigh)
        if (nrow(neigh_totals) == 0) {
            c_merge <- nearest_cluster_by_centroid_global(c_small, mapping_current_global)
        } else {
            c_merge <- neigh_totals$cluster[which.max(neigh_totals$cluster_deaths)]
        }
    } else {
        c_merge <- nearest_cluster_by_centroid_global(c_small, mapping_current_global)
    }
    if (is.na(c_merge) || c_merge == c_small) {
        largest_cluster <- totals$cluster[which.max(totals$cluster_deaths)]
        if (!is.na(largest_cluster) && largest_cluster != c_small) {
            c_merge <- largest_cluster
        } else {
            break
        }
    }
    mapping_current_global$cluster[mapping_current_global$cluster == c_small] <- c_merge
}
if (iter == max_iters) warning("Reached max iterations merging clusters in global mapping; some clusters may still be < min_deaths")

global_cluster_mapping <- mapping_current_global %>%
    group_by(cluster) %>% mutate(cluster_deaths = sum(as.numeric(deaths), na.rm = TRUE)) %>% ungroup() %>%
    select(fips, cluster, deaths, cluster_deaths)
readr::write_csv(global_cluster_mapping, file.path(diag_out_dir, "global_county_cluster_mapping_agg_deaths.csv.gz"))
message("Global cluster mapping built: ", n_distinct(global_cluster_mapping$cluster), " clusters; rows: ", nrow(global_cluster_mapping))

dbDisconnect(con_tmp, shutdown = TRUE)

# ----------------------- MAIN METRIC LOOP (uses global mapping) -----------------------
icd_root3 <- function(code){
    x <- toupper(gsub("[^A-Z0-9]","",as.character(code)))
    m <- stringr::str_match(x, "^([A-Z])([0-9])([0-9A-Z])")
    ifelse(is.na(m[,1]), NA_character_, paste0(m[,2],m[,3],m[,4]))
}
icd4_clean <- function(code){
    x <- toupper(gsub("[^A-Z0-9]","",as.character(code)))
    substr(x,1,4)
}

garbage_root_set <- readr::read_csv(garbage_csv, show_col_types = FALSE) %>%
    transmute(root3 = icd_root3(icd10)) %>% filter(!is.na(root3)) %>% distinct() %>% pull(root3)

compute_cstd_detail <- function(counts_df, std_ucr_df, us_within_df, maxH_global){
    if (nrow(counts_df) == 0) return(tibble(period=character(), cluster=character(), detail=numeric()))
    w_k <- counts_df %>%
        group_by(period, cluster, ucr39) %>%
        mutate(w = n / sum(n)) %>% ungroup() %>%
        select(period, cluster, ucr39, domain, w)
    
    base <- us_within_df %>% left_join(std_ucr_df, by=c("period","ucr39"))
    
    clusters_all <- distinct(counts_df, cluster) %>% pull(cluster)
    expanded <- base %>%
        left_join(w_k, by=c("period","ucr39","domain")) %>%
        group_by(period, ucr39, domain) %>%
        tidyr::complete(cluster = clusters_all, fill = list(w = NA_real_)) %>%
        ungroup() %>%
        mutate(w_eff = dplyr::coalesce(w, w_us))
    
    p_star <- expanded %>%
        group_by(period, cluster, domain) %>%
        summarise(p = sum(s * w_eff, na.rm=TRUE), .groups = "drop")
    
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

cluster_membership <- list()
metrics_list       <- list()
std_ucr_list       <- list()
us_within_root3    <- list()
us_within_icd4     <- list()
K_summary          <- list()

for (pname in names(periods)) {
    log_msg("Processing ", pname)
    yrs <- periods[[pname]]
    con <- dbConnect(duckdb::duckdb(), dbdir=":memory:")
    on.exit(try(dbDisconnect(con, shutdown=TRUE), silent=TRUE), add=TRUE)
    dbExecute(con, paste0("PRAGMA threads = ", threads))
    duckdb::duckdb_register(con, "garbage", tibble(root3 = garbage_root_set))
    
    yrs_range <- range(yrs)
    dbExecute(con, glue("
    CREATE TEMP TABLE base AS
    SELECT
      row_number() OVER () AS id,
      year,
      LPAD(CAST({`county_var`} AS VARCHAR), 5, '0') AS fips,
      ucod,
      {`UCR39_COL`} AS ucr39,
      record_1,record_2,record_3,record_4,record_5,
      record_6,record_7,record_8,record_9,record_10,
      record_11,record_12,record_13,record_14,record_15,
      record_16,record_17,record_18,record_19,record_20
    FROM parquet_scan('{normalizePath(parquet_dir, winslash = '/') }/*.parquet')
    WHERE year BETWEEN {yrs_range[1]} AND {yrs_range[2]}
      AND {`UCR39_COL`} IS NOT NULL
  "))
    
    dbExecute(con, "CREATE TEMP TABLE cert AS SELECT id, fips, ucod, ucr39 FROM base WHERE ucod IS NOT NULL")
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
    
    # use canonical year = min(years) object's centroid/nbrs (but mapping is global)
    sp_info <- period_spatial_list[[pname]]
    canonical_year <- as.character(min(yrs))
    if (!canonical_year %in% names(sp_info$year_objs)) {
        stop("Canonical year ", canonical_year, " not available in period_spatial_list for ", pname)
    }
    centroid_mat_use <- sp_info$year_objs[[canonical_year]]$centroid_mat
    nbrs_list_use <- sp_info$year_objs[[canonical_year]]$nbrs_list
    valid_county_geoids <- rownames(centroid_mat_use)
    
    # deaths per county for this period (used to compute per-period cluster metrics later)
    death_tbl <- dbGetQuery(con, "SELECT fips, COUNT(*) AS deaths FROM ucod_valid GROUP BY fips") %>%
        as_tibble() %>%
        mutate(fips = stringr::str_pad(as.character(fips), 5, pad = "0")) %>%
        mutate(fips = apply_fips_patch(fips)) %>%
        filter(!is.na(fips) & fips %in% valid_county_geoids)
    
    # GLOBAL mapping usage: use the pre-built 'global_cluster_mapping' but restrict to current period valid GEOIDs to avoid indexing errors
    mapping_period <- global_cluster_mapping %>%
        select(fips, cluster) %>%
        mutate(fips = as.character(fips)) %>%
        filter(fips %in% valid_county_geoids)
    
    readr::write_csv(mapping_period, file.path(diag_out_dir, paste0("mapping_for_period_", pname, ".csv")))
    dropped_period <- global_cluster_mapping %>% filter(!(fips %in% valid_county_geoids)) %>%
        distinct(fips, cluster) %>% mutate(reason = "not_stable_in_period")
    readr::write_csv(dropped_period, file.path(diag_out_dir, paste0("mapping_dropped_period_", pname, ".csv")))
    
    # build initial cluster mapping with deaths (period-specific deaths)
    clu_df <- death_tbl %>%
        left_join(mapping_period, by = c("fips" = "fips")) %>%
        select(fips, cluster, deaths) %>%
        mutate(cluster = as.character(cluster), fips = as.character(fips))
    
    cluster_membership[[pname]] <- clu_df %>% mutate(period = pname)
    
    # local helpers (period-scoped)
    valid_fips_set <- rownames(centroid_mat_use)
    nbrs_list <- nbrs_list_use
    centroid_mat <- centroid_mat_use
    
    compute_cluster_totals_period <- function(mapping_df_local) {
        mapping_df_local %>%
            filter(!is.na(cluster)) %>%
            group_by(cluster) %>%
            summarise(cluster_deaths = sum(as.numeric(deaths), na.rm = TRUE),
                      n_counties = n(), .groups = "drop") %>%
            arrange(cluster_deaths)
    }
    
    adjacent_clusters_period <- function(cluster_id, mapping_df_local) {
        fips_in_cluster_raw <- mapping_df_local$fips[mapping_df_local$cluster == cluster_id]
        fips_in_cluster <- intersect(fips_in_cluster_raw, names(nbrs_list))
        if (length(fips_in_cluster) == 0) return(character(0))
        neigh_fips <- unique(unlist(nbrs_list[fips_in_cluster], use.names = FALSE))
        neigh_fips <- intersect(neigh_fips, mapping_df_local$fips)
        neigh_clusters <- unique(mapping_df_local$cluster[mapping_df_local$fips %in% neigh_fips])
        neigh_clusters <- neigh_clusters[!is.na(neigh_clusters) & neigh_clusters != cluster_id]
        neigh_clusters
    }
    
    nearest_cluster_by_centroid_safe <- function(cluster_id, mapping_df_local) {
        members_raw <- mapping_df_local$fips[mapping_df_local$cluster == cluster_id]
        members <- intersect(members_raw, valid_fips_set)
        if (length(members) == 0) return(NA_character_)
        centroid_this <- tryCatch(colMeans(centroid_mat[members, , drop = FALSE]), error = function(e) NA_real_)
        if (any(is.na(centroid_this))) return(NA_character_)
        other_cluster_ids <- unique(mapping_df_local$cluster[!is.na(mapping_df_local$cluster) & mapping_df_local$cluster != cluster_id])
        if (length(other_cluster_ids) == 0) return(NA_character_)
        coords_list <- lapply(other_cluster_ids, function(cid) {
            memb <- intersect(mapping_df_local$fips[mapping_df_local$cluster == cid], valid_fips_set)
            if (length(memb) == 0) return(NULL)
            as.numeric(colMeans(centroid_mat[memb, , drop = FALSE]))
        })
        valid_idx <- which(!vapply(coords_list, is.null, logical(1)))
        if (length(valid_idx) == 0) return(NA_character_)
        coords_mat <- do.call(rbind, coords_list[valid_idx])
        candidate_ids <- other_cluster_ids[valid_idx]
        dists <- sqrt((coords_mat[,1] - centroid_this[1])^2 + (coords_mat[,2] - centroid_this[2])^2)
        candidate_ids[which.min(dists)]
    }
    
    mapping_current <- clu_df %>% select(fips, cluster, deaths)
    max_iters <- 10000L
    for (iter in seq_len(max_iters)) {
        totals <- compute_cluster_totals_period(mapping_current %>% filter(!is.na(cluster)))
        small <- totals$cluster[which(totals$cluster_deaths < min_deaths)]
        if (length(small) == 0) break
        c_small <- small[1]
        neigh <- adjacent_clusters_period(c_small, mapping_current)
        if (length(neigh) > 0) {
            neigh_totals <- totals %>% filter(cluster %in% neigh)
            if (nrow(neigh_totals) == 0) c_merge <- nearest_cluster_by_centroid_safe(c_small, mapping_current) else c_merge <- neigh_totals$cluster[which.max(neigh_totals$cluster_deaths)]
        } else {
            c_merge <- nearest_cluster_by_centroid_safe(c_small, mapping_current)
        }
        if (is.na(c_merge) || c_merge == c_small) break
        mapping_current$cluster[mapping_current$cluster == c_small] <- c_merge
    }
    if (iter == max_iters) warning("Reached max iterations merging clusters; some clusters may still be < min_deaths")
    
    clu_df_merged <- mapping_current %>% select(fips, cluster, deaths)
    cluster_membership[[pname]] <- clu_df_merged %>% mutate(period = pname)
    
    dbExecute(con, "DROP TABLE IF EXISTS map_fips_cluster")
    dbExecute(con, "CREATE TEMP TABLE map_fips_cluster (fips VARCHAR, cluster VARCHAR)")
    DBI::dbAppendTable(con, "map_fips_cluster", clu_df_merged %>% select(fips, cluster) %>% as.data.frame())
    
    std_ucr <- dbGetQuery(con, "
    SELECT ucr39, COUNT(*) AS n
    FROM ucod_valid
    GROUP BY 1
  ") %>% as_tibble() %>% mutate(period = pname, s = n / sum(n))
    std_ucr_list[[pname]] <- std_ucr %>% select(period, ucr39, s)
    
    us_r3 <- dbGetQuery(con, "
    SELECT ucr39, root3 AS domain, COUNT(*) AS n
    FROM ucod_valid
    GROUP BY 1,2
  ") %>% as_tibble() %>% group_by(ucr39) %>% mutate(w_us = n / sum(n)) %>% ungroup() %>% mutate(period = pname) %>% select(period, ucr39, domain, w_us)
    
    us_r4 <- dbGetQuery(con, "
    SELECT ucr39, icd4 AS domain, COUNT(*) AS n
    FROM ucod_valid
    GROUP BY 1,2
  ") %>% as_tibble() %>% group_by(ucr39) %>% mutate(w_us = n / sum(n)) %>% ungroup() %>% mutate(period = pname) %>% select(period, ucr39, domain, w_us)
    
    K_root3_global <- dplyr::n_distinct(us_r3$domain)
    K_icd4_global  <- dplyr::n_distinct(us_r4$domain)
    maxH_root3_global <- if (K_root3_global > 1) log(K_root3_global) else NA_real_
    maxH_icd4_global  <- if (K_icd4_global  > 1) log(K_icd4_global)  else NA_real_
    
    dbExecute(con, "DROP TABLE IF EXISTS mcod_codes_raw")
    dbExecute(con, "CREATE TEMP TABLE mcod_codes_raw (id BIGINT, fips VARCHAR, root3 VARCHAR)")
    dbExecute(con, "
    INSERT INTO mcod_codes_raw
    SELECT id, fips,
           SUBSTR(UPPER(regexp_replace(code,'[^A-Za-z0-9]','')),1,3) AS root3
    FROM (
      SELECT id, fips, UNNEST(LIST_VALUE(
        record_1,record_2,record_3,record_4,record_5,
        record_6,record_7,record_8,record_9,record_10,
        record_11,record_12,record_13,record_14,record_15,
        record_16,record_17,record_18,record_19,record_20
      )) AS code
      FROM base
    )
    WHERE code IS NOT NULL AND code <> '' AND LENGTH(code) >= 3
  ")
    
    dbExecute(con, "
    CREATE TEMP TABLE mcod_clean AS
    SELECT DISTINCT id, fips, root3
    FROM mcod_codes_raw
    WHERE root3 BETWEEN 'A00' AND 'Z99'
      AND SUBSTR(root3,2,1) BETWEEN '0' AND '9'
      AND root3 NOT IN (SELECT root3 FROM garbage)
  ")
    
    dbExecute(con, "
    CREATE TEMP TABLE mcod_nonunderlying AS
    SELECT m.id, m.fips, m.root3
    FROM mcod_clean m
    LEFT JOIN ucod_valid u
      ON m.id = u.id AND m.root3 = u.root3
    WHERE u.root3 IS NULL
  ")
    
    dbExecute(con, "
    CREATE TEMP TABLE mcod_with_cause AS
    SELECT n.id, n.fips, n.root3, u.ucr39
    FROM mcod_nonunderlying n
    JOIN ucod_valid u USING (id)
  ")
    
    mcod_counts <- DBI::dbGetQuery(con, "
    SELECT m.cluster, w.ucr39, w.root3 AS domain, COUNT(DISTINCT w.id) AS n
    FROM mcod_with_cause w
    JOIN map_fips_cluster m ON w.fips = m.fips
    WHERE m.cluster IS NOT NULL
    GROUP BY 1,2,3
  ") %>% as_tibble() %>% mutate(period = pname)
    
    us_mcod_r3 <- DBI::dbGetQuery(con, "
    SELECT ucr39, root3 AS domain, COUNT(DISTINCT id) AS n
    FROM mcod_with_cause
    GROUP BY 1,2
  ") %>% as_tibble() %>% group_by(ucr39) %>% mutate(w_us = n / sum(n)) %>% ungroup() %>% mutate(period = pname) %>% select(period, ucr39, domain, w_us)
    
    K_mcod_root3_global <- dplyr::n_distinct(us_mcod_r3$domain)
    maxH_mcod_root3_global <- if (K_mcod_root3_global > 1) log(K_mcod_root3_global) else NA_real_
    
    det_mcod_r3 <- compute_cstd_detail(mcod_counts %>% select(period, cluster, ucr39, domain, n),
                                       std_ucr %>% select(period, ucr39, s),
                                       us_mcod_r3,
                                       maxH_mcod_root3_global) %>% mutate(metric = "detail_mcod_root3_cstd")
    
    dbExecute(con, "
    CREATE TEMP TABLE mcod_nonunderlying_icd4 AS
    WITH expanded AS (
      SELECT
        nu.id,
        nu.fips,
        SUBSTR(UPPER(regexp_replace(code,'[^A-Za-z0-9]','')),1,4) AS icd4,
        SUBSTR(UPPER(regexp_replace(code,'[^A-Za-z0-9]','')),1,3) AS root3
      FROM mcod_nonunderlying nu
      JOIN (
        SELECT id, fips, UNNEST(LIST_VALUE(
          record_1,record_2,record_3,record_4,record_5,
          record_6,record_7,record_8,record_9,record_10,
          record_11,record_12,record_13,record_14,record_15,
          record_16,record_17,record_18,record_19,record_20
        )) AS code
        FROM base
      ) z USING (id,fips)
    )
    SELECT DISTINCT id, fips, icd4
    FROM expanded
    WHERE icd4 IS NOT NULL AND LENGTH(icd4) = 4
      AND root3 BETWEEN 'A00' AND 'Z99'
      AND SUBSTR(root3,2,1) BETWEEN '0' AND '9'
      AND root3 NOT IN (SELECT root3 FROM garbage)
  ")
    
    dbExecute(con, "
    CREATE TEMP TABLE mcod_with_cause_icd4 AS
    SELECT n.id, n.fips, n.icd4, u.ucr39
    FROM mcod_nonunderlying_icd4 n
    JOIN ucod_valid u USING (id)
  ")
    
    mcod4_counts <- DBI::dbGetQuery(con, "
    SELECT m.cluster, w.ucr39, w.icd4 AS domain, COUNT(DISTINCT w.id) AS n
    FROM mcod_with_cause_icd4 w
    JOIN map_fips_cluster m ON w.fips = m.fips
    WHERE m.cluster IS NOT NULL
    GROUP BY 1,2,3
  ") %>% as_tibble() %>% mutate(period = pname)
    
    us_mcod_icd4 <- DBI::dbGetQuery(con, "
    SELECT ucr39, icd4 AS domain, COUNT(DISTINCT id) AS n
    FROM mcod_with_cause_icd4
    GROUP BY 1,2
  ") %>% as_tibble() %>% group_by(ucr39) %>% mutate(w_us = n / sum(n)) %>% ungroup() %>% mutate(period = pname) %>% select(period, ucr39, domain, w_us)
    
    K_mcod_icd4_global <- dplyr::n_distinct(us_mcod_icd4$domain)
    maxH_mcod_icd4_global <- if (K_mcod_icd4_global > 1) log(K_mcod_icd4_global) else NA_real_
    
    det_mcod_icd4 <- compute_cstd_detail(mcod4_counts %>% select(period, cluster, ucr39, domain, n),
                                         std_ucr %>% select(period, ucr39, s),
                                         us_mcod_icd4,
                                         maxH_mcod_icd4_global) %>% mutate(metric = "detail_mcod_icd4_cstd")
    
    counts_r3 <- dbGetQuery(con, "
    SELECT m.cluster, v.ucr39, v.root3 AS domain, COUNT(*) AS n
    FROM ucod_valid v
    JOIN map_fips_cluster m ON v.fips = m.fips
    WHERE m.cluster IS NOT NULL
    GROUP BY 1,2,3
  ") %>% as_tibble() %>% mutate(period = pname)
    
    counts_r4 <- dbGetQuery(con, "
    SELECT m.cluster, v.ucr39, v.icd4 AS domain, COUNT(*) AS n
    FROM ucod_valid v
    JOIN map_fips_cluster m ON v.fips = m.fips
    WHERE m.cluster IS NOT NULL
    GROUP BY 1,2,3
  ") %>% as_tibble() %>% mutate(period = pname)
    
    det_r3 <- compute_cstd_detail(counts_r3 %>% select(period, cluster, ucr39, domain, n),
                                  std_ucr %>% select(period, ucr39, s),
                                  us_r3,
                                  maxH_root3_global) %>% mutate(metric = "detail_ucod_root3_cstd")
    
    det_r4 <- compute_cstd_detail(counts_r4 %>% select(period, cluster, ucr39, domain, n),
                                  std_ucr %>% select(period, ucr39, s),
                                  us_r4,
                                  maxH_icd4_global) %>% mutate(metric = "detail_ucod_icd4_cstd")
    
    metrics_list[[pname]] <- bind_rows(det_mcod_r3, det_mcod_icd4, det_r3, det_r4)
    us_within_root3[[pname]] <- us_r3
    us_within_icd4[[pname]]  <- us_r4
    K_summary[[pname]] <- tibble(
        period = pname,
        K_ucod_root3 = K_root3_global,
        K_ucod_icd4  = K_icd4_global,
        K_mcod_root3 = K_mcod_root3_global,
        K_mcod_icd4  = K_mcod_icd4_global
    )
    
    dbDisconnect(con, shutdown=TRUE); gc()
}

# ----------------------- outputs -----------------------
cluster_members_county <- bind_rows(cluster_membership) %>%
    transmute(period, cluster, fips, cluster_deaths = deaths)

cluster_members_cluster <- cluster_members_county %>%
    group_by(period, cluster) %>%
    summarise(
        cluster_deaths = sum(as.numeric(cluster_deaths), na.rm = TRUE),
        n_counties = n(),
        .groups = "drop"
    )

metrics_out <- bind_rows(metrics_list) %>%
    tidyr::pivot_wider(names_from = metric, values_from = detail) %>%
    arrange(period, cluster)

std_ucr_all      <- bind_rows(std_ucr_list) %>% arrange(period, ucr39)
us_within_r3_all <- bind_rows(us_within_root3) %>% arrange(period, ucr39, domain)
us_within_r4_all <- bind_rows(us_within_icd4)  %>% arrange(period, ucr39, domain)
K_summary_all    <- bind_rows(K_summary) %>% arrange(period)

write_csv(cluster_members_cluster, file.path(out_dir, "cluster_totals_by_period.csv.gz"))
write_csv(metrics_out, file.path(out_dir, "cluster_metrics_ucr39_cstd.csv.gz"))

cluster_members_all_periods <- bind_rows(cluster_membership) %>%
    transmute(period, cluster, fips, cluster_deaths = deaths)
write_csv(cluster_members_all_periods, file.path(out_dir, "county_cluster_membership.csv.gz"))

write_csv(std_ucr_all,      file.path(out_dir, "ucr39_standard_by_period.csv.gz"))
write_csv(us_within_r3_all, file.path(out_dir, "national_within_cause_shares_root3.csv.gz"))
write_csv(us_within_r4_all, file.path(out_dir, "national_within_cause_shares_icd4.csv.gz"))
write_csv(K_summary_all,    file.path(out_dir, "domain_K_counts_by_period.csv.gz"))

# console summary
if (nrow(metrics_out)) {
    cat("\n— Summary counts by period —\n")
    metrics_out %>%
        summarise(.by = period,
                  n_clusters = n_distinct(cluster),
                  nonNA_ucod_root3 = sum(!is.na(detail_ucod_root3_cstd)),
                  nonNA_ucod_icd4  = sum(!is.na(detail_ucod_icd4_cstd)),
                  nonNA_mcod_root3 = sum(!is.na(detail_mcod_root3_cstd)),
                  nonNA_mcod_icd4  = sum(!is.na(detail_mcod_icd4_cstd))) %>%
        print(n=Inf)
    cat("\n— K by period (for maxH = log K) —\n")
    print(K_summary_all, n=Inf)
}
cat("\nDone.\n")
