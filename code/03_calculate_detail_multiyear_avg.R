# ──────────────────────────────────────────────────────────────
# FAST county-cluster detail-richness metrics   ·  single-node
# ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(arrow);    library(dplyr);  library(stringr); library(purrr)
    library(sf);       library(tigris); library(igraph);  library(iNEXT)
    library(here);     library(readr)
})

log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "│", ..., "\n")

# 0 · Paths ----------------------------------------------------
parquet_dir <- if (fs::dir_exists(here("data_private","mcod"))) {
    here("data_private","mcod")
} else {
    here("data_private","mcod_sample")
}

out_csv       <- here("data","period_cluster_detail_d95.csv")
out_timingcsv <- here("data","period_cluster_timings.csv")
county_var    <- "county_ihme"        # column holding YYCCC county code

## 1 · IHME cross-walk  ----------------------------------------
obj  <- load(cw_path)          # loads and returns the object name
ihme_xwalk <- get(obj[1]) |>   # usually a data.frame with orig_fips / ihme_fips
    transmute(
        fips_raw   = str_pad(orig_fips , 5, pad = "0"),
        stable_fips = str_pad(ihme_fips, 5, pad = "0")
    )


# 2 · FIPS helpers (static) -----------------------------------
state_fips_lookup <- tibble::tribble(
    ~state_abbr, ~state_fips,
    "AL", "01", "AK", "02", "AZ", "04", "AR", "05", "CA", "06", "CO", "08", "CT", "09",
    "DE", "10", "DC", "11", "FL", "12", "GA", "13", "HI", "15", "ID", "16", "IL", "17",
    "IN", "18", "IA", "19", "KS", "20", "KY", "21", "LA", "22", "ME", "23", "MD", "24",
    "MA", "25", "MI", "26", "MN", "27", "MS", "28", "MO", "29", "MT", "30", "NE", "31",
    "NV", "32", "NH", "33", "NJ", "34", "NM", "35", "NY", "36", "NC", "37", "ND", "38",
    "OH", "39", "OK", "40", "OR", "41", "PA", "42", "RI", "44", "SC", "45", "SD", "46",
    "TN", "47", "TX", "48", "UT", "49", "VT", "50", "VA", "51", "WA", "53", "WV", "54",
    "WI", "55", "WY", "56", "PR", "72"
)

nchs_to_fips <- tibble::tribble(
    ~nchs_code, ~fips_code,
    "01","01","02","02","03","04","04","05","05","06","06","08","07","09",
    "08","10","09","11","10","12","11","13","12","15","13","16","14","17",
    "15","18","16","19","17","20","18","21","19","22","20","23","21","24",
    "22","25","23","26","24","27","25","28","26","29","27","30","28","31",
    "29","32","30","33","31","34","32","35","33","36","34","37","35","38",
    "36","39","37","40","38","41","39","42","40","44","41","45","42","46",
    "43","47","44","48","45","49","46","50","47","51","48","53","49","54",
    "50","55","51","56","52","72"
)


# 3 · Adjacency graph (once) -------------------------------
options(tigris_use_cache = TRUE, tigris_class = "sf")
counties_sf <- tigris::counties(cb = TRUE, year = 2020) |> select(GEOID)
adj <- sf::st_touches(counties_sf)
edge_df <- tibble(
    from = rep(counties_sf$GEOID, lengths(adj)),
    to   = counties_sf$GEOID[ unlist(adj) ]
) |> filter(from < to)
g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = counties_sf$GEOID)

# 4 · Richness helper  (patched) -------------------------------
REF_SIZE <- 358
TARGET_COVERAGE <- 0.95

detail_metric <- function(vec) {
    vec <- vec[!is.na(vec)]
    if (length(vec) == 0) return(c(d95 = NA_real_, cov = 0))
    
    ftab <- as.numeric(table(vec))
    N    <- sum(ftab)
    f1   <- sum(ftab == 1)
    C_hat <- 1 - f1 / N                        # observed sample coverage
    
    est <- tryCatch(
        iNEXT::estimateD(ftab, "abundance",
                         base = "coverage", level = TARGET_COVERAGE),
        error = function(e) NULL)
    
    # work with either old ("Estimator") or new ("qD") column name
    d_raw <- NA_real_
    if (!is.null(est)) {
        if ("Estimator" %in% names(est))
            d_raw <- est$Estimator[ est$Order.q == 0 ]
        else if ("qD" %in% names(est))
            d_raw <- est$qD[ est$Order.q == 0 ]
    }
    
    # fall-back: use observed richness if extrapolation failed
    if (is.na(d_raw)) d_raw <- length(ftab)
    
    c(d95 = d_raw / REF_SIZE,
      cov = C_hat)
}
# 5 · Periods --------------------------------------------------
periods <- list(
    `1999-2004` = 1999:2004,
    `2005-2010` = 2005:2010,
    `2011-2017` = 2011:2017,
    `2018-2022` = 2018:2022
)

# 6 · Fast clustering (igraph) --------------------------------
## 6 · Fast clustering with igraph (robust) --------------------
cluster_counties <- function(death_tbl, label) {
    # deaths named vector for O(1) lookup
    deaths <- setNames(death_tbl$n, death_tbl$fips)
    
    # initialise cluster id map
    ids <- setNames(rep(NA_character_, length(deaths)), names(deaths))
    
    # seeds: counties sorted by descending deaths
    seeds <- names(sort(deaths, decreasing = TRUE))
    cid   <- 1
    
    for (seed in seeds) {
        if (!is.na(ids[seed])) next                # already clustered
        if (!seed %in% V(g)$name)   next           # not in 2020 graph
        
        cluster <- seed
        total   <- deaths[seed]
        
        repeat {
            # union of neighbours of current cluster nodes
            nbrs <- unique(unlist(lapply(cluster, \(v)
                                         neighbors(g, v, mode = "all") |> names()
            )))
            
            # keep only new nodes with death info
            nbrs <- nbrs[is.na(ids[nbrs]) & !is.na(deaths[nbrs])]
            nbrs <- setdiff(nbrs, cluster)
            
            if (!length(nbrs) || total >= 50) break
            
            # pick neighbour with most deaths
            best   <- nbrs[ which.max(deaths[nbrs]) ]
            cluster <- c(cluster, best)
            total   <- total + deaths[best]
        }
        
        ids[cluster] <- paste0(label, "_C", cid)
        cid <- cid + 1
    }
    
    ids
}

# 7 · Main loop -----------------------------------------------
all_metrics <- list(); all_time <- list()

for (nm in names(periods)) {
    yrs <- periods[[nm]]
    tic <- proc.time()[3]
    
    ## 7A · Load
    ds <- open_dataset(parquet_dir)
    cert <- ds |>
        filter(year %in% yrs, !is.na(ucr358)) |>
        select(countyrs_chr = !!county_var, ucr358) |>
        collect()
    load_t <- proc.time()[3] - tic
    log_msg(nm, "load", round(load_t,1), "s")
    
    ## 7B · FIPS parse
    tic <- proc.time()[3]
    cert <- cert |>
        mutate(
            state_code  = substr(countyrs_chr, 1, 2),
            county_code = substr(countyrs_chr, 3, 5),
            pre2003_fips = paste0(
                nchs_to_fips$fips_code[ match(state_code, nchs_to_fips$nchs_code) ],
                county_code
            ),
            state_abbr = substr(countyrs_chr, 1, 2),
            county_part = substr(countyrs_chr, 3, 5)
        ) |>
        left_join(state_fips_lookup, by = "state_abbr") |>
        mutate(
            fips_raw = if_else(substr(countyrs_chr,1,1) == "0", pre2003_fips,
                               paste0(state_fips, county_part)),
            fips_raw = str_pad(fips_raw,5,pad="0")
        ) |>
        select(fips_raw, ucr358) |>
        left_join(ihme_xwalk, by = "fips_raw") |>
        transmute(fips = coalesce(stable_fips, fips_raw), ucr358)
    fips_t <- proc.time()[3] - tic
    log_msg(nm, "FIPS", round(fips_t,1), "s")
    
    ## 7C · Cluster
    tic <- proc.time()[3]
    death_tbl <- count(cert, fips, name = "n")
    clusters  <- cluster_counties(death_tbl, nm)
    clust_t <- proc.time()[3] - tic
    log_msg(nm, "cluster", round(clust_t,1), "s")
    
    ## 7D · Richness
    tic <- proc.time()[3]
    metrics <- cert |>
        mutate(cluster = clusters[fips]) |>
        group_by(cluster) |>
        summarise(
            met = list(detail_metric(ucr358)),
            n  = n(), .groups = "drop"
        ) |>
        mutate(
            detail_d95 = map_dbl(met, "d95"),
            detail_cov = map_dbl(met, "cov"),
            period     = nm
        ) |>
        select(cluster, period, n, detail_d95, detail_cov)
    rich_t <- proc.time()[3] - tic
    log_msg(nm, "richness", round(rich_t,1), "s")
    
    all_metrics[[nm]] <- metrics
    all_time[[nm]] <- tibble(
        period = nm,
        stage  = c("load","FIPS","cluster","richness"),
        secs   = c(load_t, fips_t, clust_t, rich_t)
    )
    log_msg(nm, "DONE ⏱ total", round(sum(load_t,fips_t,clust_t,rich_t),1), "s\n")
}

# 8 · Save -----------------------------------------------------
period_detail <- bind_rows(all_metrics)
timings_df    <- bind_rows(all_time)

write_csv(period_detail, out_csv)
write_csv(timings_df,   out_timingcsv)

log_msg("Finished. Results:", out_csv)
log_msg("Timings :", out_timingcsv)

