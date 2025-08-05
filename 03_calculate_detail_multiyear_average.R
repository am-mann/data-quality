# ──────────────────────────────────────────────────────────────
# County-cluster detail-richness metrics · parallel + timings
# ──────────────────────────────────────────────────────────────

library(future); library(furrr); library(tictoc)
library(dplyr); library(purrr); library(tibble); library(stringr)
library(tidyr); library(sf); library(tigris); library(arrow)
library(iNEXT); library(here); library(readr)

future::plan(multisession, workers = future::availableCores() - 1)
options(tigris_use_cache = TRUE, tigris_class = "sf")

# ── 0 · Paths and helpers ─────────────────────────────────────
parquet_dir  <- if (fs::dir_exists(here("data_private","mcod"))) {
    here("data_private","mcod")
} else {
    here("data_private","mcod_sample")
}
dictionary_dir <- here("data_raw","cause-codes")
cw_path        <- here("data_raw","ihme_fips.rda")
out_csv        <- here("data","period_cluster_detail_d95.csv")
county_var     <- "countyrs"

source(here("code/helpers","dq_entropy_helper.R"))

cw_name   <- load(cw_path)
ihme_xwalk <- get(cw_name[1]) |>
    mutate(
        fips        = stringr::str_pad(as.character(orig_fips), 5, pad = "0"),
        stable_fips = stringr::str_pad(as.character(ihme_fips), 5, pad = "0")
    ) |>
    select(fips, stable_fips) |>
    distinct()

state_fips_lookup <- tibble::tribble(
    ~state_abbr, ~state_fips,
    "AL","01","AK","02","AZ","04","AR","05","CA","06","CO","08","CT","09",
    "DE","10","DC","11","FL","12","GA","13","HI","15","ID","16","IL","17",
    "IN","18","IA","19","KS","20","KY","21","LA","22","ME","23","MD","24",
    "MA","25","MI","26","MN","27","MS","28","MO","29","MT","30","NE","31",
    "NV","32","NH","33","NJ","34","NM","35","NY","36","NC","37","ND","38",
    "OH","39","OK","40","OR","41","PA","42","RI","44","SC","45","SD","46",
    "TN","47","TX","48","UT","49","VT","50","VA","51","WA","53","WV","54",
    "WI","55","WY","56","PR","72"
)

nchs_to_fips_state_lookup <- tibble::tribble(
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

# ── 1 · Adjacency helper ──────────────────────────────────────
counties_sf <- tigris::counties(cb = TRUE, year = 2020) |> select(GEOID)
adj_list    <- sf::st_touches(counties_sf)
neigh <- function(fips) {
    idx <- match(fips, counties_sf$GEOID)
    if (is.na(idx)) character(0) else counties_sf$GEOID[ adj_list[[idx]] ]
}

# ── 2 · Richness metric ───────────────────────────────────────
REF_SIZE <- 358
TARGET_COVERAGE <- 0.95
get_detail_metric <- function(vec,
                              ref_size = REF_SIZE,
                              coverage = TARGET_COVERAGE) {
    vec <- vec[!is.na(vec)]
    if (!length(vec)) return(list(raw = NA_real_, d95 = NA_real_, cov = 0, n = 0))
    ftab <- as.numeric(table(vec))
    N <- sum(ftab); f1 <- sum(ftab == 1); C_hat <- 1 - f1 / N
    est <- tryCatch(
        iNEXT::estimateD(ftab, datatype = "abundance",
                         base = "coverage", level = coverage),
        error = function(e) NULL
    )
    D95 <- if (is.null(est) || !"Estimator" %in% names(est)) NA_real_ else
        as.numeric(est$Estimator[ est$Order.q == 0 ])
    list(raw = length(ftab) / ref_size,
         d95 = D95 / ref_size,
         cov = C_hat,
         n = N)
}

# ── 3 · Periods and loader ────────────────────────────────────
periods <- list(`1999-2004` = 1999:2004,
                `2005-2010` = 2005:2010,
                `2011-2017` = 2011:2017,
                `2018-2022` = 2018:2022)

load_period <- function(yrs, path = parquet_dir) {
    arrow::open_dataset(path) %>%
        filter(year %in% yrs, !is.na(ucr358)) %>%
        select(countyrs_chr = !!county_var, ucr358) %>%
        collect()
}

# ── 4 · Timer wrapper ─────────────────────────────────────────
time_it <- function(expr) {
    t0 <- tic()
    val <- force(expr)
    dt <- toc(quiet = TRUE)
    list(result = val, secs = dt$toc - dt$tic)
}

# ── 5 · Main worker function ─────────────────────────────────
run_period <- function(yrs, per_name) {
    
    ## loader
    loader <- time_it( load_period(yrs) )
    cert_raw <- loader$result
    
    ## FIPS
    fips <- time_it({
        cert_raw |>
            mutate(
                state_code  = substr(countyrs_chr, 1, 2),
                county_code = substr(countyrs_chr, 3, 5),
                pre2003_fips = paste0(
                    nchs_to_fips_state_lookup$fips_code[
                        match(state_code, nchs_to_fips_state_lookup$nchs_code)
                    ],
                    county_code
                ),
                state_abbr  = substr(countyrs_chr, 1, 2),
                county_part = substr(countyrs_chr, 3, 5)
            ) |>
            left_join(state_fips_lookup, by = "state_abbr") |>
            mutate(
                fips_raw = if_else(substr(countyrs_chr,1,1) == "0", pre2003_fips,
                                   paste0(state_fips, county_part)),
                fips_raw = stringr::str_pad(fips_raw, 5, pad = "0")
            ) |>
            select(fips_raw, ucr358) |>
            left_join(ihme_xwalk, by = c("fips_raw" = "fips")) |>
            mutate(fips = coalesce(stable_fips, fips_raw)) |>
            select(fips, ucr358)
    })
    cert_raw <- fips$result
    
    ## cluster
    cluster <- time_it({
        death_counts <- count(cert_raw, fips, name = "n_deaths")
        cluster_id <- setNames(rep(NA_character_, nrow(death_counts)),
                               death_counts$fips)
        next_id <- 1
        for (seed in death_counts$fips) {
            if (!is.na(cluster_id[[seed]])) next
            group <- c(seed)
            total <- death_counts$n_deaths[match(seed, death_counts$fips)]
            while (total < 50) {
                cand  <- unique(unlist(lapply(group, neigh)))
                cand  <- setdiff(cand, group)
                cand  <- cand[is.na(cluster_id[cand])]
                if (!length(cand)) break
                best  <- cand[which.max(
                    death_counts$n_deaths[match(cand, death_counts$fips)]
                )]
                group <- c(group, best)
                total <- sum(
                    death_counts$n_deaths[match(group, death_counts$fips)], na.rm = TRUE
                )
            }
            cluster_id[group] <- paste0(per_name, "_C", next_id)
            next_id <- next_id + 1
        }
        cluster_id
    })
    
    ## richness
    richness <- time_it({
        cert_raw |>
            mutate(cluster = cluster$result[fips]) |>
            group_by(cluster) |>
            summarise(
                tmp      = list(get_detail_metric(ucr358)),
                n_deaths = n(),
                .groups  = "drop"
            ) |>
            mutate(
                detail_d95 = purrr::map_dbl(tmp, "d95"),
                detail_cov = purrr::map_dbl(tmp, "cov"),
                period     = per_name
            ) |>
            select(cluster, period, n_deaths, detail_d95, detail_cov)
    })
    
    list(
        data = richness$result,
        timings = tibble(
            period  = per_name,
            stage   = c("loader", "FIPS", "cluster", "richness"),
            seconds = c(loader$secs, fips$secs, cluster$secs, richness$secs)
        )
    )
}

# ── 6 · Run all periods ───────────────────────────────────────
tic("Full run")

results <- future_map2(periods, names(periods), run_period,
                       .options = furrr::furrr_options(seed = TRUE))

period_detail <- bind_rows(map(results, "data"))
timings       <- bind_rows(map(results, "timings"))

toc()

# ── 7 · Save results ──────────────────────────────────────────
readr::write_csv(period_detail, out_csv)
print(head(period_detail, 10))
print(timings)

