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
    library(arrow);    library(dplyr);  library(stringr); library(purrr)
    library(readr);    library(tidyr);   library(tibble)
    library(sf);       library(tigris);  library(igraph);  library(here)
    library(fs);       library(iNEXT)
})
options(tigris_use_cache = TRUE, tigris_class = "sf")

# initial stuff and helpers
out_dir     <- here("output"); dir_create(out_dir) 
parquet_dir <- if (dir_exists(here("data_private","mcod"))) here("data_private","mcod") else here("data_private","mcod_sample")
garbage_csv <- here("data_raw","cause-codes","gbd_garbage_codes_without_overdose.csv")

county_var  <- "county_ihme"
min_deaths  <- 2000L
periods <- list(
    "1999_2006" = 1999:2006,
    "2007_2014" = 2007:2014,
    "2015_2022" = 2015:2022
)

M_REF        <- 2000L       # fixed reference size (pure interpolation)
EXTRAP_MULT  <- 2           # upper bound used only for diagnostics
crs_proj     <- 5070        # US Albers (meters) for distance ops

icd_root3 <- function(code) {
    x <- toupper(gsub("[^A-Z0-9]", "", as.character(code)))
    m <- stringr::str_match(x, "^([A-Z])([0-9])([0-9A-Z])")
    ifelse(is.na(m[,1]), NA_character_, paste0(m[,2], m[,3], m[,4]))
}

extract_q1 <- function(est) {
    if (inherits(est, "try-error") || is.null(est)) return(NA_real_)
    col <- intersect(c("Estimator", "qD"), names(est)); if (!length(col)) return(NA_real_)
    oq  <- suppressWarnings(as.numeric(est$Order.q))
    val <- suppressWarnings(as.numeric(est[[col[1]]]))
    out <- val[which(oq == 1)]
    ifelse(length(out) && is.finite(out) && out > 0, out, NA_real_)
}

compute_detail_row <- function(tbl_counts, max_H_period, m_ref = M_REF, extrap_mult = EXTRAP_MULT) {
    ftab <- as.numeric(tbl_counts$n); names(ftab) <- tbl_counts$domain_code
    ftab <- ftab[is.finite(ftab) & ftab > 0]
    n <- sum(ftab); S <- length(ftab)
    f1 <- sum(ftab == 1); f2 <- sum(ftab == 2)
    
    if (!is.finite(n) || n < 2 || S < 2 || !is.finite(max_H_period) || max_H_period <= 0) {
        return(tibble(
            deaths_no_garbage = n, S = S, f1 = f1, f2 = f2,
            H_raw = NA_real_, detail_phillips_raw = NA_real_,
            detail_phillips_refsize = NA_real_,
            D1_mref = NA_real_, D1_eff_used = NA_real_,
            m_ref = m_ref, m_max = NA_real_, eff_m_used = NA_real_,
            used_extrapolation = NA, fallback_reason = "degenerate_cluster_or_period"
        ))
    }
    
    p <- ftab / n
    H_raw <- -sum(p * log(p))
    detail_raw <- 100 * H_raw / max_H_period
    
    # size-standardized detail at m_ref
    m_max <- floor(extrap_mult * n)
    m_use <- max(2L, min(as.integer(m_ref), m_max))   # should equal 2000 for all clusters
    
    est <- try(iNEXT::estimateD(ftab, datatype = "abundance", base = "size",
                                level = m_use, conf = FALSE), silent = TRUE)
    D1_use <- extract_q1(est)
    fb <- "ok"
    if (!is.finite(D1_use)) {
        fb <- "estimateD_failed_expH"
        D1_use <- if (is.finite(H_raw)) exp(H_raw) else NA_real_
    }
    
    detail_refsize <- if (is.finite(D1_use) && D1_use > 0) 100 * log(D1_use) / max_H_period else NA_real_
    
    D1_at_mref <- NA_real_
    if (m_ref <= m_max) {
        est_ref <- try(iNEXT::estimateD(ftab, "abundance", base = "size",
                                        level = as.integer(m_ref), conf = FALSE), silent = TRUE)
        D1_at_mref <- extract_q1(est_ref)
    }
    
    tibble(
        deaths_no_garbage = n, S = S, f1 = f1, f2 = f2,
        H_raw = H_raw, detail_phillips_raw = detail_raw,
        detail_phillips_refsize = detail_refsize,
        D1_mref = D1_at_mref, D1_eff_used = D1_use,
        m_ref = as.integer(m_ref), m_max = m_max, eff_m_used = m_use,
        used_extrapolation = (m_use > n),
        fallback_reason = fb
    )
}

# clustering helpers

# greedy clustering starts with largest-deaths counties and expands toward min_deaths via adjacency
make_greedy_clusters <- function(death_tbl, g, min_deaths) {
    deaths <- setNames(death_tbl$deaths, death_tbl$fips)
    to_assign <- names(deaths)
    clusters <- setNames(rep(NA_character_, length(deaths)), to_assign)
    cid <- 1L
    while (length(to_assign) > 0) {
        this <- to_assign[which.max(deaths[to_assign])]
        if (deaths[this] >= min_deaths) {
            clusters[this] <- paste0("C", cid)
            to_assign <- setdiff(to_assign, this)
            cid <- cid + 1L; next
        }
        cluster <- this; total <- deaths[this]; avail <- setdiff(to_assign, this)
        repeat {
            nbrs <- unique(unlist(lapply(cluster, function(v) {
                nv <- tryCatch(neighbors(g, v, mode="all"), error=function(e) NULL)
                if (is.null(nv)) character(0) else names(nv)
            })))
            nbrs <- setdiff(intersect(nbrs, avail), cluster)
            if (!length(nbrs)) break
            best <- nbrs[which.min(abs((total + deaths[nbrs]) - min_deaths))]
            new_tot <- total + deaths[best]
            if (is.na(new_tot)) break
            if (new_tot <= min_deaths * 1.5 || total < min_deaths) {
                cluster <- c(cluster, best); total <- new_tot; avail <- setdiff(avail, best)
            } else break
        }
        clusters[cluster] <- paste0("C", cid)
        to_assign <- setdiff(to_assign, cluster)
        cid <- cid + 1L
    }
    clusters
}

# merge any undersized clusters to adjacent “best” cluster by deaths
merge_small_clusters <- function(clu, death_tbl, g, min_deaths, max_iter = 20L) {
    deaths_vec <- setNames(death_tbl$deaths, death_tbl$fips)
    unlabeled <- names(deaths_vec)[is.na(clu[names(deaths_vec)])]
    if (length(unlabeled)) {
        next_id <- suppressWarnings(max(as.integer(sub("^C","", na.omit(unique(clu)))), na.rm = TRUE))
        if (!is.finite(next_id)) next_id <- 0L
        new_ids <- paste0("C", seq(next_id + 1L, next_id + length(unlabeled)))
        clu[unlabeled] <- new_ids
    }
    safe_top <- function(x, exclude = character(0)) {
        x <- x[setdiff(names(x), exclude)]; x <- x[is.finite(x)]
        if (!length(x)) return(NA_character_)
        names(x)[which.max(x)]
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
            neigh_vertices <- unique(unlist(lapply(members, function(v) {
                tryCatch(names(neighbors(g, v, mode="all")), error=function(e) character(0))
            })))
            neigh_clusters <- setdiff(unique(clu[intersect(neigh_vertices, names(clu))]), sc)
            neigh_clusters <- neigh_clusters[neigh_clusters %in% names(sums)]
            target <- if (length(neigh_clusters)) neigh_clusters[which.max(sums[neigh_clusters])]
            else safe_top(sums, exclude = sc)
            if (!nzchar(target) || is.na(target)) next
            clu[members] <- target
        }
    }
    clu
}

# if there are any NA counties
fix_na_by_nearest <- function(clu, counties_sf) {
    cty <- counties_sf %>% select(GEOID, geometry)
    cty$cluster <- unname(clu[cty$GEOID])
    na_idx <- which(is.na(cty$cluster)); ok_idx <- which(!is.na(cty$cluster))
    if (!length(na_idx)) return(clu)
    cty_proj <- st_transform(cty, crs_proj)
    nn <- st_nearest_feature(cty_proj[na_idx, ], cty_proj[ok_idx, ])
    clu[cty$GEOID[na_idx]] <- cty$cluster[ok_idx[nn]]
    clu
}

# if any clusters remain < min_deaths put them in the nearest cluster by centroid distance
# this exists because one clusters doesn't have anything adjacent so it's a fall back
merge_small_clusters_by_distance <- function(clu, death_tbl, counties_sf, min_deaths) {
    deaths_vec <- setNames(death_tbl$deaths, death_tbl$fips)
    keep <- intersect(names(clu), counties_sf$GEOID)
    cty <- counties_sf %>% select(GEOID, geometry)
    cty$cluster <- unname(clu[cty$GEOID])
    
    sums <- tapply(deaths_vec[keep], clu[keep], sum, na.rm = TRUE)
    sums <- sums[is.finite(sums)]
    small <- names(sums)[sums < min_deaths]
    big   <- names(sums)[sums >= min_deaths]
    if (!length(small) || !length(big)) return(clu)
    
    cty_proj <- st_transform(cty, crs_proj)
    clu_geom <- cty_proj %>%
        filter(!is.na(cluster)) %>%
        group_by(cluster) %>% summarise(geometry = st_union(geometry), .groups = "drop")
    clu_cent <- st_centroid(clu_geom)
    small_sf <- clu_cent[clu_cent$cluster %in% small, ]
    big_sf   <- clu_cent[clu_cent$cluster %in% big, ]
    if (!nrow(small_sf) || !nrow(big_sf)) return(clu)
    
    nn <- st_nearest_feature(small_sf, big_sf)
    targets <- big_sf$cluster[nn]
    
    for (i in seq_len(nrow(small_sf))) {
        sc <- small_sf$cluster[i]; tg <- targets[i]
        members <- names(clu)[clu == sc]; if (!length(members)) next
        clu[members] <- tg
    }
    clu
}

# lookup tables
stopifnot(file.exists(garbage_csv))
lookup_garbage_root <- read_csv(garbage_csv, show_col_types = FALSE) %>%
    transmute(icd10_root = icd_root3(icd10)) %>%
    filter(!is.na(icd10_root)) %>% distinct()

counties_sf <- tigris::counties(cb = TRUE, year = 2020) %>% select(GEOID, geometry)
adj <- st_touches(counties_sf)
edge_df <- tibble(from = rep(counties_sf$GEOID, lengths(adj)),
                  to   = counties_sf$GEOID[unlist(adj)]) %>% filter(from < to)
g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = counties_sf$GEOID)


# main script
ds <- open_dataset(parquet_dir)
cluster_results <- list()
cluster_membership <- list()

for (pname in names(periods)) {
    message("Processing ", pname)
    
    cert <- ds %>%
        filter(year %in% periods[[pname]], !is.na(ucod)) %>%
        select(!!county_var, ucod) %>% collect() %>%
        mutate(fips = stringr::str_pad(as.character(.data[[county_var]]), 5, pad = "0")) %>%
        filter(fips %in% counties_sf$GEOID) %>%
        mutate(icd10_root = icd_root3(ucod))
    
    cert_valid <- anti_join(cert, lookup_garbage_root, by = "icd10_root") %>%
        mutate(domain_code = icd10_root)
    
    K_period <- n_distinct(cert_valid$domain_code)
    if (K_period <= 1) { warning("Period ", pname, " has <=1 domain; skipping."); next }
    max_H_period <- log(K_period)

    death_tbl <- count(cert_valid, fips, name = "deaths")
    if (nrow(death_tbl) == 0) { next }
    
    # build clusters and put any unmatched county into a cluster
    clu <- make_greedy_clusters(death_tbl, g, min_deaths)
    clu <- merge_small_clusters(clu, death_tbl, g, min_deaths)
    clu <- fix_na_by_nearest(clu, counties_sf)
    
    deaths_vec <- setNames(death_tbl$deaths, death_tbl$fips)
    sums_post <- tapply(deaths_vec, clu[names(deaths_vec)], sum, na.rm = TRUE)
    if (any(is.finite(sums_post) & sums_post < min_deaths)) {
        clu <- merge_small_clusters_by_distance(clu, death_tbl, counties_sf, min_deaths)
        clu <- fix_na_by_nearest(clu, counties_sf)
        sums_post <- tapply(deaths_vec, clu[names(deaths_vec)], sum, na.rm = TRUE)
    }
    
    clu_df <- tibble(fips = names(clu), cluster = unname(clu)) %>%
        left_join(death_tbl, by = "fips")
    cluster_membership[[pname]] <- clu_df %>% mutate(period = pname)
    
    counts_tbl <- cert_valid %>%
        mutate(cluster = unname(clu[fips])) %>%
        filter(!is.na(cluster)) %>%
        count(cluster, domain_code, name = "n")
    
    res <- counts_tbl %>% group_by(cluster) %>%
        group_modify(~ compute_detail_row(.x, max_H_period, m_ref = M_REF,
                                          extrap_mult = EXTRAP_MULT)) %>%
        ungroup() %>%
        mutate(period = pname, unit = "cluster", unit_id = cluster)
    
    cluster_results[[pname]] <- res
    }

# output as csv file
cluster_out <- bind_rows(cluster_results)
write_csv(cluster_out, file.path(out_dir, "cluster_phillips_detail.csv"))
write_csv(bind_rows(cluster_membership), file.path(out_dir, "county_cluster_membership.csv"))
message("Script finished running.")