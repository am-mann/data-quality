library(arrow)
library(dplyr)
library(stringr)
library(purrr)
library(readr)
library(readxl)
library(tidyr)
library(tibble)
library(sf)
library(tigris)
library(igraph)
library(here)
library(fs)

options(tigris_use_cache = TRUE, tigris_class = "sf")

# preliminary stuff
out_dir <- here("output"); dir_create(out_dir)
parquet_dir <- if (dir_exists(here("data_private","mcod"))) {
    here("data_private","mcod")
} else {
    here("data_private","mcod_sample")
}
gbd_icd_excel <- here("data_raw","IHME_GBD_2021_COD_CAUSE_ICD_CODE_MAP.XLSX")
gbd_hierarchy_xlsx <- here("data_raw","GBD_2021_CAUSE_HIERARCHY_Y2024M05D16.XLSX")
county_var <- "county_ihme"
min_deaths <- 200
periods <- list(
    "1999_2004" = 1999:2004,
    "2005_2010" = 2005:2010,
    "2011_2017" = 2011:2017,
    "2018_2022" = 2018:2022
)

# finds adjacent counties to merge into clusters for counties that are too small to be reliable
counties_sf <- counties(cb = TRUE, year = 2020) %>% select(GEOID, geometry)
adj <- st_touches(counties_sf)
edge_df <- tibble(from = rep(counties_sf$GEOID, lengths(adj)),
                  to = counties_sf$GEOID[unlist(adj)]) %>% filter(from < to)
g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = counties_sf$GEOID)

# Cluster builders
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
            nbrs <- intersect(setdiff(unique(unlist(lapply(cluster, function(v)
                neighbors(g, v, mode="all") %>% names()))), cluster), avail)
            if (length(nbrs) == 0) break
            best <- nbrs[which.min(abs((total + deaths[nbrs]) - min_deaths))]
            new_tot <- total + deaths[best]
            if (new_tot <= min_deaths * 1.5 || total < min_deaths) {
                cluster <- c(cluster, best)
                total <- new_tot
                avail <- setdiff(avail, best)
            } else break
        }
        clusters[cluster] <- paste0("C", cid)
        to_assign <- setdiff(to_assign, cluster)
        cid <- cid + 1L
    }
    clusters
}
merge_small_clusters <- function(clu, death_tbl, g, min_deaths = 200) {
    deaths_vec <- setNames(death_tbl$deaths, death_tbl$fips)
    repeat {
        sums <- tapply(deaths_vec, clu, sum, na.rm = TRUE)
        small <- names(sums)[sums < min_deaths]
        if (length(small) == 0) break
        for (sc in small) {
            members <- names(clu)[clu == sc]
            neigh_clusters <- unique(clu[neighbors(g, members, mode = "all") %>% names()])
            neigh_clusters <- setdiff(neigh_clusters, sc)
            target <- if (length(neigh_clusters) == 0) {
                names(sums)[which.max(sums)]
            } else {
                neigh_clusters[which.max(sums[neigh_clusters])]
            }
            clu[members] <- target
        }
    }
    clu
}

# ICD-10 to GBD lookup
raw <- read_excel(gbd_icd_excel, sheet = 1, col_names = FALSE) %>%
    select(1:3) %>% rename(cause_raw = 1, icd10_raw = 2, icd9_raw = 3) %>%
    filter(!(toupper(icd10_raw) == "ICD10" | toupper(icd9_raw) == "ICD9"))
map_icd10 <- raw %>% filter(!is.na(icd10_raw) & icd10_raw != "") %>%
    transmute(gbd_cause_lvl3 = str_squish(as.character(cause_raw)),
              icd10_list = str_squish(as.character(icd10_raw)))
canon_icd10 <- function(x) toupper(str_replace_all(x, "[^A-Z0-9\\.,\\-]", ""))
expand_token_to_roots <- function(tok) {
    tok <- str_trim(tok); if (tok == "") return(character(0))
    parts <- str_split(tok, "-", n = 2, simplify = TRUE)
    lo <- parts[1]; hi <- ifelse(ncol(parts) == 2, parts[2], NA_character_)
    root3 <- function(code) {
        m <- str_match(code, "^([A-Z])([0-9]{1,2})")
        if (is.na(m[1,1])) return(NA_character_)
        sprintf("%s%02d", m[1,2], as.integer(m[1,3]))
    }
    if (is.na(hi) || hi == "") return(root3(lo))
    rlo <- root3(lo); rhi <- root3(hi)
    if (any(is.na(c(rlo, rhi)))) return(character(0))
    letter_lo <- substr(rlo, 1, 1); letter_hi <- substr(rhi, 1, 1)
    if (letter_lo != letter_hi) return(unique(c(rlo, rhi)))
    paste0(letter_lo, sprintf("%02d", seq(as.integer(substr(rlo, 2, 3)),
                                          as.integer(substr(rhi, 2, 3)))))
}
expand_icdlist_to_roots <- function(s) unique(unlist(map(unlist(str_split(canon_icd10(s), ",")), expand_token_to_roots)))
lookup_icd10_root <- map_icd10 %>%
    mutate(icd10_roots = map(icd10_list, expand_icdlist_to_roots)) %>%
    select(gbd_cause_lvl3, icd10_roots) %>% unnest(icd10_roots) %>%
    distinct(icd10_root = icd10_roots, gbd_cause_lvl3)
icd_root3 <- function(code) {
    code <- toupper(as.character(code))
    m <- str_match(code, "^([A-Z])([0-9]{1,2})")
    ifelse(is.na(m[,1]), NA_character_, sprintf("%s%02d", m[,2], as.integer(m[,3])))
}
if (file.exists(gbd_hierarchy_xlsx)) {
    gbd_hier <- read_excel(gbd_hierarchy_xlsx)
    lvl3_tbl <- gbd_hier %>% filter(level == 3) %>%
        select(gbd_cause_lvl3 = cause_name, parent_id)
    lvl2_tbl <- gbd_hier %>% filter(level == 2) %>%
        select(parent_id = cause_id, gbd_cause_lvl2 = cause_name)
    lookup_icd10_root <- lookup_icd10_root %>%
        left_join(lvl3_tbl, by = "gbd_cause_lvl3") %>%
        left_join(lvl2_tbl, by = "parent_id")
}

# Diversity metrics
simpson_S <- function(counts) { p <- counts / sum(counts); 1 - sum(p^2) }
rao_Q <- function(counts, dmat) { p <- counts / sum(counts); sum(outer(p, p) * dmat[names(counts), names(counts)]) }
build_dmat <- function(cause_lvl3, cause_lvl2) {
    uniq <- distinct(tibble(lvl3 = cause_lvl3, lvl2 = cause_lvl2))
    nm <- uniq$lvl3; l2 <- setNames(uniq$lvl2, uniq$lvl3)
    outer(seq_along(nm), seq_along(nm), Vectorize(function(i,j) {
        if (i == j) 0 else if (l2[nm[i]] == l2[nm[j]]) 0.5 else 1
    })) %>% `dimnames<-`(list(nm, nm))
}

# Global dissimilarity matrix
ds <- open_dataset(parquet_dir)
icd_roots_union <- ds %>% filter(!is.na(ucod)) %>% select(ucod) %>% collect() %>%
    mutate(icd10_root = icd_root3(ucod)) %>% distinct(icd10_root) %>% filter(!is.na(icd10_root))
union_map <- icd_roots_union %>%
    left_join(lookup_icd10_root, by = "icd10_root") %>%
    filter(!is.na(gbd_cause_lvl3), !is.na(gbd_cause_lvl2)) %>%
    distinct(gbd_cause_lvl3, gbd_cause_lvl2)
dmat_global <- build_dmat(union_map$gbd_cause_lvl3, union_map$gbd_cause_lvl2)
saveRDS(dmat_global, file.path(out_dir, "gbd_dissimilarity_matrix_global.rds"))

# Main loop
cluster_results <- list(); diag_rows <- list(); cluster_membership <- list()
for (pname in names(periods)) {
    cert <- ds %>% filter(year %in% periods[[pname]], !is.na(ucod)) %>%
        select(!!county_var, ucod) %>% collect() %>%
        mutate(fips = str_pad(as.character(.data[[county_var]]), 5, pad = "0")) %>%
        filter(fips %in% counties_sf$GEOID)
    death_tbl <- count(cert, fips, name = "deaths")
    clu <- merge_small_clusters(make_greedy_clusters(death_tbl, g, min_deaths),
                                death_tbl, g, min_deaths)
    cert <- mutate(cert, cluster = unname(clu[fips]))
    ccm_period <- distinct(cert, fips, cluster) %>% mutate(period = pname)
    write_csv(ccm_period, file.path(out_dir, paste0("county_cluster_membership_", pname, ".csv")))
    cluster_membership[[pname]] <- ccm_period
    cert <- cert %>% mutate(icd10_root = icd_root3(ucod)) %>%
        left_join(lookup_icd10_root, by = "icd10_root")
    counts_tbl <- cert %>% filter(!is.na(gbd_cause_lvl3), !is.na(cluster)) %>%
        count(cluster, gbd_cause_lvl3, name = "n")
    res <- counts_tbl %>% group_by(cluster) %>%
        summarise(deaths = sum(n),
                  S = simpson_S(setNames(n, gbd_cause_lvl3)),
                  I = rao_Q(setNames(n, gbd_cause_lvl3), dmat_global),
                  .groups = "drop") %>% mutate(period = pname)
    cluster_results[[pname]] <- res
}
cluster_div <- bind_rows(cluster_results)
write_csv(cluster_div, file.path(out_dir, "cluster_cod_diversity.csv"))

# Bias check
bias_df <- map_dfr(names(periods), function(pname) {
    div_period <- cluster_results[[pname]]
    cor_res <- cor.test(div_period$S, div_period$deaths, method = "spearman")
    fit <- lm(S ~ deaths, data = div_period)
    tibble(period = pname,
           spearman_rho = unname(cor_res$estimate),
           spearman_p = cor_res$p.value,
           lm_slope = coef(fit)[["deaths"]],
           lm_p = summary(fit)$coefficients[2,4],
           r2 = summary(fit)$r.squared)
})
write_csv(bias_df, file.path(out_dir, "cluster_size_bias_check.csv"))
