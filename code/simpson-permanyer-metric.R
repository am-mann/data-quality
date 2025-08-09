# ─────────────────────────────────────────────────────────────────────────────
# Date: August 8, 2025
# Outputs Simpson's S & CoD Inequality from Permanyer paper
# Similar to phillips metric but considers the degree of difference between different cause codes
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(arrow);  library(dplyr);  library(stringr); library(purrr)
    library(readr);  library(readxl); library(tidyr);   library(tibble)
    library(sf);     library(tigris); library(igraph);  library(here)
    library(fs);     library(vegan)
})

options(tigris_use_cache = TRUE, tigris_class = "sf")

# paths
out_dir <- here("output"); dir_create(out_dir)

parquet_dir <- if (dir_exists(here("data_private","mcod"))) here("data_private","mcod") else here("data_private","mcod_sample")
gbd_icd_excel      <- here("data_raw","IHME_GBD_2021_COD_CAUSE_ICD_CODE_MAP.XLSX")
gbd_hierarchy_xlsx <- here("data_raw","GBD_2021_CAUSE_HIERARCHY.XLSX")
gbd_garbage_csv    <- here("data_raw","cause-codes","gbd_garbage_codes_without_overdose.csv")

stopifnot(file.exists(gbd_icd_excel), file.exists(gbd_hierarchy_xlsx), file.exists(gbd_garbage_csv))

county_var <- "county_ihme"
min_deaths_cluster <- 200
periods <- list("1999_2004" = 1999:2004, "2005_2010" = 2005:2010, "2011_2017" = 2011:2017, "2018_2022" = 2018:2022)

use_nongarbage_for_clustering <- FALSE   # you can toggle this
membership_path <- here("output", "county_cluster_membership_all_periods.csv")

# helpers
canon_icd10 <- function(x) toupper(str_replace_all(x, "[^A-Z0-9\\.,\\-]", ""))
expand_token_to_roots <- function(tok) {
    tok <- str_trim(tok); if (tok == "") return(character(0))
    parts <- str_split(tok, "-", n = 2, simplify = TRUE)
    lo <- parts[1]; hi <- ifelse(ncol(parts) == 2, parts[2], NA_character_)
    root3 <- function(code) {
        m <- str_match(code, "^([A-Z])([0-9]{1,2})"); if (is.na(m[1,1])) return(NA_character_)
        sprintf("%s%02d", m[1,2], as.integer(m[1,3]))
    }
    if (is.na(hi) || hi == "") return(root3(lo))
    rlo <- root3(lo); rhi <- root3(hi); if (any(is.na(c(rlo, rhi)))) return(character(0))
    Llo <- substr(rlo,1,1); Lhi <- substr(rhi,1,1); if (Llo != Lhi) return(unique(c(rlo,rhi)))
    paste0(Llo, sprintf("%02d", seq(as.integer(substr(rlo,2,3)), as.integer(substr(rhi,2,3)))))
}
expand_icdlist_to_roots <- function(s) unique(unlist(map(unlist(str_split(canon_icd10(s), ",")), expand_token_to_roots)))
icd_root3 <- function(code) { code <- toupper(as.character(code))
m <- str_match(code, "^([A-Z])([0-9]{1,2})"); ifelse(is.na(m[,1]), NA_character_, sprintf("%s%02d", m[,2], as.integer(m[,3])))
}

read_garbage_roots <- function(csv_path) {
    gb <- read_csv(csv_path, show_col_types = FALSE)
    candidates <- c("code","icd10","icd10_code","ICD10","ICD_10","ICD10Code","ICD10_Code","icd_code","ICD Code","ICD-10")
    present <- intersect(candidates, names(gb))
    if (!length(present)) stop("No recognized ICD-10 code columns found in: ", csv_path, call. = FALSE)
    coalesced <- do.call(dplyr::coalesce, lapply(present, function(nm) as.character(gb[[nm]])))
    gb <- tibble(code_str = stringr::str_squish(coalesced)) %>% filter(!is.na(code_str), code_str != "")
    tokens <- gb %>% transmute(tok = str_split(canon_icd10(code_str), ",")) %>% unnest_longer(tok, values_to = "tok") %>%
        mutate(tok = str_trim(tok)) %>% filter(tok != "")
    roots <- tokens %>% mutate(root = map(tok, expand_token_to_roots)) %>% unnest(root) %>%
        filter(root != "") %>% distinct(root) %>% pull(root)
    roots
}

# loads maps and hierarchy from permanyer paper
raw_map <- read_excel(gbd_icd_excel, sheet = 1, col_names = FALSE) %>%
    select(1:3) %>% rename(cause_raw = 1, icd10_raw = 2, icd9_raw = 3) %>%
    filter(!(toupper(icd10_raw) == "ICD10" | toupper(icd9_raw) == "ICD9"))

lookup_icd10_root <- raw_map %>%
    filter(!is.na(icd10_raw) & icd10_raw != "") %>%
    transmute(gbd_cause_lvl3 = str_squish(as.character(cause_raw)),
              icd10_list     = str_squish(as.character(icd10_raw))) %>%
    mutate(icd10_roots = map(icd10_list, expand_icdlist_to_roots)) %>%
    select(gbd_cause_lvl3, icd10_roots) %>% unnest(icd10_roots) %>%
    distinct(icd10_root = icd10_roots, gbd_cause_lvl3)

stopifnot(nrow(lookup_icd10_root) > 0)

gbd_hier <- readxl::read_excel(gbd_hierarchy_xlsx, guess_max = 10000) %>%
    transmute(
        level      = suppressWarnings(as.integer(readr::parse_number(as.character(`Level`)))),
        cause_id   = `Cause ID`,
        cause_name = as.character(`Cause Name`),
        parent_id  = `Parent ID`
    )

lvl3 <- gbd_hier %>% filter(level == 3) %>% select(cause_id3 = cause_id, gbd_cause_lvl3 = cause_name, parent_id3 = parent_id)
lvl2 <- gbd_hier %>% filter(level == 2) %>% select(cause_id2 = cause_id, gbd_cause_lvl2 = cause_name, parent_id2 = parent_id)
lvl1 <- gbd_hier %>% filter(level == 1) %>% select(cause_id1 = cause_id, gbd_cause_lvl1 = cause_name)

l3_l2 <- lvl3 %>% left_join(lvl2, by = c("parent_id3" = "cause_id2"))
l3_l2_l1 <- l3_l2 %>% left_join(lvl1, by = c("parent_id2" = "cause_id1")) %>%
    select(gbd_cause_lvl3, gbd_cause_lvl2, gbd_cause_lvl1)

lookup_icd10_root <- lookup_icd10_root %>%
    left_join(l3_l2_l1, by = "gbd_cause_lvl3")

stopifnot(all(c("gbd_cause_lvl2","gbd_cause_lvl1") %in% names(lookup_icd10_root)))

# d = 0 if same L3; 1/3 if different L3 but same L2; 2/3 if diff L2 same L1; 1 if diff L1
build_tree_dmat <- function(map_l3_l2_l1) {
    causes <- unique(map_l3_l2_l1$gbd_cause_lvl3)
    L2 <- setNames(map_l3_l2_l1$gbd_cause_lvl2, map_l3_l2_l1$gbd_cause_lvl3)
    L1 <- setNames(map_l3_l2_l1$gbd_cause_lvl1, map_l3_l2_l1$gbd_cause_lvl3)
    n <- length(causes)
    d <- matrix(0, n, n, dimnames = list(causes, causes))
    for (i in seq_len(n)) for (j in seq_len(n)) {
        if (i == j) { d[i,j] <- 0; next }
        ci <- causes[i]; cj <- causes[j]
        if (L2[ci] == L2[cj]) d[i,j] <- 1/3
        else if (L1[ci] == L1[cj]) d[i,j] <- 2/3
        else d[i,j] <- 1
    }
    d
}


garbage_roots <- read_garbage_roots(gbd_garbage_csv)

ds <- open_dataset(parquet_dir)

make_greedy_clusters <- function(death_tbl, g, min_deaths) {
    deaths    <- setNames(death_tbl$deaths, death_tbl$fips)
    to_assign <- names(deaths)
    clusters  <- setNames(rep(NA_character_, length(deaths)), to_assign)
    cid <- 1L
    while (length(to_assign) > 0) {
        this <- to_assign[which.max(deaths[to_assign])]
        if (deaths[this] >= min_deaths) {
            clusters[this] <- paste0("C", cid); to_assign <- setdiff(to_assign, this); cid <- cid + 1L; next
        }
        cluster <- this; total <- deaths[this]; avail <- setdiff(to_assign, this)
        repeat {
            nbrs <- intersect(setdiff(unique(unlist(lapply(cluster, function(v) neighbors(g, v) %>% names()))), cluster), avail)
            if (!length(nbrs)) break
            best <- nbrs[which.min(abs((total + deaths[nbrs]) - min_deaths))]
            new_tot <- total + deaths[best]
            if (new_tot <= min_deaths * 1.5 || total < min_deaths) {
                cluster <- c(cluster, best); total <- new_tot; avail <- setdiff(avail, best)
            } else break
        }
        clusters[cluster] <- paste0("C", cid)
        to_assign <- setdiff(to_assign, cluster); cid <- cid + 1L
    }
    clusters
}
merge_small_clusters <- function(clu, death_tbl, g, min_deaths = 200) {
    deaths_vec <- setNames(death_tbl$deaths, death_tbl$fips); deaths_vec <- deaths_vec[names(clu)]
    repeat {
        sums <- tapply(deaths_vec, clu, sum, na.rm = TRUE); sums <- sums[!is.na(names(sums))]
        small <- names(sums)[sums < min_deaths]; if (!length(small)) break
        for (sc in small) {
            members <- names(clu)[clu == sc]
            neigh_clusters <- unique(clu[neighbors(g, members) %>% names()]); neigh_clusters <- setdiff(na.omit(neigh_clusters), sc)
            candidates <- intersect(neigh_clusters, names(sums))
            target <- if (length(candidates)) { s <- sums[candidates]; s[is.na(s)] <- -Inf; candidates[which.max(s)] }
            else { bigs <- setdiff(names(sums), sc); if (length(bigs)) bigs[which.max(sums[bigs])] else NA_character_ }
            if (is.na(target) || !length(target)) next
            clu[members] <- target
        }
    }
    clu
}

get_or_build_membership <- function() {
    if (file.exists(membership_path)) {
        read_csv(membership_path, show_col_types = FALSE) %>%
            mutate(fips = str_pad(as.character(fips), 5, pad = "0")) %>%
            select(fips, period, cluster)
    } else {
        counties_sf <- counties(cb = TRUE, year = 2020) %>% select(GEOID, geometry)
        adj <- st_touches(counties_sf)
        edge_df <- tibble(from = rep(counties_sf$GEOID, lengths(adj)),
                          to   = counties_sf$GEOID[unlist(adj)]) %>% filter(from < to)
        g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = counties_sf$GEOID)
        
        out <- list()
        for (pname in names(periods)) {
            cert <- ds %>% filter(year %in% periods[[pname]], !is.na(ucod)) %>%
                select(!!county_var, ucod) %>% collect() %>%
                mutate(fips = str_pad(as.character(.data[[county_var]]), 5, pad = "0")) %>%
                filter(fips %in% counties_sf$GEOID)
            death_tbl <- if (isTRUE(use_nongarbage_for_clustering)) {
                cert %>% mutate(icd10_root = icd_root3(ucod)) %>% filter(!(icd10_root %in% garbage_roots)) %>% count(fips, name = "deaths")
            } else count(cert, fips, name = "deaths")
            clu <- merge_small_clusters(make_greedy_clusters(death_tbl, g, min_deaths_cluster),
                                        death_tbl, g, min_deaths_cluster)
            out[[pname]] <- tibble(fips = names(clu), cluster = unname(clu)) %>% mutate(period = pname)
        }
        bind_rows(out) %>% write_csv(membership_path)
    }
}

ccm <- get_or_build_membership()

icd_roots_union <- ds %>%
    filter(!is.na(ucod)) %>% select(ucod) %>% collect() %>%
    mutate(icd10_root = icd_root3(ucod)) %>%
    distinct(icd10_root) %>%
    filter(!is.na(icd10_root), !(icd10_root %in% garbage_roots))

union_map <- icd_roots_union %>%
    left_join(lookup_icd10_root, by = "icd10_root") %>%
    filter(!is.na(gbd_cause_lvl3), !is.na(gbd_cause_lvl2), !is.na(gbd_cause_lvl1)) %>%
    distinct(gbd_cause_lvl3, gbd_cause_lvl2, gbd_cause_lvl1)

dmat_global <- build_tree_dmat(union_map)

# S and I helpers
simpson_S_from_counts <- function(counts) {
    if (length(counts) == 0 || sum(counts) == 0) return(NA_real_)
    p <- counts / sum(counts); 1 - sum(p^2)
}
I_from_counts <- function(counts, dmat) {
    if (length(counts) == 0 || sum(counts) == 0) return(NA_real_)
    nm <- names(counts); nm <- nm[nm %in% rownames(dmat)]
    counts <- counts[nm]; if (!length(counts)) return(NA_real_)
    p <- counts / sum(counts); sum(outer(p, p) * dmat[names(p), names(p)])
}
I_contribs <- function(counts, dmat) {
    nm <- names(counts); nm <- nm[nm %in% rownames(dmat)]
    counts <- counts[nm]; if (!length(counts) || sum(counts)==0) return(tibble(lvl3 = character(0), C = numeric(0)))
    p <- counts / sum(counts)
    # C_c = p_c * sum_i d_ci * p_i
    C <- sapply(names(p), function(c) p[c] * sum(dmat[c, names(p)] * p))
    tibble(lvl3 = names(C), C = as.numeric(C))
}

rarefy_counts <- function(counts_vec, N_ref) {
    if (sum(counts_vec) < N_ref || N_ref <= 0) return(NA)
    mat <- matrix(counts_vec, nrow = 1); colnames(mat) <- names(counts_vec)
    as.numeric(vegan::rrarefy(mat, N_ref)[1,])
}

results <- list(); contrib_rows <- list(); meta_lines <- c()

for (pname in names(periods)) {
    yrs <- periods[[pname]]

    cert <- ds %>% filter(year %in% yrs, !is.na(ucod)) %>%
        select(!!county_var, ucod) %>% collect() %>%
        mutate(fips = str_pad(as.character(.data[[county_var]]), 5, pad = "0"),
               icd10_root = icd_root3(ucod)) %>%
        filter(!(icd10_root %in% garbage_roots)) %>%
        left_join(lookup_icd10_root, by = "icd10_root") %>%
        filter(!is.na(gbd_cause_lvl3))
    
    cert <- cert %>% left_join(ccm %>% filter(period == pname), by = "fips") %>% filter(!is.na(cluster))
    
    counts_tbl <- cert %>% count(cluster, gbd_cause_lvl3, name = "n")
    deaths_per_cluster <- counts_tbl %>% group_by(cluster) %>% summarise(deaths = sum(n), .groups = "drop")
    
    # choose rarefaction N_ref for this period (use 25th percentile of deaths among clusters with >=100 deaths)
    elig <- deaths_per_cluster$deaths[deaths_per_cluster$deaths >= 100]
    N_ref <- if (length(elig)) max(50, floor(quantile(elig, 0.25))) else NA_integer_
    meta_lines <- c(meta_lines, sprintf("%s: N_ref=%s (clusters=%d)", pname, as.character(N_ref), nrow(deaths_per_cluster)))
    
    # compute per-cluster metrics
    by_cluster <- counts_tbl %>% group_by(cluster) %>%
        summarise(
            deaths = sum(n),
            S_raw = { v <- setNames(n, gbd_cause_lvl3); simpson_S_from_counts(v) },
            I_raw = { v <- setNames(n, gbd_cause_lvl3); I_from_counts(v, dmat_global) },
            .groups = "drop"
        )
    
    add_rare <- counts_tbl %>% group_by(cluster) %>% group_map(function(df, key){
        cl <- key$cluster[[1]]
        v <- setNames(df$n, df$gbd_cause_lvl3)
        if (is.na(N_ref) || sum(v) < N_ref) {
            tibble(cluster = cl, S_rare = NA_real_, I_rare = NA_real_)
        } else {
            vr <- rarefy_counts(v, N_ref); names(vr) <- names(v)
            S_rare <- simpson_S_from_counts(vr)
            I_rare <- I_from_counts(vr, dmat_global)
            tibble(cluster = cl, S_rare = S_rare, I_rare = I_rare)
        }
    }) %>% bind_rows()
    
    # residual (size-adjusted) using log(deaths) — on raw metrics
    tmp <- by_cluster %>% left_join(add_rare, by = "cluster")
    if (nrow(tmp) >= 4) {
        fitS <- try(lm(S_raw ~ log1p(deaths), data = tmp), silent = TRUE)
        fitI <- try(lm(I_raw ~ log1p(deaths), data = tmp), silent = TRUE)
        tmp$S_resid <- if (inherits(fitS, "try-error")) NA_real_ else resid(fitS)
        tmp$I_resid <- if (inherits(fitI, "try-error")) NA_real_ else resid(fitI)
    } else {
        tmp$S_resid <- NA_real_; tmp$I_resid <- NA_real_
    }
    tmp$period <- pname
    results[[pname]] <- tmp
    
    # cause-specific contributions using rarefied probabilities when available (fallback to raw)
    con <- counts_tbl %>% group_by(cluster) %>% group_map(function(df, key){
        cl <- key$cluster[[1]]
        v <- setNames(df$n, df$gbd_cause_lvl3)
        use_v <- v
        if (!is.na(N_ref) && sum(v) >= N_ref) {
            vr <- rarefy_counts(v, N_ref); names(vr) <- names(v); use_v <- vr
        }
        I_contribs(use_v, dmat_global) %>% transmute(cluster = cl, lvl3 = lvl3, contribution_C = C)
    }) %>% bind_rows() %>% mutate(period = pname)
    contrib_rows[[pname]] <- con
}

quality_df <- bind_rows(results) %>% arrange(period, cluster)
write_csv(quality_df, file.path(out_dir, "cod_quality_indices.csv"))

contrib_df <- bind_rows(contrib_rows) %>% arrange(period, cluster, lvl3)
write_csv(contrib_df, file.path(out_dir, "cod_quality_contributions.csv"))

write_lines(meta_lines, file.path(out_dir, "cod_quality_meta.txt"))

message("Code finished running.")
