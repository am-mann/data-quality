# =============================================================================
# Mapping & Validation: One-File R Script
# Converts your multi-chunk Rmd into a single executable R script.
# Assumes your project has the same folder structure and input files referenced
# in the original notebook (data/, data_raw/, output/, figures/…).
# Run with: Rscript mapping_pipeline.R
# =============================================================================

# ──────────────────────────────────────────────────────────────
# 0) Packages & global options
# ──────────────────────────────────────────────────────────────
required <- c(
  "sf","tigris","readr","dplyr","stringr","purrr","ggplot2","here","scales","patchwork",
  "igraph","tidyr","tibble","stats"
)
missing <- setdiff(required, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
invisible(lapply(required, library, character.only = TRUE))

options(tigris_use_cache = TRUE, tigris_class = "sf")
sf::sf_use_s2(FALSE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ──────────────────────────────────────────────────────────────
# Common helpers used across sections
# ──────────────────────────────────────────────────────────────
crs_proj <- 2163

st_shift <- function(x, offset = c(0,0)) { sf::st_geometry(x) <- sf::st_geometry(x) + offset; x }
st_scale <- function(x, factor = 1) { ctr <- sf::st_centroid(sf::st_union(x)); sf::st_geometry(x) <- (sf::st_geometry(x) - ctr) * factor + ctr; x }

normalize_geoid <- function(sf, year = NA) {
  nms <- names(sf)
  geoid_col <- grep("^GEOID", nms, value = TRUE)[1]
  if (!is.na(geoid_col)) {
    sf <- dplyr::rename(sf, GEOID = !!rlang::sym(geoid_col))
    sf$GEOID <- stringr::str_pad(as.character(sf$GEOID), 5, pad = "0")
    return(sf)
  }
  state_col  <- grep("^STATE.*FP$", nms, value = TRUE)[1]
  county_col <- grep("^COUNTY.*FP$", nms, value = TRUE)[1]
  if (!is.na(state_col) && !is.na(county_col)) {
    sf$GEOID <- paste0(
      stringr::str_pad(as.character(sf[[state_col]]),  2, pad = "0"),
      stringr::str_pad(as.character(sf[[county_col]]), 3, pad = "0")
    )
    return(sf)
  }
  combo_col <- grep("^CNTYIDFP$", nms, value = TRUE)[1]
  if (!is.na(combo_col)) {
    sf <- dplyr::rename(sf, GEOID = !!rlang::sym(combo_col))
    sf$GEOID <- stringr::str_pad(as.character(sf$GEOID), 5, pad = "0")
    return(sf)
  }
  stop("No GEOID-compatible columns found in shapefile", if (!is.na(year)) paste0(" for year ", year) else "")
}

mode_val <- function(x) {
  ux <- na.omit(x)
  if (!length(ux)) return(NA_character_)
  names(sort(table(ux), decreasing = TRUE))[1]
}

safe_div <- function(num, den) ifelse(is.finite(num) & is.finite(den) & den > 0, num / den, NA_real_)
compute_limits_symmetric <- function(x) {
  x <- x[is.finite(x)]; if (!length(x)) return(c(-1,1))
  L <- as.numeric(stats::quantile(abs(x), 0.98, na.rm = TRUE, names = FALSE))
  if (!is.finite(L) || L == 0) L <- 1
  c(-L, L)
}
sanitize <- function(x) {
  x %>% stringr::str_replace_all("[^A-Za-z0-9]+","_") %>%
    stringr::str_replace_all("^_+|_+$","") %>% tolower()
}
pick_col <- function(nms, candidates) { hit <- intersect(candidates, nms); if (length(hit)) hit[1] else NA_character_ }
mode_str <- function(x){ ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
std_fips <- function(x) stringr::str_pad(as.character(x), 5, pad = "0")

period_of_year <- function(y){
  dplyr::case_when(
    y %in% 1999:2005 ~ "1999_2005",
    y %in% 2006:2012 ~ "2006_2012",
    y %in% 2013:2019 ~ "2013_2019",
    y %in% 2020:2022 ~ "2020_2022",
    TRUE ~ NA_character_
  )
}
standardize_period <- function(p){
  p <- as.character(p)
  out <- ifelse(p %in% c("1999_2005","2006_2012","2013_2019","2020_2022"), p, NA_character_)
  bad <- is.na(out) & grepl("^\\d{4}_\\d{4}$", p)
  if (any(bad)) {
    start <- suppressWarnings(as.integer(sub("_(\\d{4})$", "", p[bad])))
    out[bad] <- period_of_year(start)
  }
  out
}

# ──────────────────────────────────────────────────────────────
# INPUTS & paths
# ──────────────────────────────────────────────────────────────
metrics_path_gz <- here::here("data","county_year_quality_metrics.csv.gz")
metrics_path    <- if (file.exists(metrics_path_gz)) metrics_path_gz else here::here("data","county_year_quality_metrics.csv")

# Load county-year main dataset
dq <- readr::read_csv(metrics_path, show_col_types = FALSE) %>%
  dplyr::mutate(
    year        = suppressWarnings(as.integer(year)),
    county_ihme = stringr::str_pad(as.character(county_ihme), 5, pad = "0"),
    time_window = dplyr::case_when(
      year >= 1999 & year <= 2005 ~ "1999_2005",
      year >= 2006 & year <= 2012 ~ "2006_2012",
      year >= 2013 & year <= 2019 ~ "2013_2019",
      year >= 2020 & year <= 2022 ~ "2020_2022",
      TRUE ~ NA_character_
    )
  )

# IHME crosswalk
load(here::here("data_raw","ihme_fips.rda"))  # provides ihme_fips
ihme_map <- ihme_fips %>%
  dplyr::transmute(
    GEOID       = stringr::str_pad(as.character(orig_fips), 5, pad = "0"),
    county_ihme = stringr::str_pad(as.character(ihme_fips), 5, pad = "0")
  ) %>%
  dplyr::distinct()

# Output directories
dir.create(here::here("figures","5yr_avg_stable"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("figures","validation_zscore"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("figures","phillips_detail"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("output","cluster_scores"), recursive = TRUE, showWarnings = FALSE)

# ──────────────────────────────────────────────────────────────
# SECTION A — Temporality-stable group 5-year average maps
# ──────────────────────────────────────────────────────────────
vars_to_map <- c("prop_garbage", "DQ_rec_ig_frac_mean_garbage", "DQ_overall", "DQ_rec_ig_abs_mean", "foreman_garbage", "RI")
vars_present <- intersect(vars_to_map, names(dq))
if (!length(vars_present)) stop("None of vars_to_map exist in dq. Check column names.")
vars_missing <- setdiff(vars_to_map, vars_present)
if (length(vars_missing)) message("Skipping missing vars: ", paste(vars_missing, collapse = ", "))

dq_avg <- dq %>%
  dplyr::group_by(county_ihme, time_window) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(vars_present), ~ mean(.x, na.rm = TRUE)),
    n_years = dplyr::n(),
    .groups = "drop"
  )

build_groups_sf <- function(year) {
  counties_raw <- tigris::counties(year = year, cb = TRUE, class = "sf") %>%
    sf::st_zm(drop = TRUE, what = "ZM")
  counties_norm <- normalize_geoid(counties_raw, year) %>%
    sf::st_transform(crs_proj)

  groups_sf <- counties_norm %>%
    dplyr::left_join(ihme_map, by = "GEOID") %>%
    dplyr::mutate(county_ihme = dplyr::coalesce(county_ihme, GEOID)) %>%
    dplyr::group_by(county_ihme) %>%
    dplyr::summarise(.groups = "drop")

  alaska   <- groups_sf %>% dplyr::filter(substr(county_ihme,1,2)=="02") %>%
    { ctr <- sf::st_centroid(sf::st_union(.)); sf::st_geometry(.) <- (sf::st_geometry(.) - ctr)*0.40 + ctr; . } %>%
    { sf::st_geometry(.) <- sf::st_geometry(.) + c(1300000, -4900000); . } %>% sf::st_set_crs(crs_proj)
  hawaii   <- groups_sf %>% dplyr::filter(substr(county_ihme,1,2)=="15") %>%
    { ctr <- sf::st_centroid(sf::st_union(.)); sf::st_geometry(.) <- (sf::st_geometry(.) - ctr)*1.50 + ctr; . } %>%
    { sf::st_geometry(.) <- sf::st_geometry(.) + c(5200000, -1400000); . } %>% sf::st_set_crs(crs_proj)
  mainland <- groups_sf %>% dplyr::filter(!substr(county_ihme,1,2) %in% c("02","15"))
  dplyr::bind_rows(mainland, alaska, hawaii)
}

shapefile_years_A <- c("1999–2005"=2000, "2006–2012"=2010, "2013–2019"=2019, "2020–2022"=2020)

joined_by_period_A <- purrr::imap(
  shapefile_years_A,
  function(year, window) {
    shape_all <- build_groups_sf(year)
    dplyr::left_join(shape_all, dplyr::filter(dq_avg, time_window == gsub("–","_", window)), by = "county_ihme")
  }
)

make_map_A <- function(sf_data, var, title = NULL, palette = "Reds") {
  states <- tigris::states(cb = TRUE, class = "sf") |> st_transform(2163)
  fill_scale <- if (var == "prop_garbage") {
    scale_fill_distiller(palette = palette, limits = c(0.15, 0.4), oob = scales::squish, direction = -1, na.value = "grey90")
  } else if (var == "RI") {
    scale_fill_distiller(palette = palette, limits = c(0.15, 0.3), oob = scales::squish, direction = 1, na.value = "grey90")
  } else if (var == "DQ_overall") {
    scale_fill_distiller(palette = "RdBlu", limits = c(0,0.9), oob = scales::squish, direction = 1, na.value = "grey90")
  } else if (startsWith(var, "pct_") || var == "overall_completeness_pct") {
    scale_fill_distiller(palette = palette, limits = c(0, 100), oob = scales::squish, direction = 1, na.value = "grey90")
  } else {
    scale_fill_distiller(palette = palette, oob = scales::squish, na.value = "grey90", direction = -1)
  }

  ggplot() +
    geom_sf(data = sf_data, aes(fill = .data[[var]]), colour = NA) +
    geom_sf(data = states, fill = NA, colour = "white", linewidth = 0.3) +
    fill_scale +
    coord_sf(xlim = c(-2500000, 2500000),
             ylim = c(-2200000, 730000),
             expand = FALSE) +
    labs(title = title %||% "", fill = NULL) +
    theme_void() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 9)
    )
}

out_dir_A <- here::here("figures","5yr_avg_stable")
purrr::walk(vars_present, function(var) {
  p_1999 <- make_map_A(joined_by_period_A[["1999–2005"]], var, "1999–2005")
  p_2006 <- make_map_A(joined_by_period_A[["2006–2012"]], var, "2006–2012")
  p_2013 <- make_map_A(joined_by_period_A[["2013–2019"]], var, "2013–2019")
  p_2020 <- make_map_A(joined_by_period_A[["2020–2022"]], var, "2020–2022")
  combined <- (p_1999 | p_2006) / (p_2013 | p_2020)
  ggsave(filename = file.path(out_dir_A, paste0(var, "_4panel.png")),
         plot = combined, width = 12, height = 9, dpi = 320)
})
message("✓ [A] 4-panel maps saved to: ", out_dir_A)

# ──────────────────────────────────────────────────────────────
# SECTION B — RI maps with light clustering (n_cert proxy)
# ──────────────────────────────────────────────────────────────
periods_vec <- c("1999_2005","2006_2012","2013_2019","2020_2022")
shapefile_years_B <- c("1999_2005"=2000, "2006_2012"=2010, "2013_2019"=2019, "2020_2022"=2020)

req_cols_B <- c("RI","RI_post_only","n_cert")
stopifnot(all(req_cols_B %in% names(dq)))

has_jsd <- "RI_jsd" %in% names(dq)
if (!has_jsd) message("ℹ 'RI_jsd' not found in input; JSD maps will be skipped.")

dq$gb_n <- suppressWarnings(as.numeric(dq$n_cert))
dq$gb_n[is.na(dq$gb_n)] <- 0

gb_by_window <- dq %>%
  dplyr::filter(!is.na(time_window)) %>%
  dplyr::group_by(county_ihme, time_window) %>%
  dplyr::summarise(gb = sum(gb_n, na.rm = TRUE), .groups = "drop")

gb_wide <- gb_by_window %>%
  tidyr::pivot_wider(names_from = time_window, values_from = gb, values_fill = 0)
for (w in periods_vec) if (!w %in% names(gb_wide)) gb_wide[[w]] <- 0

robust_tbl_B <- gb_wide %>%
  dplyr::mutate(robust_gb = pmin(`1999_2005`,`2006_2012`,`2013_2019`,`2020_2022`, na.rm = TRUE)) %>%
  dplyr::transmute(county_ihme, robust_gb = as.numeric(robust_gb))

build_ihme_groups_sf <- function(year) {
  counties_raw <- tigris::counties(year = year, cb = TRUE, class = "sf") %>%
    sf::st_zm(drop = TRUE, what = "ZM") %>% sf::st_make_valid()
  counties_norm <- normalize_geoid(counties_raw, year) %>%
    sf::st_transform(crs_proj) %>% sf::st_make_valid()
  counties_norm %>%
    dplyr::left_join(ihme_map, by = "GEOID") %>%
    dplyr::mutate(county_ihme = dplyr::coalesce(county_ihme, GEOID)) %>%
    dplyr::group_by(county_ihme) %>%
    dplyr::summarise(.groups = "drop")
}

# Light clustering (min_gb = 2000)
albers_5070 <- 5070
ihme2020 <- build_ihme_groups_sf(2020) %>% sf::st_transform(albers_5070)
adj_list <- sf::st_touches(ihme2020)
neighbors_of <- setNames(lapply(seq_len(nrow(ihme2020)), function(i) ihme2020$county_ihme[adj_list[[i]]]),
                         ihme2020$county_ihme)
min_gb <- 2000L
nodes_tbl <- ihme2020 %>%
  sf::st_drop_geometry() %>%
  dplyr::select(county_ihme) %>%
  dplyr::left_join(robust_tbl_B, by = "county_ihme") %>%
  dplyr::mutate(robust_gb = dplyr::coalesce(robust_gb, 0))
gb_vec    <- setNames(nodes_tbl$robust_gb, nodes_tbl$county_ihme)
to_assign <- names(gb_vec)
clusters  <- setNames(rep(NA_character_, length(gb_vec)), to_assign)
cid <- 1L
while (length(to_assign) > 0) {
  seed  <- to_assign[which.max(gb_vec[to_assign])]
  cur   <- seed
  total <- gb_vec[seed]
  avail <- setdiff(to_assign, seed)
  repeat {
    nbrs <- unique(unlist(neighbors_of[cur], use.names = FALSE))
    nbrs <- setdiff(intersect(nbrs, avail), cur)
    if (!length(nbrs)) break
    cand_totals <- total + gb_vec[nbrs]
    best_idx <- if (any(cand_totals < min_gb)) which.max(cand_totals) else which.min(abs(cand_totals - min_gb))
    best <- nbrs[best_idx]
    new_total <- total + gb_vec[best]
    if (is.na(new_total)) break
    if (new_total <= min_gb * 1.6 || total < min_gb) {
      cur   <- c(cur, best)
      total <- new_total
      avail <- setdiff(avail, best)
      if (total >= min_gb) break
    } else break
  }
  clusters[cur] <- paste0("CL", cid)
  to_assign <- setdiff(to_assign, cur)
  cid <- cid + 1L
}
cluster_map <- tibble::tibble(county_ihme = names(clusters), cluster_id = unname(clusters))

agg_metric <- function(varname) {
  dq %>%
    dplyr::filter(!is.na(time_window)) %>%
    dplyr::inner_join(cluster_map, by = "county_ihme") %>%
    dplyr::group_by(cluster_id, time_window) %>%
    dplyr::summarise(val = mean(.data[[varname]], na.rm = TRUE),
                     n_counties = dplyr::n_distinct(county_ihme),
                     .groups = "drop") %>%
    tidyr::complete(cluster_id, time_window = periods_vec,
                    fill = list(val = NA_real_, n_counties = 0L))
}
dq_cluster_RI   <- agg_metric("RI") %>% dplyr::rename(RI = val)
dq_cluster_post <- agg_metric("RI_post_only") %>% dplyr::rename(RI_post_only = val)
dq_cluster_jsd  <- if (has_jsd) agg_metric("RI_jsd") %>% dplyr::rename(RI_jsd = val) else NULL

build_clusters_sf <- function(year) {
  base <- build_ihme_groups_sf(year)
  base %>%
    dplyr::left_join(cluster_map, by = "county_ihme") %>%
    dplyr::filter(!is.na(cluster_id)) %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(.groups = "drop") %>%
    sf::st_transform(crs_proj) %>%
    sf::st_make_valid()
}
join_shapes <- function(metric_df) {
  lapply(names(shapefile_years_B), function(window) {
    year <- unname(shapefile_years_B[[window]])
    shape_all <- build_clusters_sf(year)
    dplyr::left_join(shape_all, dplyr::filter(metric_df, time_window == window), by = "cluster_id")
  }) %>% setNames(names(shapefile_years_B))
}
joined_by_period_RI   <- join_shapes(dq_cluster_RI)
joined_by_period_post <- join_shapes(dq_cluster_post)
joined_by_period_jsd  <- if (!is.null(dq_cluster_jsd)) join_shapes(dq_cluster_jsd) else NULL

states_outline <- tigris::states(cb = TRUE, class = "sf") |>
  sf::st_transform(crs_proj) |>
  sf::st_make_valid()

make_map_numeric <- function(sf_data, fill_col, limits, title) {
  ggplot() +
    geom_sf(data = sf_data, aes(fill = .data[[fill_col]]), colour = NA) +
    geom_sf(data = states_outline, fill = NA, colour = "white", linewidth = 0.3) +
    scale_fill_distiller(palette = "Reds", limits = limits,
                         oob = scales::squish, direction = -1,
                         na.value = "grey90") +
    coord_sf(xlim = c(-2500000, 2500000), ylim = c(-2200000, 730000), expand = FALSE) +
    labs(title = title, fill = NULL) +
    theme_void() +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 9))
}
make_map_RI   <- function(sf_data, title) make_map_numeric(sf_data, "RI", c(0.02, 0.05), title)
make_map_post <- function(sf_data, title) make_map_numeric(sf_data, "RI_post_only", c(0.08, 0.12), title)
make_map_jsd  <- function(sf_data, title) make_map_numeric(sf_data, "RI_jsd", c(0.02, 0.04), title)

out_dir_B <- here::here("figures", "5yr_avg_stable")
p_1999_RI <- make_map_RI(joined_by_period_RI[["1999_2005"]], "1999_2005")
p_2006_RI <- make_map_RI(joined_by_period_RI[["2006_2012"]], "2006_2012")
p_2013_RI <- make_map_RI(joined_by_period_RI[["2013_2019"]], "2013_2019")
p_2020_RI <- make_map_RI(joined_by_period_RI[["2020_2022"]], "2020_2022")
combined_RI <- (p_1999_RI | p_2006_RI) / (p_2013_RI | p_2020_RI)
ggsave(file.path(out_dir_B, "RI_clustered_garbage150_4panel.png"),
       plot = combined_RI, width = 12, height = 9, dpi = 320)

p_1999_post <- make_map_post(joined_by_period_post[["1999_2005"]], "1999_2005")
p_2006_post <- make_map_post(joined_by_period_post[["2006_2012"]], "2006_2012")
p_2013_post <- make_map_post(joined_by_period_post[["2013_2019"]], "2013_2019")
p_2020_post <- make_map_post(joined_by_period_post[["2020_2022"]], "2020_2022")
combined_post <- (p_1999_post | p_2006_post) / (p_2013_post | p_2020_post)
ggsave(file.path(out_dir_B, "RI_post_only_clustered_garbage150_4panel.png"),
       plot = combined_post, width = 12, height = 9, dpi = 320)

if (!is.null(joined_by_period_jsd)) {
  p_1999_jsd <- make_map_jsd(joined_by_period_jsd[["1999_2005"]], "1999_2005")
  p_2006_jsd <- make_map_jsd(joined_by_period_jsd[["2006_2012"]], "2006_2012")
  p_2013_jsd <- make_map_jsd(joined_by_period_jsd[["2013_2019"]], "2013_2019")
  p_2020_jsd <- make_map_jsd(joined_by_period_jsd[["2020_2022"]], "2020_2022")
  combined_jsd <- (p_1999_jsd | p_2006_jsd) / (p_2013_jsd | p_2020_jsd)
  ggsave(file.path(out_dir_B, "RI_jsd_clustered_garbage150_4panel.png"),
         plot = combined_jsd, width = 12, height = 9, dpi = 320)
} else {
  message("⚠ [B] Skipping JSD maps: 'RI_jsd' not present in input file.")
}

# ──────────────────────────────────────────────────────────────
# SECTION C — Average K-L divergence (RI) time series
# ──────────────────────────────────────────────────────────────
.wmean <- function(x, w) {
  w <- ifelse(is.finite(w), w, NA_real_)
  if (all(is.na(w)) || sum(w, na.rm = TRUE) == 0) return(NA_real_)
  stats::weighted.mean(x, w, na.rm = TRUE)
}
ri_ts <- dq %>%
  dplyr::filter(!is.na(year), year >= 1999, year <= 2022, is.finite(RI)) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    wmean_RI  = .wmean(RI, N_garb_g1g2g4g5g6g7),  # adjust if your weight col differs
    n_ctyyr   = dplyr::n(),
    .groups = "drop"
  )
readr::write_csv(ri_ts, file.path(here::here("figures","5yr_avg_stable"), "RI_time_series_1999_2022.csv"))
p_ri_ts <- ggplot(ri_ts, aes(x = year, y = wmean_RI)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.7) +
  labs(x = NULL, y = "Normalized K-L Divergence") +
  scale_x_continuous(breaks = seq(2000, 2022, by = 2)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
ggsave(file.path(here::here("figures","5yr_avg_stable"), "RI_time_series_1999_2022.png"),
       plot = p_ri_ts, width = 8.5, height = 4.8, dpi = 320)

# ──────────────────────────────────────────────────────────────
# SECTION D — % missing overdoses map with clusters (4 panels)
# ──────────────────────────────────────────────────────────────
periods_list_D <- list(
  "1999_2005" = 1999:2005,
  "2006_2012" = 2006:2012,
  "2013_2019" = 2013:2019,
  "2020_2022" = 2020:2022
)
periods_vec_D <- names(periods_list_D)
stopifnot(all(c("overd_n","pct_overd_miss") %in% names(dq)))

overd_by_window <- dq %>%
  dplyr::filter(!is.na(time_window)) %>%
  dplyr::group_by(county_ihme, time_window) %>%
  dplyr::summarise(overd = sum(overd_n, na.rm = TRUE), .groups = "drop")
overd_wide <- overd_by_window %>%
  tidyr::pivot_wider(names_from = time_window, values_from = overd, values_fill = 0)
for (w in periods_vec_D) if (!w %in% names(overd_wide)) overd_wide[[w]] <- 0
robust_tbl_D <- overd_wide %>%
  dplyr::mutate(robust_overd = pmin(`1999_2005`,`2006_2012`,`2013_2019`,`2020_2022`, na.rm = TRUE)) %>%
  dplyr::transmute(county_ihme, robust_overd = as.numeric(robust_overd))

build_ihme_groups_sf_state <- function(year) {
  counties_raw <- tigris::counties(year = year, cb = TRUE, class = "sf") %>% sf::st_zm(drop = TRUE, what = "ZM")
  counties_norm <- normalize_geoid(counties_raw, year) %>% sf::st_transform(crs_proj)
  counties_norm %>%
    dplyr::left_join(ihme_map, by = "GEOID") %>%
    dplyr::mutate(
      county_ihme = dplyr::coalesce(county_ihme, GEOID),
      STATEFP = as.character(STATEFP)
    ) %>%
    dplyr::group_by(county_ihme) %>%
    dplyr::summarise(state = mode_val(STATEFP), .groups = "drop")
}

ihme2020_D <- build_ihme_groups_sf_state(2020) %>% sf::st_transform(5070)
adj_list_D <- sf::st_touches(ihme2020_D)
neighbors_of_D <- setNames(lapply(seq_len(nrow(ihme2020_D)), function(i) ihme2020_D$county_ihme[adj_list_D[[i]]]),
                           ihme2020_D$county_ihme)
state_of <- setNames(ihme2020_D$state, ihme2020_D$county_ihme)
min_overd <- 50L
nodes_tbl_D <- ihme2020_D %>% sf::st_drop_geometry() %>%
  dplyr::select(county_ihme, state) %>%
  dplyr::left_join(robust_tbl_D, by = "county_ihme") %>%
  dplyr::mutate(robust_overd = dplyr::coalesce(robust_overd, 0))
overd_vec <- setNames(nodes_tbl_D$robust_overd, nodes_tbl_D$county_ihme)
to_assign <- names(overd_vec)
clustersD  <- setNames(rep(NA_character_, length(overd_vec)), to_assign)
cid <- 1L
while (length(to_assign) > 0) {
  seed       <- to_assign[which.max(overd_vec[to_assign])]
  home_state <- state_of[seed]
  cur        <- seed
  total      <- overd_vec[seed]
  avail      <- setdiff(to_assign, seed)
  repeat {
    nbrs_all <- unique(unlist(neighbors_of_D[cur], use.names = FALSE))
    nbrs_all <- setdiff(intersect(nbrs_all, avail), cur)
    if (!length(nbrs_all)) break
    nbrs_in_state   <- nbrs_all[state_of[nbrs_all] == home_state]
    nbrs_candidates <- if (length(nbrs_in_state)) nbrs_in_state else nbrs_all
    cand_totals <- total + overd_vec[nbrs_candidates]
    best <- if (all(cand_totals < min_overd)) nbrs_candidates[which.max(cand_totals)]
            else nbrs_candidates[which.min(abs(cand_totals - min_overd))]
    new_total <- total + overd_vec[best]
    if (is.na(new_total)) break
    allow_cross <- (!length(nbrs_in_state)) && (total < min_overd)
    if (new_total <= min_overd * 1.6 || total < min_overd || allow_cross) {
      cur   <- c(cur, best)
      total <- new_total
      avail <- setdiff(avail, best)
      if (total >= min_overd) break
    } else break
  }
  clustersD[cur] <- paste0("CL", cid)
  to_assign <- setdiff(to_assign, cur)
  cid <- cid + 1L
}
cluster_map_D <- tibble::tibble(county_ihme = names(clustersD), cluster_id  = unname(clustersD))

dq_cluster_D <- dq %>%
  dplyr::filter(!is.na(time_window)) %>%
  dplyr::inner_join(cluster_map_D, by = "county_ihme") %>%
  dplyr::group_by(cluster_id, time_window) %>%
  dplyr::summarise(
    pct_overd_miss = mean(tidyr::replace_na(pct_overd_miss, 0)),
    n_counties     = dplyr::n_distinct(county_ihme),
    .groups = "drop"
  ) %>%
  tidyr::complete(cluster_id, time_window = periods_vec_D,
                  fill = list(pct_overd_miss = 0, n_counties = 0L))

build_clusters_sf_D <- function(year) {
  base <- build_ihme_groups_sf_state(year)
  base %>%
    dplyr::left_join(cluster_map_D, by = "county_ihme") %>%
    dplyr::filter(!is.na(cluster_id)) %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(.groups = "drop") %>%
    sf::st_transform(crs_proj)
}
shapefile_years_D <- c("1999_2005"=2000, "2006_2012"=2010, "2013_2019"=2019, "2020_2022"=2020)
joined_by_period_D <- purrr::imap(shapefile_years_D, function(year, window) {
  shape_all <- build_clusters_sf_D(year)
  dplyr::left_join(shape_all, dplyr::filter(dq_cluster_D, time_window == window), by = "cluster_id")
})
make_map_D <- function(sf_data, var = "pct_overd_miss", title = NULL, palette = "Reds") {
  states <- tigris::states(cb = TRUE, class = "sf") |> sf::st_transform(crs_proj)
  fill_scale <- scale_fill_distiller(palette = palette, limits = c(0, 0.6), oob = scales::squish, direction = -1, na.value = "grey90")
  ggplot() +
    geom_sf(data = sf_data, aes(fill = .data[[var]]), colour = NA) +
    geom_sf(data = states, fill = NA, colour = "white", linewidth = 0.3) +
    fill_scale +
    coord_sf(xlim = c(-2500000, 2500000), ylim = c(-2200000, 730000), expand = FALSE) +
    labs(title = title %||% "", fill = NULL) +
    theme_void() +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 9))
}
p_1999_D <- make_map_D(joined_by_period_D[["1999_2005"]], "pct_overd_miss", "1999_2005")
p_2006_D <- make_map_D(joined_by_period_D[["2006_2012"]], "pct_overd_miss", "2006_2012")
p_2013_D <- make_map_D(joined_by_period_D[["2013_2019"]], "pct_overd_miss", "2013_2019")
p_2020_D <- make_map_D(joined_by_period_D[["2020_2022"]], "pct_overd_miss", "2020_2022")
combined_D <- (p_1999_D | p_2006_D) / (p_2013_D | p_2020_D)
ggsave(filename = file.path(here::here("figures","5yr_avg_stable"), "pct_overd_miss_4panel.png"),
       plot = combined_D, width = 12, height = 9, dpi = 320)


# ──────────────────────────────────────────────────────────────
# SECTION E — Z-scores & maps (overdose-unspecified, prop_garbage, RI, Philips)
# ──────────────────────────────────────────────────────────────
data_file       <- here::here("data",   "county_year_quality_metrics.csv")
membership_file <- here::here("output", "county_cluster_membership.csv.gz")
philips_file    <- here::here("output", "cluster_metrics_ucr39_cstd.csv.gz")
out_dir_figs    <- here::here("figures","5yr_avg_stable")
out_dir_scores  <- here::here("output", "cluster_scores")
dir.create(out_dir_figs,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_scores, recursive = TRUE, showWarnings = FALSE)

period_levels <- c("1999_2005","2006_2012","2013_2019","2020_2022")

# counties & states AK/HI shifted
build_counties_2163_resized <- function() {
  us_counties <- tigris::counties(year = 2020, cb = TRUE, class = "sf") |>
    sf::st_transform(crs_proj) |>
    dplyr::select(GEOID, STATEFP, geometry)
  mainland <- us_counties |> dplyr::filter(!STATEFP %in% c("02","15"))
  alaska   <- us_counties |> dplyr::filter(STATEFP == "02") |> st_scale(0.33) |> st_shift(c(1100000, -4700000)) |> sf::st_set_crs(crs_proj)
  hawaii   <- us_counties |> dplyr::filter(STATEFP == "15") |> st_scale(1.00)  |> st_shift(c(5000000, -1100000)) |> sf::st_set_crs(crs_proj)
  dplyr::bind_rows(mainland, alaska, hawaii) |> dplyr::rename(fips = GEOID, statefp = STATEFP) |> sf::st_make_valid()
}
build_states_2163_resized <- function() {
  us_states <- tigris::states(year = 2020, cb = TRUE, class = "sf") |>
    sf::st_transform(crs_proj) |>
    dplyr::select(STATEFP, geometry)
  mainland <- us_states |> dplyr::filter(!STATEFP %in% c("02","15"))
  alaska   <- us_states |> dplyr::filter(STATEFP == "02") |> st_scale(0.33) |> st_shift(c(1100000, -4700000)) |> sf::st_set_crs(crs_proj)
  hawaii   <- us_states |> dplyr::filter(STATEFP == "15") |> st_scale(1.00)  |> st_shift(c(5000000, -1100000)) |> sf::st_set_crs(crs_proj)
  dplyr::bind_rows(mainland, alaska, hawaii) |> sf::st_make_valid()
}
load(here::here("data_raw","ihme_fips.rda"))
ihme_map2 <- ihme_fips %>%
  dplyr::transmute(
    fips        = stringr::str_pad(as.character(orig_fips), 5, pad = "0"),
    county_ihme = stringr::str_pad(as.character(ihme_fips), 5, pad = "0")
  ) %>% dplyr::distinct()

membership_raw <- readr::read_csv(membership_file, show_col_types = FALSE)
membership <-
  if ("county_ihme" %in% names(membership_raw)) {
    membership_raw %>%
      dplyr::mutate(
        county_ihme = stringr::str_pad(as.character(county_ihme), 5, pad = "0"),
        period      = standardize_period(as.character(period)),
        cluster     = as.character(cluster)
      )
  } else if ("fips" %in% names(membership_raw)) {
    membership_raw %>%
      dplyr::mutate(
        fips   = stringr::str_pad(as.character(fips), 5, pad = "0"),
        period = standardize_period(as.character(period)),
        cluster= as.character(cluster)
      ) %>%
      dplyr::left_join(ihme_map2, by = "fips") %>%
      dplyr::mutate(county_ihme = dplyr::coalesce(county_ihme, fips)) %>%
      dplyr::select(county_ihme, period, cluster)
  } else stop("membership_file must contain 'county_ihme' or 'fips'.")

membership <- membership %>%
  dplyr::filter(!is.na(period), !is.na(county_ihme), !is.na(cluster)) %>%
  dplyr::distinct(county_ihme, period, cluster) %>%
  dplyr::group_by(county_ihme, period) %>%
  dplyr::summarise(cluster = mode_str(cluster), .groups="drop")

dy0 <- readr::read_csv(data_file, show_col_types = FALSE) %>%
  dplyr::mutate(year = suppressWarnings(as.integer(year)),
                period = period_of_year(year))
if ("county_ihme" %in% names(dy0)) {
  dy <- dy0 %>% dplyr::mutate(county_ihme = stringr::str_pad(as.character(county_ihme), 5, pad = "0"))
} else if ("fips" %in% names(dy0)) {
  dy <- dy0 %>% dplyr::mutate(fips = stringr::str_pad(as.character(fips), 5, pad = "0")) %>%
    dplyr::left_join(ihme_map2, by = "fips") %>%
    dplyr::mutate(county_ihme = dplyr::coalesce(county_ihme, fips))
} else stop("data_file must contain 'county_ihme' or 'fips'.")

nms <- names(dy)
col_n_cert <- pick_col(nms, c("n_cert","deaths","total","N"))
col_prop_g <- pick_col(nms, c("foreman_garbage","prop_garbage"))
if (!("N_garbage" %in% nms)) {
  if (!is.na(col_prop_g) && !is.na(col_n_cert)) {
    message("[E] Creating N_garbage = ", col_prop_g, " * ", col_n_cert)
    dy <- dy %>% dplyr::mutate(N_garbage = .data[[col_prop_g]] * .data[[col_n_cert]])
  } else if ("garb_k" %in% nms) {
    message("[E] Using garb_k as N_garbage"); dy <- dy %>% dplyr::mutate(N_garbage = garb_k)
  } else stop("[E] Missing N_garbage and cannot construct it.")
}
if (!("garb_k" %in% nms)) {
  if (!is.na(col_prop_g) && !is.na(col_n_cert)) {
    message("[E] Creating garb_k = ", col_prop_g, " * ", col_n_cert)
    dy <- dy %>% dplyr::mutate(garb_k = .data[[col_prop_g]] * .data[[col_n_cert]])
  } else if ("N_garbage" %in% names(dy)) {
    message("[E] Creating garb_k from N_garbage"); dy <- dy %>% dplyr::mutate(garb_k = N_garbage)
  } else stop("[E] Cannot create garb_k.")
}
need_vars <- c("county_ihme","period","overd_n","overd_miss_k","garb_k","n_cert","N_garbage","RI_post_only")
miss <- setdiff(need_vars, names(dy))
if (length(miss)) stop("[E] Missing required columns in data_file: ", paste(miss, collapse=", "))

cluster_core <- dy %>%
  dplyr::select(dplyr::all_of(need_vars)) %>%
  dplyr::filter(!is.na(period)) %>%
  dplyr::inner_join(membership, by = c("county_ihme","period")) %>%
  dplyr::group_by(period, cluster) %>%
  dplyr::summarise(
    sum_overd_miss_k = sum(overd_miss_k, na.rm = TRUE),
    sum_overd_n      = sum(overd_n,      na.rm = TRUE),
    sum_garb_k       = sum(garb_k,       na.rm = TRUE),
    sum_n_cert       = sum(n_cert,       na.rm = TRUE),
    sum_N_garbage    = sum(N_garbage,    na.rm = TRUE),
    num_RI_weighted  = sum(n_cert * RI_post_only, na.rm = TRUE),
    prop_overd_unspec = safe_div(sum_overd_miss_k, sum_overd_n),
    prop_garbage      = safe_div(sum_garb_k,       sum_n_cert),
    RI_cluster        = safe_div(num_RI_weighted,  sum_N_garbage),
    .groups = "drop"
  )

ph0 <- readr::read_csv(philips_file, show_col_types = FALSE)
names(ph0) <- tolower(names(ph0))
if (!("period" %in% names(ph0))) {
  if ("year" %in% names(ph0)) ph0 <- ph0 %>% dplyr::mutate(period = period_of_year(as.integer(year)))
  else stop("[E] Philips file needs 'period' or 'year'.")
} else {
  ph0 <- ph0 %>% dplyr::mutate(period = standardize_period(period))
}
if ("cluster" %in% names(ph0)) {
  ph_norm <- ph0 %>% dplyr::mutate(cluster = as.character(cluster))
} else if ("county_ihme" %in% names(ph0) || "fips" %in% names(ph0)) {
  tmp <- ph0
  if ("fips" %in% names(tmp) && !("county_ihme" %in% names(tmp))) {
    tmp <- tmp %>% dplyr::mutate(fips = stringr::str_pad(as.character(fips), 5, pad = "0")) %>%
      dplyr::left_join(ihme_map2, by = "fips")
  }
  if (!("county_ihme" %in% names(tmp))) stop("[E] Philips file lacks id columns (cluster/county_ihme/fips).")
  ph_norm <- tmp %>%
    dplyr::mutate(county_ihme = stringr::str_pad(as.character(county_ihme), 5, pad = "0")) %>%
    dplyr::inner_join(membership, by = c("county_ihme","period"))
} else stop("[E] Philips file must contain 'cluster' or county id (county_ihme/fips).")
ph_col <- pick_col(names(ph_norm), c("detail_ucod_icd4_cstd","cstd_ucr39","cstd","phillips_detail","philips_detail","cod_cstd","cod_cstd_ucr39"))
if (is.na(ph_col)) stop("[E] Philips source missing CSTD column.")
w_col  <- pick_col(names(ph_norm), c("n_cert","deaths","total","n"))
philips_clu <- ph_norm %>%
  dplyr::filter(!is.na(period)) %>%
  dplyr::mutate(cluster = as.character(cluster)) %>%
  dplyr::group_by(period, cluster) %>%
  dplyr::summarise(
    w   = if (!is.na(w_col)) sum(.data[[w_col]], na.rm=TRUE) else NA_real_,
    num = if (!is.na(w_col)) sum(.data[[w_col]] * .data[[ph_col]], na.rm=TRUE) else sum(.data[[ph_col]], na.rm=TRUE),
    den = if (!is.na(w_col)) w else sum(!is.na(.data[[ph_col]])),
    philips_cstd = safe_div(num, den),
    .groups = "drop"
  ) %>% dplyr::filter(!is.na(philips_cstd))

metrics_clu <- cluster_core %>%
  dplyr::inner_join(philips_clu, by = c("period","cluster")) %>%
  dplyr::filter(period %in% period_levels) %>%
  dplyr::arrange(period, cluster) %>%
  dplyr::mutate(cluster = as.character(cluster), period = as.character(period))
score_vars  <- c("prop_overd_unspec","philips_cstd","RI_cluster","prop_garbage")
global_means <- vapply(score_vars, function(v) mean(metrics_clu[[v]], na.rm = TRUE), numeric(1))
global_sds   <- vapply(score_vars, function(v)  sd(metrics_clu[[v]], na.rm = TRUE), numeric(1))
global_sds[!is.finite(global_sds) | global_sds == 0] <- NA_real_

z_tbl <- metrics_clu %>%
  dplyr::mutate(dplyr::across(
    dplyr::all_of(score_vars),
    \(x, nm=cur_column()) if (is.na(global_sds[[nm]])) NA_real_ else (x - global_means[[nm]]) / global_sds[[nm]],
    .names = "z_{.col}"
  )) %>%
  dplyr::mutate(
    z_prop_garbage      = -z_prop_garbage,
    z_prop_overd_unspec = -z_prop_overd_unspec
  ) %>%
  dplyr::mutate(direction_score = rowMeans(dplyr::across(starts_with("z_")), na.rm = TRUE)) %>%
  dplyr::arrange(period, cluster)

readr::write_csv(z_tbl, file.path(out_dir_scores, "cluster_zscores_overdose_philips_RI_propgarbage.csv"))

counties_tf <- build_counties_2163_resized()
states_tf   <- build_states_2163_resized()
counties_aug <- counties_tf %>% dplyr::left_join(ihme_map2, by = "fips") %>% dplyr::mutate(county_ihme = dplyr::coalesce(county_ihme, fips))

build_cluster_sf <- function(period_key, membership_df, counties_sf_aug) {
  mm <- membership_df %>% dplyr::filter(.data$period == period_key) %>% dplyr::transmute(county_ihme = as.character(county_ihme), cluster = as.character(cluster)) %>% dplyr::distinct()
  if (!nrow(mm)) stop("No membership rows for period: ", period_key)
  counties_sf_aug %>%
    dplyr::inner_join(mm, by = "county_ihme") %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
    sf::st_make_valid()
}
make_faceted_map <- function(stacked_sf, fill_col, lims, states_sf) {
  ggplot2::ggplot() +
    geom_sf(data = stacked_sf, aes(fill = .data[[fill_col]]), colour = NA) +
    geom_sf(data = states_sf, fill = NA, colour = "white", linewidth = 0.3) +
    scale_fill_distiller(palette = "RdBu", direction = 1, limits = lims, oob = scales::squish, na.value = "grey90") +
    coord_sf(xlim = c(-2500000, 2500000), ylim = c(-2200000, 730000), expand = FALSE) +
    labs(fill = NULL) + facet_wrap(~ period, ncol = 2) + theme_void() +
    theme(strip.text = element_text(size = 12, face = "bold"), legend.text = element_text(size = 9), plot.title = element_blank())
}
to_map <- c("z_RI_cluster","z_prop_overd_unspec","z_prop_garbage","z_philips_cstd","direction_score")
for (mcol in to_map) {
  message("[E] Assembling 4-panel for metric: ", mcol)
  stacked_sf <- purrr::map_dfr(period_levels, function(win){
    cl_sf <- build_cluster_sf(win, membership, counties_aug)
    dat   <- z_tbl %>% dplyr::filter(period==win) %>% dplyr::select(cluster, !!rlang::sym(mcol)) %>% dplyr::rename(value = !!rlang::sym(mcol))
    out   <- cl_sf %>% dplyr::left_join(dat, by="cluster"); out$period <- win; out
  }) %>% dplyr::mutate(period=factor(period, levels=period_levels))
  vals <- stacked_sf$value[is.finite(stacked_sf$value)]
  if (!length(vals) || length(unique(vals))==1) { message("  Skipping ", mcol, " (no variance)"); next }
  lims_sym <- compute_limits_symmetric(vals)
  p <- make_faceted_map(stacked_sf, "value", lims_sym, states_tf)
  base <- if (mcol == "direction_score") "direction_score" else sanitize(mcol)
  ggsave(file.path(out_dir_figs, paste0(base, "_4panel.png")), p, width=10, height=7.5, dpi=320)
}

# ──────────────────────────────────────────────────────────────
# SECTION F — Two maps + box plot (aggregate index) 
# ──────────────────────────────────────────────────────────────
metric_to_plot <- "direction_score"
metric_label <- dplyr::case_when(
  metric_to_plot == "direction_score"      ~ "Aggregate data quality index (z)",
  metric_to_plot == "z_prop_garbage"       ~ "Proportion garbage (z, higher is better)",
  metric_to_plot == "z_prop_overd_unspec"  ~ "Overdose unspecified (z, higher is better)",
  metric_to_plot == "z_philips_cstd"       ~ "Philips CSTD (z, higher is better)",
  metric_to_plot == "z_RI_cluster"         ~ "RI (z, higher is better)",
  TRUE ~ metric_to_plot
)

make_single_map_two <- function(period_key, mcol, lims, states_sf, membership_df, counties_sf_aug) {
  cl_sf <- build_cluster_sf(period_key, membership_df, counties_sf_aug)
  dat   <- z_tbl %>% dplyr::filter(period == period_key) %>% dplyr::select(cluster, value = !!rlang::sym(mcol))
  map_sf <- cl_sf %>% dplyr::left_join(dat, by = "cluster")
  ggplot2::ggplot() +
    geom_sf(data = map_sf, aes(fill = value), colour = NA) +
    geom_sf(data = states_sf, fill = NA, colour = "white", linewidth = 0.3) +
    scale_fill_distiller(palette = "RdBu", direction = 1, limits = lims, oob = scales::squish, na.value = "grey90") +
    coord_sf(xlim = c(-2500000, 2500000), ylim = c(-2200000, 730000), expand = FALSE) +
    labs(title = gsub("_", "–", period_key), fill = NULL) +
    theme_void(base_size = 11) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}
periods_two <- c("1999_2005","2020_2022")
map_vals <- purrr::map_dfr(periods_two, function(win){
  cl_sf <- build_cluster_sf(win, membership, counties_aug)
  dat   <- z_tbl %>% dplyr::filter(period == win) %>% dplyr::select(cluster, value = !!rlang::sym(metric_to_plot))
  cl_sf %>% dplyr::left_join(dat, by="cluster") %>% dplyr::transmute(value)
})
lims_sym <- compute_limits_symmetric(map_vals$value)
pA <- make_single_map_two("1999_2005", metric_to_plot, lims_sym, states_tf, membership, counties_aug)
pB <- make_single_map_two("2020_2022", metric_to_plot, lims_sym, states_tf, membership, counties_aug)

states_info <- tigris::states(year = 2020, cb = TRUE, class = "sf") |>
  sf::st_drop_geometry() |>
  dplyr::select(STATEFP, STUSPS) %>% dplyr::rename(statefp = STATEFP, state = STUSPS)
counties_lite <- counties_aug %>% sf::st_drop_geometry() %>% dplyr::select(county_ihme, statefp)
build_county_values <- function(per){
  mm <- membership %>% dplyr::filter(period == per) %>% dplyr::select(county_ihme, cluster)
  zv <- z_tbl %>% dplyr::filter(period == per) %>% dplyr::select(cluster, value = !!rlang::sym(metric_to_plot))
  mm %>% dplyr::left_join(zv, by = "cluster") %>% dplyr::left_join(counties_lite, by = "county_ihme") %>%
    dplyr::left_join(states_info, by = "statefp") %>% dplyr::mutate(period = per) %>%
    dplyr::filter(!is.na(state), is.finite(value)) %>% dplyr::select(state, county_ihme, period, value)
}
bx_df <- dplyr::bind_rows(build_county_values("1999_2005"), build_county_values("2020_2022"))
state_order <- bx_df %>%
  dplyr::filter(period == "2020_2022") %>% dplyr::group_by(state) %>%
  dplyr::summarise(med = median(value, na.rm=TRUE), .groups="drop") %>%
  dplyr::arrange(dplyr::desc(med)) %>% dplyr::pull(state)
territories <- c("PR","GU","VI","AS","MP")
bx_df <- bx_df %>% dplyr::filter(!state %in% territories) %>%
  dplyr::mutate(
    state = factor(state, levels = state_order),
    period = factor(period, levels = c("1999_2005","2020_2022"),
                    labels = c("1999–2005","2020–2022"))
  )
ylims_box <- c(-2, 2)
pC <- ggplot2::ggplot(bx_df, aes(x = state, y = value, fill = period)) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
  geom_boxplot(width = 0.7, outlier.alpha = 0.4,
               position = position_dodge2(width = 0.75, preserve = "single")) +
  coord_cartesian(ylim = ylims_box) +
  scale_fill_manual(values = c("#94bedf","#f2c879")) +
  labs(x = "States (ranked by declining median county-level 2020–2022 aggregate data quality index)",
       y = metric_label) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.minor = element_blank())
combined_two <- (pA | pB) / pC + patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(tag_levels = "A")
base <- if (metric_to_plot == "direction_score") "direction_score" else sanitize(metric_to_plot)
outfile <- file.path(out_dir_figs, paste0("two_maps_box_", base, ".png"))
ggsave(outfile, combined_two, width = 11, height = 8.5, dpi = 320)

# ──────────────────────────────────────────────────────────────
# SECTION F_RI — Figure 5: Two RI maps + box plot (1999–2005 vs 2020–2022)
# Depends on: dq, cluster_map, dq_cluster_RI, build_ihme_groups_sf(), crs_proj
# Output: figures/5yr_avg_stable/two_maps_box_RI.png
# ──────────────────────────────────────────────────────────────

out_dir_figs <- here::here("figures","5yr_avg_stable")
dir.create(out_dir_figs, recursive = TRUE, showWarnings = FALSE)

# Fallbacks if needed
if (!exists("states_outline")) {
    states_outline <- tigris::states(cb = TRUE, class = "sf") |>
        sf::st_transform(crs_proj) |>
        sf::st_make_valid()
}
if (!exists("counties_aug")) {
    # Build counties with AK/HI repositioned (same style as Section E)
    st_shift <- function(x, offset = c(0,0)) { sf::st_geometry(x) <- sf::st_geometry(x) + offset; x }
    st_scale <- function(x, factor = 1) { ctr <- sf::st_centroid(sf::st_union(x)); sf::st_geometry(x) <- (sf::st_geometry(x) - ctr) * factor + ctr; x }
    counties_tf <- tigris::counties(year = 2020, cb = TRUE, class = "sf") |>
        sf::st_transform(crs_proj) |>
        dplyr::select(GEOID, STATEFP, geometry)
    mainland <- counties_tf |> dplyr::filter(!STATEFP %in% c("02","15"))
    alaska   <- counties_tf |> dplyr::filter(STATEFP == "02") |> st_scale(0.33) |> st_shift(c(1100000, -4700000))
    hawaii   <- counties_tf |> dplyr::filter(STATEFP == "15") |> st_scale(1.00)  |> st_shift(c(5000000, -1100000))
    counties_tf <- dplyr::bind_rows(mainland, alaska, hawaii) |> dplyr::rename(fips = GEOID, statefp = STATEFP) |> sf::st_make_valid()
    if (!exists("ihme_map2")) {
        load(here::here("data_raw","ihme_fips.rda"))
        ihme_map2 <- ihme_fips %>%
            dplyr::transmute(
                fips        = stringr::str_pad(as.character(orig_fips), 5, pad = "0"),
                county_ihme = stringr::str_pad(as.character(ihme_fips), 5, pad = "0")
            ) %>% dplyr::distinct()
    }
    counties_aug <- counties_tf %>% dplyr::left_join(ihme_map2, by = "fips") %>%
        dplyr::mutate(county_ihme = dplyr::coalesce(county_ihme, fips))
}

# Cluster polygons (union counties by cluster_id from Section B)
clusters_aug_sf_RI <- counties_aug %>%
    dplyr::inner_join(cluster_map, by = "county_ihme") %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
    sf::st_make_valid()

# Keep only the two requested periods
periods_two_RI <- c("1999_2005","2020_2022")

# Map data: cluster-level average RI for the two periods
ri_cluster_two <- dq_cluster_post %>%
    dplyr::filter(time_window %in% periods_two_RI) %>%
    dplyr::select(cluster_id, time_window, RI_post_only)

# Choose sensible, fixed color limits for RI (adjust if needed)
lims_RI <- c(0.07, 0.12)

make_map_RI_two <- function(map_sf, title) {
    ggplot() +
        geom_sf(data = map_sf, aes(fill = RI_post_only), colour = NA) +
        geom_sf(data = states_outline, fill = NA, colour = "white", linewidth = 0.3) +
        scale_fill_distiller(palette = "Reds", limits = lims_RI,
                             oob = scales::squish, direction = -1, na.value = "grey90") +
        coord_sf(xlim = c(-2500000, 2500000), ylim = c(-2200000, 730000), expand = FALSE) +
        labs(title = title, fill = "RI") +
        theme_void(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

map_1999_RI <- clusters_aug_sf_RI %>%
    dplyr::left_join(dplyr::filter(ri_cluster_two, time_window == "1999_2005"), by = "cluster_id")
map_2020_RI <- clusters_aug_sf_RI %>%
    dplyr::left_join(dplyr::filter(ri_cluster_two, time_window == "2020_2022"), by = "cluster_id")

pA_RI <- make_map_RI_two(map_1999_RI, "1999–2005")
pB_RI <- make_map_RI_two(map_2020_RI, "2020–2022")

# ---- Box plot with the same two periods ----
# Copy cluster RI values down to member counties, then summarise by state
states_info_RI <- tigris::states(year = 2020, cb = TRUE, class = "sf") %>%
    sf::st_drop_geometry() %>%
    dplyr::select(STATEFP, STUSPS) %>%
    dplyr::rename(statefp = STATEFP, state = STUSPS) %>%
    dplyr::mutate(statefp = stringr::str_pad(as.character(statefp), 2, pad = "0"))

counties_lite_RI <- counties_aug %>%
    sf::st_drop_geometry() %>%
    dplyr::select(county_ihme, statefp)

build_county_values_post <- function(per) {
    zv <- ri_cluster_two %>%
        dplyr::filter(time_window == per) %>%
        dplyr::select(cluster_id, RI_post_only) %>%
        dplyr::rename(value = RI_post_only)
    
    cluster_map %>%
        dplyr::left_join(zv, by = "cluster_id") %>%
        dplyr::left_join(counties_lite_RI, by = "county_ihme") %>%
        dplyr::left_join(states_info_RI, by = "statefp") %>%
        dplyr::mutate(period = per) %>%
        dplyr::filter(!is.na(state), is.finite(value)) %>%
        dplyr::select(state, county_ihme, period, value)
}

bx_df_RI <- dplyr::bind_rows(
    build_county_values_post("1999_2005"),
    build_county_values_post("2020_2022")
)

# Order states by 2020–2022 median RI
territories_RI <- c("PR","GU","VI","AS","MP")
state_order_RI <- bx_df_RI %>%
    dplyr::filter(period == "2020_2022", !state %in% territories_RI) %>%
    dplyr::group_by(state) %>%
    dplyr::summarise(med = stats::median(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(med)) %>% dplyr::pull(state)

bx_df_RI <- bx_df_RI %>%
    dplyr::filter(!state %in% territories_RI) %>%
    dplyr::mutate(
        state  = factor(state, levels = state_order_RI),
        period = factor(period, levels = c("1999_2005","2020_2022"),
                        labels = c("1999–2005","2020–2022"))
    )

# Y-limits based on central 98% of data (robust to outliers)
ri_vals <- bx_df_RI$value[is.finite(bx_df_RI$value)]
ylims_box_RI <- as.numeric(stats::quantile(ri_vals, c(0.01, 0.99), na.rm = TRUE))
if (!all(is.finite(ylims_box_RI))) ylims_box_RI <- c(min(ri_vals, na.rm=TRUE), max(ri_vals, na.rm=TRUE))

pC_RI <- ggplot2::ggplot(bx_df_RI, aes(x = state, y = value, fill = period)) +
    geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
    geom_boxplot(width = 0.7, outlier.alpha = 0.4,
                 position = position_dodge2(width = 0.75, preserve = "single")) +
    coord_cartesian(ylim = ylims_box_RI) +
    labs(
        x = "States (ranked by declining median county-level 2020–2022 RI)",
        y = "Average RI (normalized K–L divergence)"
    ) +
    theme_bw(base_size = 10) +
    theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.minor = element_blank()
    )

# Combine (two maps on top, box below)
combined_RI_two <- (pA_RI | pB_RI) / pC_RI +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A")

outfile_RI_two <- file.path(out_dir_figs, "two_maps_box_RI.png")
ggsave(outfile_RI_two, combined_RI_two, width = 11, height = 8.5, dpi = 320)
message("✓ [Figure 5] Two RI maps + box plot (1999–2005 vs 2020–2022) saved to: ", outfile_RI_two)


# ──────────────────────────────────────────────────────────────
# SECTION G — Two maps + box plot: pct_overd_miss (cluster_map from B)
# ──────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────
# Build `rep_lookup` (no auto-detect)
# Requires one CSV with county FIPS + reporting system type
# ──────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(readr); library(dplyr); library(stringr); library(here) })

# If you don't already have these:
if (!exists("std_fips"))  std_fips  <- function(x) stringr::str_pad(as.character(x), 5, pad = "0")
if (!exists("mode_val")) mode_val <- function(x){ ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

# 1) Point to your file (pick ONE and set explicitly)
reporting_file <- here::here("data_raw","County-Death-Investigation-System-2018-1-9-2024.csv")
# reporting_file <- here::here("data_raw","County-Death-Investigation-System-2018.csv")

# 2) Explicit column names in that file (EDIT these to match your headers)
fips_col <- "FIPS_CODE"           
# --------------------------
# Robust read of reporting file and build rep_lookup (replacement for original block)
# --------------------------
# File is already set above as `reporting_file`
fips_col <- "FIPS_CODE"   # stays the same

# Read file (don't assume a particular "type" column name)
raw_rep <- readr::read_csv(reporting_file, show_col_types = FALSE)

# Helper to find the first column matching any of a set of regex patterns (case-insensitive)
find_col <- function(df, patterns) {
    nm <- names(df)
    i <- which(vapply(patterns, function(p) any(stringr::str_detect(nm, regex(p, ignore_case = TRUE))), logical(1)))
    if (length(i) == 0L) return(NA_character_)
    # return the first matching column name from the NAME vector (in order of patterns)
    for (p in patterns) {
        match_idx <- which(stringr::str_detect(nm, regex(p, ignore_case = TRUE)))
        if (length(match_idx)) return(nm[match_idx[1]])
    }
    NA_character_
}

# Patterns for likely useful columns (ordered by priority)
# 1) direct 'Medicolegal Death Investigation Type' or variants -> contains 'Coroner' / 'Medical Examiner' (preferred)
# 2) 'Type of Coroner Office' -> e.g. 'Elected Coroner'
# 3) 'Medical Examiner System' -> may contain 'Non-Medical Examiner ...' (use lower priority)
# 4) 'Other Official System' -> fallback
# 5) 'State Medical Examiner and County Official' -> fallback
col_medicolegal <- find_col(raw_rep, c("Medicolegal|Medicolegal_Death|death investigation type"))
col_coroner_type <- find_col(raw_rep, c("Type of Coroner|Type_of_Coroner|coroner type"))
col_me_system <- find_col(raw_rep, c("^Medical Examiner System$|Medical_Examiner|medical examiner system"))
col_other_system <- find_col(raw_rep, c("Other Official System|Other_Official|other official"))
col_state_med <- find_col(raw_rep, c("State Medical Examiner|State_Medical|state medical"))

# Compose reporting_type0 by taking the first non-missing, non-empty value from the priority list
raw_rep2 <- raw_rep %>%
    dplyr::mutate(across(everything(), ~ as.character(.x))) %>%
    dplyr::mutate(
        # helper: NA if value is empty string or only whitespace
        across(dplyr::all_of(intersect(names(.), c(col_medicolegal, col_coroner_type, col_me_system, col_other_system, col_state_med))),
               ~ ifelse(is.na(.x) | str_squish(.x) == "", NA_character_, .x))
    )

reporting_source_vec <- c(col_medicolegal, col_coroner_type, col_me_system, col_other_system, col_state_med)
# keep only the ones that exist
reporting_source_vec <- reporting_source_vec[!is.na(reporting_source_vec)]

# Build rep_lookup
rep_lookup <- raw_rep2 %>%
    transmute(
        county_ihme_raw = .data[[fips_col]],
        # pull a reporting string using the first available (priority order)
        reporting_type0 = purrr::pmap_chr(
            list(!!!rlang::syms(reporting_source_vec)),
            function(...) {
                vals <- c(...)
                # return first non-NA, non-empty
                v <- vals[!is.na(vals) & vals != ""]
                if (length(v)) return(v[[1]])
                return(NA_character_)
            }
        )
    ) %>%
    # standardize fips
    dplyr::mutate(county_ihme = std_fips(county_ihme_raw),
                  reporting_type0 = stringr::str_squish(as.character(reporting_type0))) %>%
    # normalize labels into the four classes used elsewhere
    dplyr::mutate(
        reporting_type = dplyr::case_when(
            # exact/common coroner indicators
            stringr::str_detect(stringr::str_to_lower(reporting_type0), "\\bcoroner\\b") ~ "Coroner",
            # mixed/hybrid indicators
            stringr::str_detect(stringr::str_to_lower(reporting_type0), "mixed|hybrid")  ~ "Mixed",
            # "other" indicator
            stringr::str_detect(stringr::str_to_lower(reporting_type0), "\\bother\\b") ~ "Other County Official",
            # If the priority column explicitly says "Medical Examiner" (avoid matching 'Non-Medical Examiner' incorrectly)
            stringr::str_detect(stringr::str_to_lower(reporting_type0), "\\bmedical examiner\\b|\\bmedical_examiner\\b") ~ "Medical Examiner",
            # fallback: if the string contains 'examiner' but also 'non' immediately before, treat as NOT medical examiner.
            TRUE ~ stringr::str_to_title(stringr::str_squish(reporting_type0))
        )
    ) %>%
    # remove empty / bad fips
    dplyr::filter(!is.na(county_ihme), nchar(county_ihme) == 5, !is.na(reporting_type)) %>%
    # if multiple rows per fips keep modal reporting type
    dplyr::group_by(county_ihme) %>%
    dplyr::summarise(reporting_type = mode_val(reporting_type), .groups = "drop") %>%
    dplyr::mutate(
        reporting_type = factor(reporting_type,
                                levels = c("Coroner","Other County Official","Mixed","Medical Examiner"))
    )

# quick sanity check
rep_lookup %>% dplyr::count(reporting_type, sort = TRUE) %>% print(n = 99)

periods_two_G <- c("1999_2005","2020_2022")
lims_map_G <- c(0, 0.60)

# Build a robust ihme_map if it doesn't exist yet
if (!exists("ihme_map") || !all(c("GEOID","county_ihme") %in% names(ihme_map))) {
    if (exists("ihme_fips")) {
        ihme_map <- ihme_fips %>%
            dplyr::transmute(
                GEOID       = stringr::str_pad(as.character(orig_fips), 5, pad = "0"),
                county_ihme = stringr::str_pad(as.character(ihme_fips), 5, pad = "0")
            ) %>%
            dplyr::distinct()
    } else {
        ihme_map <- tibble::tibble(GEOID = character(), county_ihme = character())
    }
}

# Ensure counties_aug has fips + statefp + county_ihme
counties_aug_fixed <- counties_aug
if ("fips" %in% names(counties_aug_fixed)) {
    counties_aug_fixed <- counties_aug_fixed %>%
        dplyr::mutate(fips = stringr::str_pad(as.character(fips), 5, pad = "0"))
} else if ("GEOID" %in% names(counties_aug_fixed)) {
    counties_aug_fixed <- counties_aug_fixed %>%
        dplyr::mutate(fips = stringr::str_pad(as.character(GEOID), 5, pad = "0"))
} else {
    stop("`counties_aug` must have either 'fips' or 'GEOID'.")
}
if (!"statefp" %in% names(counties_aug_fixed)) {
    if ("STATEFP" %in% names(counties_aug_fixed)) {
        counties_aug_fixed <- counties_aug_fixed %>% dplyr::mutate(statefp = as.character(STATEFP))
    } else stop("`counties_aug` must have 'statefp' or 'STATEFP'.")
}
counties_aug_fixed <- counties_aug_fixed %>%
    dplyr::mutate(statefp = stringr::str_pad(as.character(statefp), 2, pad = "0")) %>%
    dplyr::left_join(ihme_map, by = c("fips" = "GEOID")) %>%
    { 
        if (!"county_ihme" %in% names(.)) dplyr::mutate(., county_ihme = .$fips)
        else dplyr::mutate(., county_ihme = dplyr::coalesce(.data$county_ihme, .data$fips))
    }

# Cluster polygons (union counties by cluster_id from Section B)
clusters_aug_sf <- counties_aug_fixed %>%
    dplyr::inner_join(cluster_map, by = "county_ihme") %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
    sf::st_make_valid()

# Cluster-level averages of pct_overd_miss for the two target periods
all_clusters <- sort(unique(cluster_map$cluster_id))
pct_cluster <- dq %>%
    dplyr::filter(!is.na(time_window)) %>%
    dplyr::select(county_ihme, time_window, pct_overd_miss) %>%
    dplyr::inner_join(cluster_map, by = "county_ihme") %>%
    dplyr::group_by(cluster_id, time_window) %>%
    dplyr::summarise(pct_overd_miss = mean(pct_overd_miss, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(cluster_id = all_clusters, time_window = periods_two_G,
                    fill = list(pct_overd_miss = 0))

map_1999 <- clusters_aug_sf %>%
    dplyr::left_join(dplyr::filter(pct_cluster, time_window == "1999_2005"),
                     by = "cluster_id") %>%
    dplyr::mutate(pct_overd_miss = dplyr::coalesce(pct_overd_miss, 0))

map_2020 <- clusters_aug_sf %>%
    dplyr::left_join(dplyr::filter(pct_cluster, time_window == "2020_2022"),
                     by = "cluster_id") %>%
    dplyr::mutate(pct_overd_miss = dplyr::coalesce(pct_overd_miss, 0))

make_map_G <- function(map_sf, title, lims) {
    ggplot() +
        geom_sf(data = map_sf, aes(fill = pct_overd_miss), colour = NA) +
        geom_sf(data = states_tf, fill = NA, colour = "white", linewidth = 0.3) +
        scale_fill_distiller(palette = "Reds", limits = lims, oob = scales::squish,
                             na.value = "grey90", direction = 1,
                             labels = scales::percent_format(accuracy = 1)) +
        coord_sf(xlim = c(-2500000, 2500000), ylim = c(-2200000, 730000), expand = FALSE) +
        labs(title = title, fill = NULL) +
        theme_void(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

pA_G <- make_map_G(map_1999, "1999–2005", lims_map_G)
pB_G <- make_map_G(map_2020, "2020–2022", lims_map_G)

# Build county-level values (for boxplot) by copying cluster values to member counties
states_info_G <- tigris::states(year = 2020, cb = TRUE, class = "sf") %>%
    sf::st_drop_geometry() %>%
    dplyr::select(STATEFP, STUSPS) %>%
    dplyr::rename(statefp = STATEFP, state = STUSPS) %>%
    dplyr::mutate(statefp = stringr::str_pad(as.character(statefp), 2, pad = "0"))

counties_lite_G <- counties_aug_fixed %>%
    sf::st_drop_geometry() %>%
    dplyr::select(county_ihme, statefp)

build_county_values_G <- function(per) {
    zv <- pct_cluster %>%
        dplyr::filter(time_window == per) %>%
        dplyr::select(cluster_id, value = pct_overd_miss)
    
    cluster_map %>%
        dplyr::left_join(zv, by = "cluster_id") %>%
        dplyr::left_join(counties_lite_G, by = "county_ihme") %>%
        dplyr::left_join(states_info_G, by = "statefp") %>%
        dplyr::mutate(period = per) %>%
        dplyr::filter(!is.na(state), is.finite(value)) %>%
        dplyr::select(state, county_ihme, period, value)
}

bx_df_G <- dplyr::bind_rows(
    build_county_values_G("1999_2005"),
    build_county_values_G("2020_2022")
)

# Order states by 2020–2022 median
state_order_G <- bx_df_G %>%
    dplyr::filter(period == "2020_2022") %>%
    dplyr::group_by(state) %>%
    dplyr::summarise(med = stats::median(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(med)) %>%
    dplyr::pull(state)

territories_G <- c("PR","GU","VI","AS","MP")
bx_df_G <- bx_df_G %>%
    dplyr::filter(!state %in% territories_G) %>%
    dplyr::mutate(
        state  = factor(state, levels = state_order_G),
        period = factor(period,
                        levels = c("1999_2005","2020_2022"),
                        labels = c("1999–2005","2020–2022"))
    )

pC_G <- ggplot2::ggplot(bx_df_G, aes(x = state, y = value, fill = period)) +
    geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
    geom_boxplot(width = 0.7, outlier.alpha = 0.4,
                 position = position_dodge2(width = 0.75, preserve = "single")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
        x = "States (ranked by declining median county-level 2020–2022 % overdose unspecified)",
        y = "% overdose unspecified"
    ) +
    theme_bw(base_size = 10) +
    theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.minor = element_blank()
    )

combined_G <- (pA_G | pB_G) / pC_G +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A")

outfile_G <- file.path(out_dir_figs, "two_maps_box_pct_overd_miss.png")
ggsave(outfile_G, combined_G, width = 11, height = 8.5, dpi = 320)
message("✓ [G] Two maps + box plot for pct_overd_miss saved to: ", outfile_G)

# ──────────────────────────────────────────────────────────────
# SECTION G2 — Basic summary stats for Results text
# ──────────────────────────────────────────────────────────────

out_dir_tab  <- here::here("output")
dir.create(out_dir_tab, recursive = TRUE, showWarnings = FALSE)

# --- Helper for median (IQR) by period
summarise_metric <- function(df, metric, periods = c("1999_2005","2020_2022")) {
    df %>%
        dplyr::filter(period %in% periods, is.finite(.data[[metric]])) %>%
        dplyr::group_by(period) %>%
        dplyr::summarise(
            median = median(.data[[metric]], na.rm = TRUE),
            p25    = quantile(.data[[metric]], 0.25, na.rm = TRUE),
            p75    = quantile(.data[[metric]], 0.75, na.rm = TRUE),
            n      = dplyr::n(),
            .groups = "drop"
        )
}

# Use z_tbl (cluster-level metrics) for stable comparisons
summary_RI      <- summarise_metric(z_tbl, "RI_cluster")
summary_garbage <- summarise_metric(z_tbl, "prop_garbage")
summary_overd   <- summarise_metric(z_tbl, "prop_overd_unspec")
summary_philips <- summarise_metric(z_tbl, "philips_cstd")

# --- Correlations between metrics (county-level, latest period)
bx_latest <- bx_df %>% dplyr::filter(period == "2020–2022")

cor_RI_garbage <- cor(z_tbl$RI_cluster, z_tbl$prop_garbage, use="pairwise", method="pearson")
cor_RI_overd   <- cor(z_tbl$RI_cluster, z_tbl$prop_overd_unspec, use="pairwise", method="pearson")
cor_g_overd    <- cor(z_tbl$prop_garbage, z_tbl$prop_overd_unspec, use="pairwise", method="pearson")
cor_overd_phil <- cor(z_tbl$prop_overd_unspec, z_tbl$philips_cstd, use="pairwise", method="pearson")
summary_dir <- summarise_metric(z_tbl, "direction_score")

# --- Export to CSV for insertion into manuscript
summary_all <- list(
    RI        = summary_RI,
    garbage   = summary_garbage,
    overd_uns = summary_overd,
    philips   = summary_philips,
    aggregate = summary_dir 
)

purrr::iwalk(summary_all, ~ readr::write_csv(.x,
                                             file.path(out_dir_tab, paste0("summary_", .y, ".csv"))))

# Save correlations
cor_tbl <- tibble::tibble(
    pair = c("RI vs garbage","RI vs overd_unspec","garbage vs overd_unspec","overd_unspec vs philips"),
    rho  = c(cor_RI_garbage, cor_RI_overd, cor_g_overd, cor_overd_phil)
)
readr::write_csv(cor_tbl, file.path(out_dir_tab, "correlations.csv"))

message("✓ [G2] Summary tables written to: ", out_dir_tab)

# ──────────────────────────────────────────────────────────────
# SECTION G3 — Numbers for two specific manuscript claims
#   1) "Median garbage coding increased from G (IQR) in 2013–2019 to G′ (IQR′) in 2020–2022"
#   2) "ME vs Coroner aggregate index in 2020–2022 (median diff Δ; ME: M (IQR), Coroner: C (IQR))"
# ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(stringr) })

# ---- guards ----
stopifnot(exists("z_tbl"), exists("membership"), exists("rep_lookup"))
if (!all(c("period","cluster") %in% names(z_tbl))) stop("[G3] z_tbl must include period and cluster.")
if (!("prop_garbage" %in% names(z_tbl))) stop("[G3] z_tbl must include prop_garbage (cluster-level).")
if (!("direction_score" %in% names(z_tbl))) stop("[G3] z_tbl must include direction_score (cluster-level).")
if (!all(c("county_ihme","period","cluster") %in% names(membership))) stop("[G3] membership missing required cols.")
if (!all(c("county_ihme","reporting_type") %in% names(rep_lookup))) stop("[G3] rep_lookup missing required cols.")

out_dir_tab <- here::here("output")
dir.create(out_dir_tab, recursive = TRUE, showWarnings = FALSE)

# helpers
std_fips <- function(x) stringr::str_pad(as.character(x), 5, pad = "0")
period_pretty <- function(p) dplyr::recode(p,
                                           "2013_2019"="2013–2019", "2020_2022"="2020–2022",
                                           "1999_2005"="1999–2005", "2006_2012"="2006–2012", .default = p
)
summarise_m_iqr <- function(df, value_col, group_col) {
    df %>%
        dplyr::filter(is.finite(.data[[value_col]])) %>%
        dplyr::group_by(.data[[group_col]]) %>%
        dplyr::summarise(
            median = stats::median(.data[[value_col]], na.rm = TRUE),
            p25    = stats::quantile(.data[[value_col]], 0.25, na.rm = TRUE),
            p75    = stats::quantile(.data[[value_col]], 0.75, na.rm = TRUE),
            n      = dplyr::n(),
            .groups = "drop"
        ) %>%
        dplyr::rename(!!group_col := 1)
}
