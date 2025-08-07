# ---- PERIODIC CLUSTERING + SIMPSON DIVERSITY ON ucr358 (WITH MERGES) ----

# 0. Packages and paths -----------------------------------------------------
pkgs <- c("sf","dplyr","arrow","here","igraph","stringr","purrr",
          "tigris","readr","data.table","fs")
for(p in pkgs) if(!requireNamespace(p,quietly=TRUE)) install.packages(p)
lapply(pkgs, library, character.only=TRUE)

parquet_dir <- if (fs::dir_exists(here("data_private","mcod"))) {
    here("data_private","mcod")
} else {
    here("data_private","mcod_sample")
}
county_var <- "county_ihme"
min_deaths  <- 200

# 1. Precompute county adjacency graph --------------------------------------
counties_sf <- tigris::counties(cb = TRUE, year = 2020) %>% select(GEOID)
adj <- sf::st_touches(counties_sf)
edge_df <- tibble(
    from = rep(counties_sf$GEOID, lengths(adj)),
    to   = counties_sf$GEOID[ unlist(adj) ]
) %>% filter(from < to)
g <- graph_from_data_frame(edge_df, directed = FALSE,
                           vertices = counties_sf$GEOID)

# 2. Greedy clustering (initial pass) ---------------------------------------
make_greedy_clusters <- function(death_tbl, g, min_deaths) {
    deaths    <- setNames(death_tbl$deaths, death_tbl$fips)
    to_assign <- names(deaths)
    clusters  <- setNames(rep(NA_character_, length(deaths)), to_assign)
    cid <- 1L
    
    while(length(to_assign) > 0) {
        this <- to_assign[which.max(deaths[to_assign])]
        if(deaths[this] >= min_deaths) {
            clusters[this] <- paste0("C", cid)
            to_assign <- setdiff(to_assign, this)
            cid <- cid + 1L
            next
        }
        
        cluster <- this
        total   <- deaths[this]
        avail   <- setdiff(to_assign, this)
        
        repeat {
            nbrs <- unique(unlist(lapply(cluster, function(v)
                neighbors(g, v, mode="all") %>% names())))
            nbrs <- intersect(setdiff(nbrs, cluster), avail)
            if(length(nbrs)==0) break
            
            gains <- deaths[nbrs]
            diffs <- abs((total + gains) - min_deaths)
            best  <- nbrs[which.min(diffs)]
            new_tot <- total + deaths[best]
            
            if(new_tot <= min_deaths * 1.5 || total < min_deaths) {
                cluster <- c(cluster, best)
                total   <- new_tot
                avail   <- setdiff(avail, best)
            } else break
        }
        
        clusters[cluster] <- paste0("C", cid)
        to_assign <- setdiff(to_assign, cluster)
        cid <- cid + 1L
    }
    
    clusters
}

# --- Improved post-hoc merge: guarantees no cluster < min_deaths ---
merge_small_clusters <- function(clu, death_tbl, g, min_deaths = 200) {
    deaths_vec <- setNames(death_tbl$deaths, death_tbl$fips)
    
    cluster_size <- function(cl) sum(deaths_vec[names(clu)[clu == cl]])
    
    repeat {
        # recompute cluster sizes each iteration
        sums   <- tapply(deaths_vec, clu, sum, na.rm = TRUE)
        small  <- names(sums)[sums < min_deaths]
        if (length(small) == 0) break         # all clusters are big enough
        
        for (sc in small) {
            members <- names(clu)[clu == sc]
            
            ## 1 · look for *valid* neighbour clusters ----------------------------
            neigh_clusters <- unique(unlist(lapply(members, function(f) {
                nb <- neighbors(g, f, mode = "all")          # igraph vertex sequence
                nb_names <- names(nb)
                unique(clu[nb_names])
            })))
            
            neigh_clusters <- setdiff(neigh_clusters, sc)          # exclude itself
            neigh_clusters <- neigh_clusters[neigh_clusters %in% names(sums)]  # valid
            
            if (length(neigh_clusters) == 0) {
                # 2 · no spatial neighbour left → merge into overall biggest cluster
                target <- names(sums)[which.max(sums)]
            } else {
                # 3 · merge into the largest neighbouring cluster
                target <- neigh_clusters[which.max(sums[neigh_clusters])]
            }
            
            clu[members] <- target
        }
        # loop again in case new merges created fresh <200 clusters
    }
    
    clu
}

# 4. Simpson diversity function ---------------------------------------------
simpson_diversity <- function(codes) {
    tab <- as.numeric(table(codes))
    p   <- tab / sum(tab)
    1 - sum(p^2)
}

# 5. Define periods ---------------------------------------------------------
periods <- list(
    "1999_2004" = 1999:2004,
    "2005_2010" = 2005:2010,
    "2011_2017" = 2011:2017,
    "2018_2022" = 2018:2022
)

# 6. Loop over periods ------------------------------------------------------
ds <- open_dataset(parquet_dir)
all_county_div  <- list()
all_cluster_div <- list()

for(pname in names(periods)) {
    yrs <- periods[[pname]]
    
    # 6A · load & prep
    cert <- ds %>%
        filter(year %in% yrs, !is.na(ucr358)) %>%
        select(county_ihme = !!county_var, ucr358) %>%
        collect() %>%
        mutate(
            fips   = str_pad(as.character(county_ihme), 5, pad = "0"),
            ucr358 = as.character(ucr358)
        ) %>%
        filter(fips %in% counties_sf$GEOID)
    
    # 6B · county death counts
    death_tbl <- cert %>%
        group_by(fips) %>%
        summarise(deaths = n(), .groups = "drop")
    
    # 6C · initial greedy clusters
    clu <- make_greedy_clusters(death_tbl, g, min_deaths)
    # 6D · post‐hoc merges
    clu <- merge_small_clusters(clu, death_tbl, g, min_deaths)
    
    cert <- cert %>% mutate(cluster = clu[fips])
    
    # 6E · county‐level Simpson
    cd <- cert %>%
        group_by(fips) %>%
        summarise(
            period  = pname,
            deaths  = n(),
            simpson = simpson_diversity(ucr358),
            .groups = "drop"
        )
    all_county_div[[pname]] <- cd
    
    # 6F · cluster‐level Simpson
    cld <- cert %>%
        group_by(cluster) %>%
        summarise(
            period  = pname,
            deaths  = n(),
            simpson = simpson_diversity(ucr358),
            .groups = "drop"
        )
    all_cluster_div[[pname]] <- cld
}

# 7. Bind & export ----------------------------------------------------------
all_county_div  <- bind_rows(all_county_div)
all_cluster_div <- bind_rows(all_cluster_div)

print(head(all_county_div,10))
print(head(all_cluster_div,10))

# Optionally:
write_csv(all_county_div,  "all_county_simpson.csv")
write_csv(all_cluster_div, "all_cluster_simpson.csv")
