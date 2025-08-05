library(tidyverse)
library(here)
library(fs)
library(arrow)
library(future)
library(furrr)

wonder_df <- read_csv(here("data_raw", "Multiple Cause of Death, 1999-2020.csv"),
    col_names = c("note", "county_name", "county_fips", "nchs_deaths", "population", "crude_rate")
) |>
    filter(is.na(note)) |>
    select(-note)

future::plan(future::multisession(workers = availableCores() - 2))
deaths_county_year <- future_map_dfr(
    .x = dir_ls(here("data_private", "mcod"), glob = "*.parquet"),
    .f = ~ {
        read_parquet(.x) |>
            group_by(county_ihme, year) |>
            summarize(n_deaths = n()) |>
            ungroup()
    }
)
future::plan(future::sequential())

summarized_county <- deaths_county_year |>
    arrange(county_ihme, year) |>
    group_by(county_ihme) |>
    mutate(norm_deaths = (n_deaths - mean(n_deaths)) / sd(n_deaths)) |>
    mutate(rel_change = (n_deaths - lag(n_deaths)) / lag(n_deaths)) |>
    ungroup()

ggplot(
    summarized_county |>
        filter(year < 2020),
    aes(
        y = county_ihme,
        x = year,
        fill = norm_deaths
    )
) +
    geom_tile() +
    scale_y_discrete(NULL) +
    scale_fill_viridis_c() +
    theme(axis.text.y = element_blank())

summarized_county |>
    slice_max(abs(norm_deaths), n = 20)

future::plan(future::multisession(workers = availableCores() - 2))
deaths_county <- future_map_dfr(
    .x = dir_ls(here("data_private", "mcod"), glob = "*.parquet"),
    .f = ~ {
        read_parquet(.x) |>
            group_by(county_fips, year) |>
            summarize(n_deaths = n()) |>
            ungroup()
    }
)
future::plan(future::sequential())

deaths_county <- deaths_county |>
    filter(year <= 2020) |>
    group_by(county_fips) |>
    summarize(n_deaths = sum(n_deaths)) |>
    left_join(wonder_df |>
        mutate(nchs_deaths = as.integer(nchs_deaths))) |>
    mutate(rel_diff = (n_deaths - nchs_deaths) / nchs_deaths)

ggplot(
    deaths_county,
    aes(
        x = nchs_deaths,
        y = n_deaths,
        size = rel_diff
    )
) +
    geom_point(alpha = .5) +
    geom_abline(yintercept = 0, slope = 1, alpha = .2) +
    coord_equal() +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    scale_size_area() +
    theme_bw()
