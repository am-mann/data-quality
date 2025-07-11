## 01_ingest_death_data.R ----
##
## Reads in the NCHS restricted access raw text files, subsets columns of
## interest, and then saves as parquets in `data_private`.
##
## NOTE: I use the {narcan} fixed width dictionaries rather than the whole
## package just to minimize dependencies. The file is `mcod_fwf_dicts.rda` and
## saved under `./data_raw` but can be downloaded on the {narcan} repo at:
##
## https://github.com/mkiang/narcan/blob/master/data/mcod_fwf_dicts.rda

## Imports ----
library(tidyverse)
library(here)
library(fs)
library(arrow)

## Constants ----
RAW_DIR <- here("data_raw", "mcod")
YEARS <- 1999:2022
SAVE_DIR <- here("data_private", "mcod")

## Data ----
load(here("data_raw", "mcod_fwf_dicts.rda"))

for (y in YEARS) {
    save_path <- here(SAVE_DIR, sprintf("mcod_%i.parquet", y))
    
    if (file_exists(save_path)) {
        next
    }
    
    ## Get file path and column positions
    f_path <- dir_ls(RAW_DIR,
        type = "file",
        regexp = sprintf("\\<MULT%i.*\\.zip\\>", y))
    fwf_col_pos <- mcod_fwf_dicts |>
        dplyr::filter(year == y) |>
        dplyr::select("start", "end", col_names = "name")
    c_types <- mcod_fwf_dicts |>
        dplyr::filter(year == y) |>
        dplyr::pull("type") |>
        paste(collapse = "")

    ## Read in
    temp_df <- readr::read_fwf(
        file = f_path,
        col_positions = readr::fwf_positions(
            start = fwf_col_pos$start,
            end = fwf_col_pos$end,
            col_names = fwf_col_pos$col_names
        ),
        col_types = c_types,
        na = c("", "NA", " ")
    )
    
    ## Clean up the columns that are inconsistently coded across years
    ## From 1999 to 2002, countyrs used NCHS coding instead of FIPS, so
    ## switch it back
    if (y %in% 1999:2002) {
        temp_df <- temp_df %>%
            rename(countyrs_old = countyrs) %>%
            rename(countyrs = fipsctyr)
    }
    
    ## Fix education
    ## In 2003, education switched to the 2003 system but retained a column
    ## called edu89. In 2014, educ89 was dropped. 
    if (y %in% 2003:2013) {
        temp_df <- temp_df |> 
            rename(educ_old = educ) |> 
            rename(educ = educ89)
    }

    ## Fix marital category
    temp_df <- temp_df %>%
        dplyr::mutate(
            marstat = dplyr::case_when(
                marstat == 1 ~ "single_never_married",
                marstat == 2 ~ "married",
                marstat == 3 ~ "widowed",
                marstat == 4 ~ "divorced",
                marstat == 8 ~ "not_on_certificate",
                marstat == 9 ~ "not_stated",
                is.na(marstat) ~ "unknown",

                marstat == "S" ~ "single_never_married",
                marstat == "M" ~ "married",
                marstat == "W" ~ "widowed",
                marstat == "D" ~ "divorced",
                marstat == "N" ~ "not_on_certificate",
                marstat == "U" ~ "not_stated",
                is.na(marstat) ~ "unknown",

                TRUE ~ as.character(marstat)
            )
        )

    ## Recode detailed age age_years
    ## Also make a variable called age_years, which will be more consistent
    if (y %in% 1989:2002) {
        temp_df <- temp_df %>%
            dplyr::mutate(age_years = case_when(
                age %in% c(299, 399, 499, 599, 699, 999) ~ NA_real_,
                age > 199 ~ 0,
                age <= 199 ~ age,
                TRUE ~ NA_real_
            ))
    } else {
        temp_df <- temp_df %>%
            dplyr::mutate(age_years = case_when(
                age %in% c(9999, 1999, 2999, 4999, 5999, 6999) ~ NA_real_,
                age > 1999 ~ 0,
                age < 1999 ~ age - 1000,
                TRUE ~ NA_real_
            ))
    }

    ## Fix sex category
    temp_df <- temp_df %>%
        dplyr::mutate(
            sex = dplyr::case_when(
                sex == 1 ~ "male",
                sex == "M" ~ "male",
                sex == 2 ~ "female",
                sex == "F" ~ "female",
                sex == "Male" ~ "male",
                sex == "Female" ~ "female",
                TRUE ~ as.character(sex)
            )
        )
    
    ## Subset to columns that match Amy's columns
    ## NOTES: 
    ##  - Race_Recode_40 is racer40 and doesn't exist for all years 
    ##  - occupational codes are new so not included. 
    ##  - educ changes over time but is renamed above.
    ##  - methdisp is new so not included
    ##  - autopsy is new so not included
    ##  - race, brace, racer3, racer5 stopped after 2020 I think
    ##  - similarly, hispanic was not recoded after 2020 so hispanicr is gone
    ##  - raceimp is not consistent across all years
    ##  - ageflag, ager22 are not across all years
    ##  - ucr130 not consistent across all years 
    ##  - we don't use the econd columns so I drop them
    temp_df <- temp_df |>
        select(
            restatus,
            countyrs, 
            educ,
            monthdth,
            sex,
            age_years,
            age,
            ager52,
            ager27,
            ager12,
            placdth,
            marstat,
            weekday,
            year,
            injwork,
            mandeath,
            activity,
            injury,
            ucod,
            ucr358,
            ucr113,
            ucr39,
            # eanum,
            # starts_with("econd"),
            ranum,
            starts_with("record_"),
            race,
            racer3,
            any_of(c("racer40")), 
            hispanic, 
            hspanicr
        )
    
    ## Save
    dir_create(dirname(save_path))
    write_parquet(temp_df,
                  save_path, 
                  compression = "snappy",
                  use_dictionary = TRUE, 
                  write_statistics = TRUE)
}

## Make some fake debugging files ----
for (f in dir_ls(SAVE_DIR, type = "file", glob = "*.parquet")) {
    temp_df <- read_parquet(f)
    
    sub_df <- temp_df |> 
        sample_frac(.2)
    
    sub_df$countyrs <- sample(temp_df$countyrs, NROW(sub_df), replace = TRUE)
    
    dir_create(here("data_private", "mcod_sample"))
    write_parquet(sub_df, here("data_private", "mcod_sample", basename(f)))
}
