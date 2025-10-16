# merge_sheets_to_csv.R
# =====================
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")

library(readxl)
library(dplyr)
library(stringr)
library(here)

xlsx_path <- here("data_raw", "fanewcsareas2.xlsx")

# === Read sheet names ===
sheets <- excel_sheets(xlsx_path)

dfs <- lapply(sheets, function(sh) {
    df <- read_excel(xlsx_path, sheet = sh, col_types = "text")
    state_name <- sub(" \\(.*\\)$", "", sh)
    df$State <- state_name
    
    df[] <- lapply(df, function(x) {
        if (is.factor(x)) as.character(x) else as.character(x)
    })

    if ("FIPS State Code" %in% names(df) || "FIPS County Code" %in% names(df)) {
        # Ensure columns exist as character (create if missing)
        if (!("FIPS State Code" %in% names(df))) df[["FIPS State Code"]] <- NA_character_
        if (!("FIPS County Code" %in% names(df))) df[["FIPS County Code"]] <- NA_character_
        df[["FIPS State Code"]] <- str_trim(df[["FIPS State Code"]])
        df[["FIPS County Code"]] <- str_trim(df[["FIPS County Code"]])
        df[["FIPS State Code"]][df[["FIPS State Code"]] == ""] <- NA_character_
        df[["FIPS County Code"]][df[["FIPS County Code"]] == ""] <- NA_character_
        
        df <- df %>%
            mutate(
                `FIPS State Code` = str_pad(str_trim(`FIPS State Code`), width = 2, side = "left", pad = "0"),
                `FIPS County Code` = str_pad(str_trim(`FIPS County Code`), width = 3, side = "left", pad = "0"),
                FIPS = if_else(
                    !is.na(`FIPS State Code`) & !is.na(`FIPS County Code`),
                    paste0(`FIPS State Code`, `FIPS County Code`),
                    NA_character_
                )
            )
    } else {
        df$FIPS <- NA_character_
    }
    
    cols_now <- names(df)
    front <- intersect(c("State", "FIPS"), cols_now)
    df <- dplyr::select(df, all_of(front), everything())
    
    df
})

combined <- dplyr::bind_rows(dfs)

output_csv <- here("data_raw", "all_county_groupings.csv")
write.csv(combined, file = output_csv, row.names = FALSE, na = "")