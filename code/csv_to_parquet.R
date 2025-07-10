# Code to convert data from .csv form to .parquet form which uses less storage 
# and is much faster to run. Code assumes data is in format mortYYYY.csv

library(arrow)
library(stringr)

path = "/Users/amymann/Documents/Data Quality Project/data/"
setwd(path)

csv_files <- list.files(pattern = "^mort\\d{4}\\.csv$")

for (csv_file in csv_files) {
  year <- str_extract(csv_file, "\\d{4}")
    parquet_file <- paste0("mort", year, ".parquet")
  
  if (file.exists(parquet_file)) {
    message("Skipping year ", year, ": parquet already exists.")
    next
  }
  
  message("Processing year ", year, "...")
  df <- readr::read_csv(csv_file, show_col_types = FALSE)
  write_parquet(df, parquet_file)
}