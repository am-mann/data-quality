# Adds descriptions of codes from ICD10 to the garbage code list 

library(tidyverse)
library(icd)

garbage <- read_csv("/Users/amymann/Documents/Data Quality Project/data_quality/gbd_garbage_codes.csv")

chapter_fallback <- tribble(
  ~prefix, ~fallback,
  "U",   "U-codes reserved / national use",
  "X59", "Exposure to unspecified factor",
  "Y99", "External cause – activity unspecified",
  "W76", "Other accidental hanging and strangulation",
  "V99", "Transport accident – unspecified"
)

label_icd <- function(code) {
  code <- toupper(str_trim(code))
  
  # (1) exact
  exact <- icd::explain_code(code, short_desc = TRUE)
  if (!is.na(exact[1])) return(exact[1])
  
  # (2) 3-character category
  cat3  <- icd::explain_code(substr(code, 1, 3), short_desc = TRUE)
  if (!is.na(cat3[1]))   return(cat3[1])
  
  # (3) fallback table (chapter or custom)
  fb <- chapter_fallback$fallback[match(TRUE, startsWith(code, chapter_fallback$prefix))]
  if (!is.na(fb[1])) return(fb[1])
  
  NA_character_
}

garbage <- garbage %>%
  mutate(description = map_chr(icd10, label_icd))

#  Inspect remaining NAs
garbage %>% filter(is.na(description)) %>% count(icd10)

write_csv(garbage, "gbd_garbage_codes_with_descr.csv")