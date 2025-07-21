import pandas as pd
import re

# Load the raw ICD-10 file
file_path = "/Users/amymann/Documents/Data Quality Project/DataQuality/data_raw/cause-codes/icd10cm-codes-April-2025.txt"

# Read lines manually (since it's fixed-format text, not CSV)
with open(file_path, "r", encoding="utf-8") as f:
    lines = f.readlines()

# Parse into (code, description) pairs
records = []
for line in lines:
    if len(line.strip()) < 8:
        continue
    code = line[:7].strip()
    desc = line[8:].strip().lower()
    records.append((code, desc))

df = pd.DataFrame(records, columns=["code", "description"])

# Garbage keywords to include
include_keywords = [
    "unspecified", "ill-defined", "not otherwise specified"
    "not elsewhere classified", "unknown"
]
include_keywords = [
     "ill-defined", "not otherwise specified"
    "not elsewhere classified", "unknown"
]

# More specific phrases we may want to exclude
exclude_keywords = ["specific", "due to", "with", "secondary to"]

# Create regex masks
inc_mask = df["description"].str.contains("|".join(include_keywords), na=False, regex=True)
exc_mask = df["description"].str.contains("|".join(exclude_keywords), na=False, regex=True)

# Final filtered light garbage codes
light_garbage_df = df[inc_mask & ~exc_mask].copy()

# Optional: manually add high-priority known garbage codes
manual_add = ["I64", "I50.9", "J18.9", "R99"]
manual_df = df[df["code"].isin(manual_add)]
light_garbage_df = pd.concat([light_garbage_df, manual_df]).drop_duplicates()

# Export to CSV
output_path = "/Users/amymann/Documents/Data Quality Project/DataQuality/data_raw/cause-codes/full_light_garbage_codes.csv"
light_garbage_df.to_csv(output_path, index=False)

output_path
