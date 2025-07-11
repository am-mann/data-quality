import pandas as pd
import pyarrow.parquet as pq
import pyarrow as pa

# === 1. Define widths and column names ===

widths = [19,1,40,2,1,1,2,2,1,4,1,2,2,2,2,1,1,1,16,4,1,1,1,
          1,34,1,1,4,3,1,3,3,2,1,2,7,7,7,7,7,7,7,7,7,7,7,7,
          7,7,7,7,7,7,7,7,36,2,1,5,5,5,5,5,5,5,5,5,5,5,5,5,
          5,5,5,5,5,5,5,1,2,1,1,1,1,33,3,1,1,2,315,4,2,4,2]

names = ["blank1", "restatus", "blank2", "educ1989", "educ2003", "educflag",
         "monthdth", "blank3", "sex", "age", "ageflag", "ager52",
         "ager27", "ager12", "ager22", "placdth", "marstat",
         "weekday", "blank4", "year", "injwork", "mandeath", "methdisp",
         "autopsy", "blank5", "activity", "injury", "ucod", "ucr358", "blank6",
         "ucr113", "ucr130", "ucr39", "blank7", "eanum"] + \
        [f"econdp_{i+1}" for i in range(20)] + \
        ["blank8", "ranum", "blank9"] + \
        [f"record_{i+1}" for i in range(20)] + \
        ["blank10", "race", "brace", "raceimp", "racer3", "racer5", "blank11",
         "hispanic", "blank12", "hspanicr", "Race_Recode_40", "blank13",
         "CensusOcc", "Occ_26", "CensusInd", "Ind_23"]

# === 2. Set file paths ===

txt_path = "/Users/amymann/Documents/Data Quality Project/data/txt/VS19MORT.DUSMCPUB_r20210304"
parquet_path = "/Users/amymann/Documents/Data Quality Project/data/mort2019.parquet"

# === 3. Determine non-blank columns ===

cols_to_keep = [i for i, n in enumerate(names) if not n.startswith("blank")]

# === 4. Read file in chunks, drop blanks, save to parquet ===

arrow_tables = []
chunksize = 50000  # you can adjust this if you run into memory issues

for chunk in pd.read_fwf(
    txt_path,
    widths=widths,
    names=names,
    chunksize=chunksize,
    dtype=str):
    chunk = chunk.iloc[:, cols_to_keep]
    # Force all columns to be string, even if all are NA in this chunk
    for col in chunk.columns:
        chunk[col] = chunk[col].astype(str)
    arrow_tables.append(pa.Table.from_pandas(chunk, preserve_index=False))
    print(f"Processed {len(chunk)} rows...")

# Concatenate all Arrow tables
full_table = pa.concat_tables(arrow_tables)

# Write to a single Parquet file
pq.write_table(full_table, parquet_path)

print(f"Saved as {parquet_path}")
