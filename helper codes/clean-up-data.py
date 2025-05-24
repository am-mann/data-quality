from pathlib import Path
import pyarrow as pa, pyarrow.compute as pc
import pyarrow.parquet as pq, pyarrow.dataset as ds
import os

DATA_DIR = Path("/Users/amymann/Documents/Data Quality Project/data/parquet")

# ─────────────────────────────────────────────────────────────────────────────
def to_int32(col: pa.ChunkedArray | pa.Array) -> pa.ChunkedArray:
    if pa.types.is_integer(col.type):
        return pc.cast(col, pa.int32(), safe=False)
    if pa.types.is_floating(col.type):
        clean = pc.if_else(pc.is_nan(col), pa.scalar(None, col.type), col)
        return pc.cast(clean, pa.int32(), safe=False)
    col_bin = pc.cast(col, pa.binary(), safe=False)
    for bad in (b"\x00", b"nan", b"NAN", b"na", b"NA"):
        col_bin = pc.replace_substring(col_bin, bad, b"")
    col_utf8 = pc.cast(col_bin, pa.utf8(), safe=False)
    digits   = pc.match_substring_regex(col_utf8, r"^\d+$")
    clean    = pc.if_else(digits, col_utf8, pa.scalar(None, pa.string()))
    return pc.cast(clean, pa.int32(), safe=False)

# ─────────────────────────────────────────────────────────────────────────────
def map_codes(raw, mapping, *, default=None):
    # cast any input (int / bytes / utf8) to utf8
    raw_utf8 = pc.cast(raw, pa.utf8(), safe=False)

    keys   = list(mapping.keys())
    values = [mapping[k] for k in keys]

    idx = pc.index_in(raw_utf8, pa.array([str(k) for k in keys]))
    out = pc.take(pa.array(values, pa.string()), idx)

    if default is not None:
        out = pc.if_else(pc.equal(idx, -1),
                         pa.scalar(default, pa.string()),
                         out)

    return pc.cast(out, pa.string(), safe=False)

# ─────────────────────────────────────────────────────────────────────────────
SEX_MAP = {1: "M", 2: "F", "M": "M", "F": "F"}
MARSTAT_MAP = {1: "S", 2: "M", 3: "W", 4: "D", 8: "U", 9: "U",
               "S": "S", "M": "M", "W": "W", "D": "D", "U": "U"}
PLACDTH_MAP = {
    1: "Hospital inpatient",
    2: "Hospital outpatient/ER",
    3: "Hospital DOA",
    4: "Home",
    5: "Hospice facility",
    6: "Nursing home/LTC",
    7: "Other",
    9: "Unknown"
}

def clean_sex(raw, *_):     return map_codes(raw, SEX_MAP)
def clean_marstat(raw, *_): return map_codes(raw, MARSTAT_MAP)
def clean_placdth(raw):     return map_codes(raw, PLACDTH_MAP)

# ─────────────────────────────────────────────────────────────────────────────
for year in range(1999, 2024):
    src = DATA_DIR / f"mort{year}.parquet"
    if not src.exists():
        print(f"⚠️  {src.name} missing — skipped"); continue
    print(f"→ {src.name}")

    tmp = src.with_suffix(".fixdemo.tmp.parquet")
    try:
        schema0 = pq.read_schema(src)

        new_fields = [pa.field(f.name, pa.string())
                      if f.name in {"sex", "marstat", "placdth"} else f
                      for f in schema0]
        out_schema = pa.schema(new_fields)

        ds_year = ds.dataset(src, format="parquet", schema=schema0)
        writer  = pq.ParquetWriter(tmp, out_schema, compression="snappy")

        for b in ds_year.to_batches(batch_size=200_000):
            cols = dict(zip(b.schema.names, b.columns))
            cols["sex"]     = clean_sex(cols["sex"])
            cols["marstat"] = clean_marstat(cols["marstat"])
            cols["placdth"] = clean_placdth(cols["placdth"])
            writer.write_table(pa.Table.from_pydict(cols, schema=out_schema))

        writer.close()
        os.replace(tmp, src)
        print("   ✓ overwritten")

    except Exception as e:
        if tmp.exists():
            tmp.unlink(missing_ok=True)
        print(f"   ✗ ERROR: {e}")
