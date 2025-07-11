"""
Overwrite mort21999–mort22023.parquet with a clean `race5` column.

 • Source column: 1999–2020 ⇒ 'race', 2021–2023 ⇒ 'racer5'
 • Hispanic codes: 1–5  (≤2002), 200–299 (≥2003)
 • Removes *every* pre-existing race5 field before writing
 • Works even if the file already contained duplicate race5 columns
"""

from pathlib import Path
import pyarrow as pa, pyarrow.compute as pc
import pyarrow.parquet as pq
import pyarrow.dataset as ds
import os

DATA_DIR = Path("/Users/amymann/Documents/Data Quality Project/data/parquet")

# ── helper: robust str → int32 -----------------------------------------------
# ── robust str/float/int → int32 ------------------------------------
def to_int32(col: pa.ChunkedArray | pa.Array) -> pa.ChunkedArray:
    """Return an int32 array. Handles integer, floating, and string columns."""
    if pa.types.is_integer(col.type):
        return pc.cast(col, pa.int32(), safe=False)

    # handle numeric floats like 1.0, 99.0, NaN
    if pa.types.is_floating(col.type):
        # turn NaN → null before cast
        clean = pc.if_else(pc.is_nan(col),
                           pa.scalar(None, col.type),
                           col)
        return pc.cast(clean, pa.int32(), safe=False)

    # otherwise treat as string/binary ------------------------------------------------
    col_bin = pc.cast(col, pa.binary(), safe=False)           # safe for utf8-cleaning
    for bad in (b"\x00", b"nan", b"NAN", b"na", b"NA"):
        col_bin = pc.replace_substring(col_bin, bad, b"")
    col_utf8 = pc.cast(col_bin, pa.utf8(), safe=False)

    digits = pc.match_substring_regex(col_utf8, r"^\d+$")
    clean  = pc.if_else(digits, col_utf8,
                        pa.scalar(None, pa.string()))
    return pc.cast(clean, pa.int32(), safe=False)


# ── build race5 ---------------------------------------------------------------
def make_race5(batch, race_col, year):
    race = to_int32(batch.column(race_col))
    hisp = to_int32(batch.column("hispanic"))

    if year <= 2002:
        is_hisp = pc.is_in(hisp, pa.array([1,2,3,4,5], pa.int32()))
    else:
        is_hisp = pc.and_(pc.greater_equal(hisp, 200), pc.less_equal(hisp, 299))

    is_unknown = pc.or_(pc.is_null(race), pc.is_null(hisp))
    is_white, is_black, is_native = [pc.equal(race, v) for v in (1, 2, 3)]

    return pc.if_else(
        is_unknown, pa.scalar(None, pa.string()),
        pc.if_else(
            is_hisp, pa.scalar("Hispanic"),
            pc.if_else(
                is_white, pa.scalar("White"),
                pc.if_else(
                    is_black, pa.scalar("Black"),
                    pc.if_else(
                        is_native, pa.scalar("American Indian"),
                        pa.scalar("Other")
                    )
                )
            )
        )
    )

# ── main loop -----------------------------------------------------------------
for year in range(1999, 2024):
    src = DATA_DIR / f"mort{year}.parquet"
    if not src.exists():
        print(f"⚠️  {src.name} missing — skipped"); continue

    race_col = "race" if year <= 2020 else "racer5"
    print(f"→ {src.name}   (race source = {race_col})")

    tmp = src.with_suffix(".race5tmp.parquet")

    try:
        # 1.  original file-level schema
        file_schema = pq.read_schema(src)

        # 2.  drop *all* race5 fields (there may be duplicates)
        base_fields = [f for f in file_schema if f.name != "race5"]
        base_schema = pa.schema(base_fields)

        # 3.  open dataset *with the cleaned schema*  ➜ bypass duplicate error
        ds_year = ds.dataset(src, format="parquet", schema=base_schema)
        scan    = ds_year.scanner(batch_size=200_000)

        # 4.  final output schema = base + new race5
        out_schema = base_schema.append(pa.field("race5", pa.string()))
        writer     = pq.ParquetWriter(tmp, out_schema, compression="snappy")

        for b in scan.to_batches():
            cols = {n: c for n, c in zip(b.schema.names, b.columns) if n != "race5"}
            cols["race5"] = make_race5(b, race_col, year)
            writer.write_table(pa.Table.from_pydict(cols, schema=out_schema))

        writer.close()
        os.replace(tmp, src)               # overwrite in place atomically
        print("   ✓ overwritten")

    except Exception as err:
        if tmp.exists():
            tmp.unlink(missing_ok=True)
        print(f"   ✗ ERROR: {err}")
