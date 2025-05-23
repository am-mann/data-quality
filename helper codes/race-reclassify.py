#!/usr/bin/env python3
"""
Add `race5` to mort21999 … mort22023.

Logic
------
race source:
    1999-2020 : column 'race'
    2021-2023 : column 'racer5'  (single-race codes 1-6)

race5 mapping (Hispanic takes precedence):
    hispanic 1-5      -> Hispanic
    race 1            -> White
    race 2            -> Black
    race 3            -> American Indian
    everything else   -> Other
    NULL (either race or hispanic missing) -> NaN in parquet
"""

from pathlib import Path
import pyarrow as pa, pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq

DATA_DIR = Path("/Users/amymann/Documents/Data Quality Project/data/parquet")

# ── helper: robust string → int32  (same as before) ──────────────────
# ── robust string → int32  (Arrow-agnostic) ──────────────────────────
def to_int32(col):
    """
    Cast an Arrow column to int32.
    • Works whether the source is int, string, binary, or mixed.
    • Invalid UTF-8, NUL bytes, 'nan', blanks → null.
    """
    if pa.types.is_integer(col.type):
        return pc.cast(col, pa.int32(), safe=False)

    # 1)  Get raw bytes so we can clean them first
    col_bin = pc.cast(col, pa.binary(), safe=False)

    # 2)  Strip embedded NULs *and* the literal bytes for 'nan'/'NA'
    #     (case-insensitive); operates safely on binary type
    for bad in (b"\x00", b"nan", b"NAN", b"na", b"NA"):
        col_bin = pc.replace_substring(col_bin, bad, b"")

    # 3)  Now try to read as UTF-8.  Any remaining invalid sequence
    #     will raise, so we wrap in a try/except and fall back to null.
    try:
        col_utf8 = pc.cast(col_bin, pa.utf8(), safe=False)
    except pa.ArrowInvalid:
        # Rare path: some chunks still contain bad bytes.
        # Replace the whole offending slot with null.
        # (ChunkedArray → iterate chunk-wise.)
        new_chunks = []
        for chunk in col_bin.chunks:
            try:
                new_chunks.append(pc.cast(chunk, pa.utf8(), safe=False))
            except pa.ArrowInvalid:
                new_chunks.append(pa.array([None] * len(chunk),
                                            type=pa.string()))
        col_utf8 = pa.chunked_array(new_chunks, type=pa.string())

    # 4)  Keep only all-digit strings → else null
    digits = pc.match_substring_regex(col_utf8, r"^\d+$")
    clean  = pc.if_else(digits, col_utf8,
                        pa.scalar(None, type=pa.string()))

    # 5)  Final cast to int32
    return pc.cast(clean, pa.int32(), safe=False)

# ── build race5 for a RecordBatch ────────────────────────────────────
def make_race5(batch, race_col):
    race = to_int32(batch.column(race_col))
    hisp = to_int32(batch.column("hispanic"))

    is_unknown  = pc.or_(pc.is_null(race), pc.is_null(hisp))
    is_hispanic = pc.is_in(hisp, pa.array([1,2,3,4,5], pa.int32()))
    is_white, is_black, is_native = [pc.equal(race, v) for v in (1,2,3)]

    return pc.if_else(
        is_unknown, pa.scalar(None, pa.string()),          # NaN
        pc.if_else(
            is_hispanic, pa.scalar("Hispanic"),
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

# ── main loop 1999-2023 ─────────────────────────────────────────────
for year in range(1999, 2024):
    src = DATA_DIR / f"mort{year}.parquet"
    if not src.exists():
        print(f"⚠️  {src.name} missing — skipped")
        continue

    race_col = "race" if year <= 2020 else "racer5"
    print(f"\n→ Processing {src.name}  (race source = {race_col})")

    tmp    = src.with_name(src.stem + "_race5tmp.parquet")
    backup = src.with_name(src.stem + "_backup.parquet")

    try:
        ds_year = ds.dataset(src, format="parquet")
        if race_col not in ds_year.schema.names:
            raise KeyError(f"Column {race_col!r} not found in {src.name}")

        scan = ds_year.scanner(batch_size=200_000)
        writer = pq.ParquetWriter(
            tmp,
            ds_year.schema.append(pa.field("race5", pa.string())),
            compression="snappy"
        )

        for b in scan.to_batches():
            rb = pa.record_batch(
                list(b.columns) + [make_race5(b, race_col)],
                names=b.schema.names + ["race5"]
            )
            writer.write_batch(rb)
        writer.close()

        src.rename(backup)
        tmp.rename(src)
        print(f"   ✓ updated; backup → {backup.name}")

    except Exception as e:
        if tmp.exists():
            tmp.unlink()
        print(f"   ✗ ERROR in {src.name}: {e}")
