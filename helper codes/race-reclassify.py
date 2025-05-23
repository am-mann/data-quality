#!/usr/bin/env python3
"""
Update mort2YYYY.parquet for YYYY = 1999…2023
Adds `race5`:
    • Hispanic (hispanic 1-5)
    • White    (race 1)
    • Black    (race 2)
    • American Indian (race 3)
    • Other    (all remaining codes)
    • NULL     (either race or hispanic is null → shows up as NaN)
Each processed file is backed up as mort2YYYY_backup.parquet.
"""

from pathlib import Path
import pyarrow as pa, pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq

DATA_DIR = Path("/Users/amymann/Documents/Data Quality Project/data/parquet")

def race5_array(batch: pa.RecordBatch) -> pa.Array:
    race = pc.cast(batch.column("race"), pa.int32())
    hisp = pc.cast(batch.column("hispanic"), pa.int32())

    is_unknown  = pc.or_(pc.is_null(race), pc.is_null(hisp))
    is_hispanic = pc.is_in(hisp, pa.array([1,2,3,4,5], type=pa.int32()))
    is_white    = pc.equal(race, 1)
    is_black    = pc.equal(race, 2)
    is_native   = pc.equal(race, 3)

    return pc.if_else(
        is_unknown, pa.scalar(None, type=pa.string()),
        pc.if_else(
            is_hispanic, pa.scalar("Hispanic"),
            pc.if_else(
                is_white,  pa.scalar("White"),
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

for year in range(1999, 2024):
    src = DATA_DIR / f"mort{year}.parquet"
    if not src.exists():
        print(f"⚠️  {src.name} not found – skipping")
        continue

    print(f"\n→ Updating {src.name}")
    tmp    = src.with_name(src.stem + "_race5tmp.parquet")
    backup = src.with_name(src.stem + "_backup.parquet")

    try:
        ds_year = ds.dataset(src, format="parquet")
        scan    = ds_year.scanner(batch_size=200_000)

        out_schema = ds_year.schema.append(pa.field("race5", pa.string()))
        writer = pq.ParquetWriter(tmp, out_schema, compression="snappy")

        for b in scan.to_batches():
            writer.write_batch(
                pa.record_batch(
                    list(b.columns) + [race5_array(b)],
                    names=b.schema.names + ["race5"]
                )
            )
        writer.close()

        src.rename(backup)
        tmp.rename(src)
        print(f"   ✓ {src.name} done (backup: {backup.name})")

    except Exception as err:
        if tmp.exists():
            tmp.unlink()
        print(f"   ✗ ERROR in {src.name}: {err}")

