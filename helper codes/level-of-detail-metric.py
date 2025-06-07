import os
from pathlib import Path
from io import StringIO
import re
import math
import pandas as pd
import openai
import tiktoken   # counts tokens for any GPT-4* model

# ── 1 · Paths ────────────────────────────────────────────────────────────────
DATA_DIR  = Path("/Users/amymann/Documents/Data Quality Project/data_quality/")
SRC_FILE  = DATA_DIR / "all_secondary_by_ucod_2021.csv"
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")

# ── 2 · Read & compress the data before sending ─────────────────────────────
df = pd.read_csv(SRC_FILE, dtype=str)

# if the file already has a 'ucod' and many secondary columns ---------------
sec_cols = [c for c in df.columns if c != "ucod"]

# collapse into:  ucod | code1|code2|code3|...
agg = (
    df.melt(id_vars="ucod", value_vars=sec_cols, value_name="icd10")
      .dropna(subset=["icd10"])
      .query("icd10 != ''")
      .groupby("ucod")["icd10"]
      .apply(lambda s: "|".join(sorted(set(s))))
      .reset_index()
)

compressed_csv = agg.to_csv(index=False)
print(f"Compressed CSV size: {len(compressed_csv)/1024:.1f} KB")

# ── 3 · Helper to count tokens so we don’t exceed limits ────────────────────
enc = tiktoken.encoding_for_model("gpt-4.1")  # same rates as gpt-4o
def n_tokens(txt: str) -> int:
    return len(enc.encode(txt))

# Split the big mapping into ≤ 6000-token chunks (input+output well < 30k/min)
CHUNK_SIZE = 6000
lines = compressed_csv.splitlines()
header, rows = lines[0], lines[1:]

chunks = []
current = [header]
for r in rows:
    current.append(r)
    if n_tokens("\n".join(current)) > CHUNK_SIZE:
        # drop the last row (too much) and start a new chunk
        current.pop()
        chunks.append("\n".join(current))
        current = [header, r]
chunks.append("\n".join(current))

print(f"Prompt will be sent in {len(chunks)} chunk(s).")

# ── 4 · OpenAI client ────────────────────────────────────────────────────────
client = openai.OpenAI(api_key=OPENAI_API_KEY)

system_prompt = (
    "You are a medical coding and public-health expert in ICD-10 mortality data. "
    "For each chunk you receive, perform three tasks:\n"
    
    "Step 1: For each underlying cause of death (UCOD), identify contributing cause ICD-10 codes that provide more "
    "specific information but do not indicate background illness or general morbidity. For example, obesity should not "
    "be considered useful. Output a CSV with one row per UCOD, and separate columns for each valid contributing cause.\n\n"

    "Step 2: Make a list of 'light garbage' contributing causes — codes that are less specific than alternatives. "
    "For example, T40.9 is less informative than T40.1 or T40.6. Output a single-column CSV: 'light_garbage_codes'.\n\n"

    "Step 3: For each UCOD, list any contributing causes that are redundant, inappropriate, or misleading. "
    "These might be improperly coded or not meaningful in that context. Improper codes refer to when there is a duplication of information"
    "but there is NO harm to data quality. If there is a non-specific code that provides no additional information, this "
    "is light garbage and NOT an improper code is there exists a more specific code that would provide more information. "
    "Output a CSV with one row per UCOD, "
    "and separate columns for each improper code.\n\n"

    "STEP 1 CSV\n<raw CSV>\nSTEP 2 CSV\n<raw CSV>\nSTEP 3 CSV\n<raw CSV>\n"
    "No markdown, no commentary."
)

STEP1, STEP2, STEP3 = [], [], []

for idx, chunk in enumerate(chunks, 1):
    user_prompt = f"{chunk}\n\n(Reminder: output ONLY the three raw CSV blocks.)"
    print(f"Sending chunk {idx}/{len(chunks)} – {n_tokens(user_prompt)} tokens")

    resp = client.chat.completions.create(
        model="gpt-4o",
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt}
        ],
        temperature=0,
        max_tokens=4096,   # generous but safe
    )

    text = resp.choices[0].message.content

    # ── 5 · Extract the three CSVs from the returned text ────────────────────
    def grab(step_tag: str) -> str | None:
        pat = rf"{step_tag}\n(.*?)(?=(STEP \d CSV|$))"
        m = re.search(pat, text, re.DOTALL)
        return m.group(1).strip() if m else None

    STEP1.append(grab("STEP 1 CSV"))
    STEP2.append(grab("STEP 2 CSV"))
    STEP3.append(grab("STEP 3 CSV"))

# ── 6 · Concatenate & save the final results ────────────────────────────────
def save_csv(blocks, fname):
    joined = "\n".join(b for b in blocks if b)
    out_df = pd.read_csv(StringIO(joined))
    out_df.to_csv(DATA_DIR / fname, index=False)
    print(f"Saved {fname} ({len(out_df):,} rows).")

save_csv(STEP1, "step1_output.csv")
save_csv(STEP2, "step2_output.csv")
save_csv(STEP3, "step3_output.csv")