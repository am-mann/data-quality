"""
Prerequisites
-------------
pip install --upgrade openai>=1.14

Environment
-----------
export OPENAI_API_KEY="sk-..."
"""

import os
import openai
import pandas as pd

accident_list = [
    "V01", "V02", "V03", "V04", "V05", "V06", "V09", "V10", "V11", "V12", "V13", "V14", "V19", "V20",
    "V24", "V28", "V30", "V34", "V38", "V40", "V44", "V48", "V50", "V54", "V58", "V60", "V64", "V68",
    "V70", "V74", "V78", "V80", "V82", "V84", "V86", "V88", "V90", "V92", "V94", "V96", "V98",
    "W00", "W01", "W02", "W03", "W04", "W05", "W06", "W07", "W08", "W09", "W10", "W11", "W12", "W13",
    "W14", "W15", "W16", "W17", "W18", "W19", "W20", "W21", "W22", "W23", "W24", "W25", "W26", "W27",
    "W28", "W29", "W30", "W31", "W32", "W34", "W35", "W36", "W37", "W38", "W39", "W41", "W42", "W43",
    "W44", "W45", "W46", "W47", "W48", "W49", "W50", "W51", "W52", "W53", "W54", "W55", "W56", "W57",
    "W58", "W59", "W60", "W61", "W62", "W63", "W64", "W65", "W66", "W67", "W68", "W69", "W70", "W73",
    "W74", "W75", "W76", "W77", "W78", "W79", "W80", "W81", "W83", "W84", "W85", "W86", "W87", "W88",
    "W89", "W90", "W91", "W92", "W93", "W94", "W99", "X00", "X01", "X02", "X03", "X04", "X05", "X06",
    "X08", "X09", "X10", "X11", "X12", "X13", "X14", "X15", "X16", "X17", "X18", "X19", "X20", "X21",
    "X22", "X23", "X24", "X25", "X26", "X27", "X30", "X31", "X32", "X33", "X34", "X35", "X36", "X37",
    "X38", "X39", "X40", "X41", "X42", "X43", "X44", "X45", "X46", "X47", "X48", "X49", "X50", "X51",
    "X52", "X53", "X54", "X57", "X58", "X59", "X60", "X61", "X62", "X63", "X64", "X65", "X66", "X67",
    "X68", "X69", "X70", "X71", "X72", "X73", "X74", "X75", "X76", "X77", "X78", "X79", "X80", "X81",
    "X82", "X83", "X84", "X85", "X86", "X87", "X88", "X89", "X90", "X91", "X92", "X93", "X94", "X95",
    "X96", "X97", "X98", "X99", "Y85", "Y86", "Y87", "Y88", "Y89", "Y90", "Y91"
]

secondary_factors = pd.read_csv("/Users/amymann/Documents/Data Quality Project/data_quality/all_secondary_by_ucod_2021.csv")
filtered = secondary_factors[secondary_factors['ucod'].str[:3].isin(accident_list)]
filtered.to_csv("/Users/amymann/Documents/Data Quality Project/data_quality/all_secondary_accident_codes.csv", index=False)
# ───────────────────────────────────────────────────────────────────
#  0 · Client
# ───────────────────────────────────────────────────────────────────
client = openai.OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

system_msg = {
    "role": "system",
    "content": (
        "You are an expert ICD-10 mortality coder. You may use web searches "
        "to find up-to-date coding guidance. Follow ICD-10 conventions and "
        "output exactly in the format requested—no extra commentary."
    ),
}

# ───────────────────────────────────────────────────────────────────
#  Step 1 · CSV of accident UCODs + useful detail codes
# ───────────────────────────────────────────────────────────────────
prompt_step1 = f"""{filtered} 
Step 1:
Attached is a list of accident codes along with all the possible contributing causes that could be listed for any given accident code. Note that the file 
lacks periods but they should be there before the last character. 
Commonly (all though not always), physicians are required to list a contributing cause with many accident codes to give context about the accident.
For example, in an overdose, T codes are used to specify which drug was overdosed on. My goal is to check the percentage of death certificates
that correctly add this context when it is required. Please create a table with the accident UCODs and possible context details. 
Each row should follow this logic: IF NONE OF THE LISTED DETAILS ARE INCLUDED AS A CONTRIBUTING CAUSE CODE THEN THE DEATH WASN'T PROPERLY SPECIFIED.
I.e., one of these details ought to appear. If the accident is already specific enough then list no contributing cause codes.

For example, another code that requires context is W01. An example of a correct coding is:
Immediate cause	e.g. Subdural haemorrhage	S06.5* – nature of injury
b. Intermediate injury	e.g. Fracture, right parietal bone	S02.0* – nature of injury
c. Triggering event	e.g. Fall on same level, slipped on ice while walking	W01 – external-cause code

Make sure to include only specific contributing cause codes, don't include generic, garbage ones like T40.9.
It should put the injury as the ucod and then the possible specific injuries as the contributing cause detail codes.
This is also the way overdoses work. 

Also be sure that you include valid ICD-10 codes that are specific with periods when needed, e.g., T40.4. 

Make sure you go through every accident code in the above list.

Return *only* a CSV table with header:
ucod,detail1,detail2,detail3, detail4, ...
(one row per UCOD; leave empty cells if fewer detail codes exist) No words!!! only ICD-10 codes.
"""

resp1 = client.chat.completions.create(
    model="gpt-4o",
    temperature=0,
    messages=[system_msg, {"role": "user", "content": prompt_step1.strip()}],
)

step1_csv = resp1.choices[0].message.content.strip()

# ----> NEW: save Step 1 output to a .csv file
save_path = (
    "/Users/amymann/Documents/Data Quality Project/data_quality/"
    "accident_detail_codes.csv"
)
with open(save_path, "w", encoding="utf-8") as f:
    f.write(step1_csv + "\n")   # final newline is conventional

print(f"Step 1 CSV saved to: {save_path}")
# ───────────────────────────────────────────────────────────────────
#  Step 2 · R vector of “light-garbage” contributor codes
# ───────────────────────────────────────────────────────────────────
prompt_step2 = """
Step 2:
Create a comprehensive list of ICD-10-CM contributing-cause codes to flag as
“light garbage” because more specific codes exist (many end in .8 or .9, but
include other generic buckets such as T40.2, T43.6, X59, I63.9). Most unspecified codes like .9 should be included as light garbage - make it as thorough as possible.

Return the list on a single line in R syntax, e.g.:

c('T40.9','I63.9','T50.90', ...)

No line breaks, no explanations. Please include ALL such contributing cause codes.
"""

resp2 = client.chat.completions.create(
    model="gpt-4.1",
    temperature=0,
    messages=[system_msg, {"role": "user", "content": prompt_step2.strip()}],
)
# Extract codes from Step 2 R vector
r_vector_raw = resp2.choices[0].message.content.strip()

# Strip c(...) and split by comma
import re
codes = re.findall(r"'([A-Z0-9\.]+)'", r_vector_raw)

# Save to CSV
light_csv_path = (
    "/Users/amymann/Documents/Data Quality Project/data_quality/"
    "light_garbage_codes.csv"
)

import csv
with open(light_csv_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["icd10"])  # header
    for code in codes:
        writer.writerow([code])

print(f"✅ Step 2 CSV saved to: {light_csv_path}")
