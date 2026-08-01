import csv

INPUT_FILE = "variant_summary.txt/variant_summary.txt"
OUTPUT_FILE = "alzheimers_clinvar.csv"

# Established Alzheimer's disease-associated genes
ALZ_GENES = {
    "APOE",
    "APP",
    "PSEN1",
    "PSEN2",
    "TREM2",
    "SORL1",
    "ABCA7",
    "BIN1",
    "CLU",
    "CR1",
    "PICALM",
    "CD33",
    "CD2AP",
    "MS4A6A",
    "EPHA1"
}

FIELDNAMES = [
    "GeneSymbol",
    "Chromosome",
    "Start",
    "Type",
    "Name",
    "ClinicalSignificance",
    "PhenotypeList",
    "ReviewStatus",
    "RS# (dbSNP)",
    "AlleleID",
    "Assembly"
]

rows_scanned = 0
rows_written = 0

seen = set()

with open(INPUT_FILE, "r", encoding="utf-8", newline="") as infile, \
     open(OUTPUT_FILE, "w", encoding="utf-8", newline="") as outfile:

    reader = csv.DictReader(infile, delimiter="\t")
    writer = csv.DictWriter(outfile, fieldnames=FIELDNAMES)

    writer.writeheader()

    for row in reader:

        rows_scanned += 1

        if rows_scanned % 100000 == 0:
            print(f"Scanned {rows_scanned:,} rows...")

        gene = row.get("GeneSymbol", "").strip().upper()

        if gene not in ALZ_GENES:
            continue

        if row.get("Assembly") != "GRCh38":
            continue

        # Skip exact duplicate variants
        key = (
            gene,
            row.get("Start"),
            row.get("Name")
        )

        if key in seen:
            continue

        seen.add(key)

        writer.writerow({
            field: row.get(field, "")
            for field in FIELDNAMES
        })

        rows_written += 1

print("\nExtraction complete!")
print(f"Rows scanned : {rows_scanned:,}")
print(f"Rows written : {rows_written:,}")
print(f"Saved to      : {OUTPUT_FILE}")

import pandas as pd

df = pd.read_csv("alzheimers_clinvar.csv")

print(df.shape)
print(df["GeneSymbol"].value_counts())
print(df["ClinicalSignificance"].value_counts().head(20))