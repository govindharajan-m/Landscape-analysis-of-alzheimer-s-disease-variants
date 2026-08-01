import pandas as pd
import numpy as np
import csv
import os


df = pd.read_csv("alzheimers_clinvar.csv")
print(df.shape)
print(df.info())
print(df.head())
#df = df.drop_duplicates(subset="AlleleID")
df = df.rename(columns={
    "GeneSymbol": "Gene",
    "Start": "Position",
    "RS# (dbSNP)": "dbSNP"
})
df["ClinicalSignificance"] = (
    df["ClinicalSignificance"]
    .str.strip()
    .str.title()
)
df.to_csv("alzheimers_clean.csv", index=False)

gene_counts = (
    df["Gene"]
    .value_counts()
)

print(gene_counts.head(20))

import matplotlib.pyplot as plt

gene_counts.head(15).plot(kind="bar", figsize=(10,5))

plt.title("Top Alzheimer's Genes by Variant Count")
plt.xlabel("Gene")
plt.ylabel("Number of Variants")

plt.tight_layout()
plt.show()

chrom = (
    df["Chromosome"]
    .value_counts()
    .sort_index()
)

chrom.plot(kind="bar")

plt.title("Variant Distribution Across Chromosomes")

plt.show()

types = df["Type"].value_counts()

types.plot(kind="bar")

plt.title("Variant Types")

plt.show()

weights = {
    "Pathogenic": 5,
    "Likely Pathogenic": 4,
    "Uncertain Significance": 2,
    "Likely Benign": 1,
    "Benign": 0
}

df["RiskScore"] = (
    df["ClinicalSignificance"]
    .map(weights)
    .fillna(0)
)
gene_stats = (
    df.groupby("Gene")
      .agg(
          VariantCount=("AlleleID", "count"),
          BurdenScore=("RiskScore", "sum"),
          SeverityIndex=("RiskScore", "mean")
      )
)
gene_stats["NormalizedBurden"] = (
    gene_stats["BurdenScore"] /
    gene_stats["BurdenScore"].max()
) * 100
gene_stats["CompositeScore"] = (
    0.6 * gene_stats["NormalizedBurden"] +
    0.4 * (gene_stats["SeverityIndex"] * 20)
)
gene_stats = gene_stats.sort_values(
    "CompositeScore",
    ascending=False
)

print(gene_stats.head(20))
gene_stats.sort_values(
    "BurdenScore",
    ascending=False
).head(15)["BurdenScore"].plot.bar(figsize=(10,5))

plt.title("Gene Burden Score")
plt.ylabel("Weighted Clinical Significance")
plt.tight_layout()
plt.show()
gene_stats.sort_values(
    "SeverityIndex",
    ascending=False
).head(15)["SeverityIndex"].plot.bar(figsize=(10,5))

plt.title("Gene Severity Index")
plt.ylabel("Mean Variant Severity")
plt.tight_layout()
plt.show()
gene_stats.sort_values(
    "CompositeScore",
    ascending=False
).head(15)["CompositeScore"].plot.bar(figsize=(10,5))

plt.title("Composite Alzheimer's Gene Risk Score")
plt.ylabel("Composite Score (0–100)")
plt.tight_layout()
plt.show()

gene_scores = (
    df.groupby("Gene")["RiskScore"]
      .sum()
      .sort_values(ascending=False)
)

print(gene_scores.head(20))

gene_scores.head(15).plot(
    kind="bar",
    figsize=(10,5),
    color="crimson"
)

plt.title("Alzheimer's Gene Risk Score")

plt.ylabel("Risk Score")

plt.tight_layout()

plt.show()

gene_stats = (
    df.groupby("Gene")
      .agg(
          TotalScore=("RiskScore","sum"),
          Variants=("RiskScore","count")
      )
)

gene_stats["SeverityIndex"] = (
    gene_stats["TotalScore"] /
    gene_stats["Variants"]
)

gene_stats = gene_stats.sort_values(
    "SeverityIndex",
    ascending=False
)

print(gene_stats.head(20))

import seaborn as sns

pivot = (
    df.groupby(
        ["Gene","ClinicalSignificance"]
    )
    .size()
    .unstack(fill_value=0)
)
sns.heatmap(
    pivot,
    cmap="viridis",
    annot=True,
    fmt="d"
)

plt.title("Clinical Significance by Gene")

plt.show()

pivot = pd.crosstab(
    df["Gene"],
    df["ClinicalSignificance"]
)

pivot.plot(
    kind="bar",
    stacked=True,
    figsize=(12,6)
)

plt.title("Clinical Significance Distribution Across Alzheimer's Genes")
plt.xlabel("Gene")
plt.ylabel("Number of Variants")

plt.tight_layout()
plt.show()