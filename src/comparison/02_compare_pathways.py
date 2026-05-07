import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

# ── Paths ─────────────────────────────────────────────────────────────────────
INPUT_DIR  = "../GWAS_Sex_Specific_Project/3_tools_results/gprofiler"
OUTPUT_DIR = "../results/02_pathways"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Load gProfiler files ───────────────────────────────────────────────────────
female_df = pd.read_csv(f"{INPUT_DIR}/female_pathways.tsv", sep=",", quotechar='"')
male_df   = pd.read_csv(f"{INPUT_DIR}/male_pathways.tsv",   sep=",", quotechar='"')
shared_df = pd.read_csv(f"{INPUT_DIR}/shared_pathways.tsv", sep=",", quotechar='"')

# Strip any leftover quotes from column names
female_df.columns = female_df.columns.str.strip('"')
male_df.columns   = male_df.columns.str.strip('"')
shared_df.columns = shared_df.columns.str.strip('"')

print(f"Female pathways loaded:  {female_df.shape[0]} rows")
print(f"Male pathways loaded:    {male_df.shape[0]} rows")
print(f"Shared pathways loaded:  {shared_df.shape[0]} rows")
print(f"\nColumns: {female_df.columns.tolist()}")

# ── Filter significant pathways (adjusted p-value < 0.05) ────────────────────
female_sig = female_df[female_df["adjusted_p_value"] < 0.05].copy()
male_sig   = male_df[male_df["adjusted_p_value"]     < 0.05].copy()
shared_sig = shared_df[shared_df["adjusted_p_value"] < 0.05].copy()

print(f"\nSignificant (adj.p < 0.05):")
print(f"  Female: {len(female_sig)}")
print(f"  Male:   {len(male_sig)}")
print(f"  Shared: {len(shared_sig)}")

# ── Keep only 'highlighted' pathways (top enriched per gProfiler) ─────────────
female_hl = female_sig[female_sig["highlighted"] == True].copy()
male_hl   = male_sig[male_sig["highlighted"]     == True].copy()

print(f"\nHighlighted significant pathways:")
print(f"  Female: {len(female_hl)}")
print(f"  Male:   {len(male_hl)}")

# ── Identify sex-specific vs shared pathway term IDs ─────────────────────────
female_terms = set(female_sig["term_id"])
male_terms   = set(male_sig["term_id"])
shared_terms = set(shared_sig["term_id"])

female_specific_ids = female_terms - male_terms - shared_terms
male_specific_ids   = male_terms   - female_terms - shared_terms
truly_shared_ids    = female_terms & male_terms

print(f"\nPathway overlap:")
print(f"  Female-specific:  {len(female_specific_ids)}")
print(f"  Male-specific:    {len(male_specific_ids)}")
print(f"  Shared (F & M):   {len(truly_shared_ids)}")

# ── Extract sex-specific dataframes ───────────────────────────────────────────
female_specific = female_sig[female_sig["term_id"].isin(female_specific_ids)].copy()
male_specific   = male_sig[male_sig["term_id"].isin(male_specific_ids)].copy()

# ── Save outputs ──────────────────────────────────────────────────────────────
female_specific.to_csv(f"{OUTPUT_DIR}/female_specific_pathways.tsv", sep="\t", index=False)
male_specific.to_csv(f"{OUTPUT_DIR}/male_specific_pathways.tsv",     sep="\t", index=False)

print(f"\nSaved:")
print(f"  {OUTPUT_DIR}/female_specific_pathways.tsv")
print(f"  {OUTPUT_DIR}/male_specific_pathways.tsv")

# ── Top pathways for plotting ─────────────────────────────────────────────────
def get_top_pathways(df, n=15):
    df = df.copy()
    df["-log10_p"] = -np.log10(df["adjusted_p_value"])
    df = df.sort_values("-log10_p", ascending=False)
    # Truncate long names
    df["term_name_short"] = df["term_name"].str[:50]
    return df.head(n)

top_female = get_top_pathways(female_specific, n=15)
top_male   = get_top_pathways(male_specific,   n=15)

# ── Plot: Side-by-side barplots ───────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(18, 8))
fig.suptitle("Sex-Specific Pathway Enrichment in CAD (gProfiler)", fontsize=14, fontweight="bold", y=1.02)

colors = {"female": "#E8537A", "male": "#4A90D9"}

for ax, df, sex, color in [
    (axes[0], top_female, "Female-Specific", colors["female"]),
    (axes[1], top_male,   "Male-Specific",   colors["male"])
]:
    if df.empty:
        ax.text(0.5, 0.5, f"No specific pathways\nfor {sex}", ha="center", va="center", fontsize=12)
        ax.set_title(sex)
        continue

    bars = ax.barh(
        df["term_name_short"][::-1],
        df["-log10_p"][::-1],
        color=color, alpha=0.85, edgecolor="white", linewidth=0.5
    )
    ax.axvline(x=-np.log10(0.05), color="black", linestyle="--", linewidth=1, label="p=0.05")
    ax.set_xlabel("-log10(adjusted p-value)", fontsize=11)
    ax.set_title(f"{sex} Pathways\n(n={len(df)})", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", labelsize=9)

plt.tight_layout()
plot_path = f"{OUTPUT_DIR}/pathway_comparison_plot.png"
plt.savefig(plot_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"  {plot_path}")

# ── Plot: Source breakdown (GO:BP, GO:MF, KEGG, etc.) ────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle("Pathway Source Distribution", fontsize=13, fontweight="bold")

for ax, df, sex, color in [
    (axes[0], female_specific, "Female-Specific", colors["female"]),
    (axes[1], male_specific,   "Male-Specific",   colors["male"])
]:
    if df.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.set_title(sex)
        continue
    source_counts = df["source"].value_counts()
    ax.bar(source_counts.index, source_counts.values, color=color, alpha=0.85, edgecolor="white")
    ax.set_title(sex, fontweight="bold")
    ax.set_ylabel("Number of pathways")
    ax.tick_params(axis="x", rotation=30)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

plt.tight_layout()
source_plot_path = f"{OUTPUT_DIR}/pathway_source_distribution.png"
plt.savefig(source_plot_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"  {source_plot_path}")

# ── Summary report ────────────────────────────────────────────────────────────
with open(f"{OUTPUT_DIR}/pathway_summary.txt", "w", encoding="utf-8") as f:
    f.write("=== gProfiler Pathway Comparison Summary ===\n\n")
    f.write(f"Total significant pathways (adj.p < 0.05):\n")
    f.write(f"  Female: {len(female_sig)}\n")
    f.write(f"  Male:   {len(male_sig)}\n")
    f.write(f"  Shared: {len(shared_sig)}\n\n")
    f.write(f"Sex-specific pathways:\n")
    f.write(f"  Female-only: {len(female_specific)}\n")
    f.write(f"  Male-only:   {len(male_specific)}\n")
    f.write(f"  Shared (F & M overlap): {len(truly_shared_ids)}\n\n")

    f.write("--- Top 10 Female-Specific Pathways ---\n")
    for _, row in top_female.head(10).iterrows():
        f.write(f"  [{row['source']}] {row['term_name']} (adj.p={row['adjusted_p_value']:.2e})\n")

    f.write("\n--- Top 10 Male-Specific Pathways ---\n")
    for _, row in top_male.head(10).iterrows():
        f.write(f"  [{row['source']}] {row['term_name']} (adj.p={row['adjusted_p_value']:.2e})\n")

print(f"  {OUTPUT_DIR}/pathway_summary.txt")
print("\n✅ Pathway comparison complete.")