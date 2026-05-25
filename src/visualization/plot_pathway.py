from pathlib import Path
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.utils.cli import parse_phenotype
from src.utils.config import get_paths, create_output_dirs


TOP_N = 15

REQUIRED_COLUMNS = [
    "term_name",
    "adjusted_p_value",
    "intersection_size",
]

GROUPS = ["male", "female", "shared"]


def load_enrichment_results(file_path: Path) -> pd.DataFrame:
    if not file_path.exists():
        raise FileNotFoundError(f"Missing enrichment file: {file_path}")

    df = pd.read_csv(file_path, sep="\t", low_memory=False)

    missing_columns = [
        col for col in REQUIRED_COLUMNS
        if col not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Missing required columns in {file_path.name}: {missing_columns}"
        )

    df = df.copy()

    df["term_name"] = df["term_name"].astype(str).str.strip()

    df["adjusted_p_value"] = pd.to_numeric(
        df["adjusted_p_value"],
        errors="coerce"
    )

    df["intersection_size"] = pd.to_numeric(
        df["intersection_size"],
        errors="coerce"
    )

    df = df.dropna(
        subset=["term_name", "adjusted_p_value", "intersection_size"]
    )

    df = df[
        (df["term_name"] != "")
        & (df["term_name"].str.lower() != "nan")
        & (df["adjusted_p_value"] > 0)
        & (df["adjusted_p_value"] <= 1)
        & (df["intersection_size"] > 0)
    ].copy()

    if df.empty:
        raise ValueError(f"No valid enrichment rows found in {file_path.name}")

    return df


def prepare_top_terms(df: pd.DataFrame, top_n: int = TOP_N) -> pd.DataFrame:
    df = df.sort_values(
        by="adjusted_p_value",
        ascending=True
    ).head(top_n).copy()

    df["minus_log10_adjusted_p"] = -np.log10(df["adjusted_p_value"])

    df["term_label"] = df["term_name"].apply(
        lambda x: "\n".join(textwrap.wrap(x, width=40))
    )

    df = df.sort_values(
        by="minus_log10_adjusted_p",
        ascending=True
    )

    return df


def normalize_bubble_sizes(values, min_size=60, max_size=500):
    values = np.array(values)

    if values.max() == values.min():
        return np.full_like(values, (min_size + max_size) / 2)

    normalized = (
        (values - values.min())
        / (values.max() - values.min())
    )

    return min_size + normalized * (max_size - min_size)


def plot_pathway_dotplot(
    df: pd.DataFrame,
    title: str,
    output_file: Path,
) -> None:
    if df.empty:
        print(f"Skipped empty dotplot: {title}")
        return

    bubble_sizes = normalize_bubble_sizes(
        df["intersection_size"]
    )

    plt.figure(figsize=(11, 7))

    scatter = plt.scatter(
        df["minus_log10_adjusted_p"],
        df["term_label"],
        s=bubble_sizes,
        c=df["adjusted_p_value"],
        alpha=0.75,
        edgecolors="black",
        linewidths=0.5,
    )

    plt.xlabel("-log10(adjusted p-value)")
    plt.ylabel("Enriched pathway / term")
    plt.title(title)

    cbar = plt.colorbar(scatter)
    cbar.set_label("Adjusted p-value")

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved pathway dotplot: {output_file}")


def main() -> None:
    phenotype = parse_phenotype()
    paths = get_paths(phenotype)

    create_output_dirs(paths)

    for group_name in GROUPS:
        file_path = (
            paths["gprofiler_dir"]
            / f"{group_name}_pathways.tsv"
        )

        print(f"Processing {group_name} enrichment results...")

        df = load_enrichment_results(file_path)
        top_terms = prepare_top_terms(df)

        output_file = (
            paths["figures_dir"]
            / f"{phenotype}_{group_name}_pathway_enrichment_dotplot.png"
        )

        plot_pathway_dotplot(
            top_terms,
            f"Top Enriched Pathways - {group_name.capitalize()} {phenotype} Genes",
            output_file,
        )

    print("-" * 60)
    print(f"Pathway enrichment dotplots completed successfully for: {phenotype}")


if __name__ == "__main__":
    main()
