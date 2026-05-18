from pathlib import Path
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[2]

GPROFILER_DIR = PROJECT_ROOT / "3_tools_results" / "gprofiler"
FIGURES_DIR = PROJECT_ROOT / "4_results" / "figures"

TOP_N = 15


ENRICHMENT_FILES = {
    "male": GPROFILER_DIR / "male_pathways.tsv",
    "female": GPROFILER_DIR / "female_pathways.tsv",
    "shared": GPROFILER_DIR / "shared_pathways.tsv",
}


REQUIRED_COLUMNS = [
    "term_name",
    "adjusted_p_value",
    "intersection_size",
]


def load_enrichment_results(file_path: Path) -> pd.DataFrame:
    if not file_path.exists():
        raise FileNotFoundError(f"Missing enrichment file: {file_path}")

    df = pd.read_csv(file_path, sep=",", low_memory=False)

    missing_columns = [col for col in REQUIRED_COLUMNS if col not in df.columns]

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
        lambda x: "\n".join(textwrap.wrap(x, width=45))
    )

    df = df.sort_values(
        by="minus_log10_adjusted_p",
        ascending=True
    )

    return df


def plot_pathway_dotplot(
    df: pd.DataFrame,
    title: str,
    output_file: Path,
) -> None:
    if df.empty:
        print(f"Skipped empty dotplot: {title}")
        return

    plt.figure(figsize=(10, 7))

    scatter = plt.scatter(
        df["minus_log10_adjusted_p"],
        df["term_label"],
        s=df["intersection_size"] * 35,
        c=df["adjusted_p_value"],
        alpha=0.8,
    )

    plt.xlabel("-log10(adjusted p-value)")
    plt.ylabel("Enriched pathway / term")
    plt.title(title)

    cbar = plt.colorbar(scatter)
    cbar.set_label("Adjusted p-value")

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Saved pathway dotplot: {output_file}")


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    for group_name, file_path in ENRICHMENT_FILES.items():
        print(f"Processing {group_name} enrichment results...")

        df = load_enrichment_results(file_path)
        top_terms = prepare_top_terms(df)

        output_file = FIGURES_DIR / f"{group_name}_pathway_enrichment_dotplot.png"

        plot_pathway_dotplot(
            top_terms,
            f"Top Enriched Pathways - {group_name.capitalize()} CAD Genes",
            output_file
        )

    print("-" * 60)
    print("Pathway enrichment dotplots completed successfully.")


if __name__ == "__main__":
    main()
