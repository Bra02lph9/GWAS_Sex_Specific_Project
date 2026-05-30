from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn2

from src.utils.cli import parse_phenotype
from src.utils.config import get_paths, create_output_dirs


REQUIRED_SUMMARY_COLUMNS = ["group", "count"]

LOCUS_GROUPS = {
    "male_total_loci": "Male total",
    "female_total_loci": "Female total",
    "shared_loci": "Shared",
    "male_only_loci": "Male-only",
    "female_only_loci": "Female-only",
}

BARPLOT_ORDER = [
    "Male total",
    "Female total",
    "Shared",
    "Male-only",
    "Female-only",
]


def save_figure(output_file: Path, dpi: int = 300) -> None:
    """
    Save the current Matplotlib figure as PNG and PDF.
    """

    output_file.parent.mkdir(parents=True, exist_ok=True)

    plt.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.savefig(output_file.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()

    print(f"Saved: {output_file}")
    print(f"Saved: {output_file.with_suffix('.pdf')}")


def load_locus_summary(file_path: Path) -> pd.DataFrame:
    """
    Load and validate the locus summary table.
    """

    if not file_path.exists():
        raise FileNotFoundError(f"Missing locus summary file: {file_path}")

    df = pd.read_csv(file_path, sep="\t")

    missing_columns = [
        col for col in REQUIRED_SUMMARY_COLUMNS if col not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Invalid summary file {file_path.name}. Missing columns: {missing_columns}"
        )

    df = df[REQUIRED_SUMMARY_COLUMNS].copy()
    df["group"] = df["group"].astype(str).str.strip()
    df["count"] = pd.to_numeric(df["count"], errors="coerce")

    df = df.dropna(subset=["group", "count"])
    df = df[df["group"] != ""].copy()
    df["count"] = df["count"].astype(int)

    return df


def get_count(summary_df: pd.DataFrame, group_name: str) -> int:
    """
    Return the count associated with one summary group.
    """

    match = summary_df.loc[summary_df["group"] == group_name, "count"]

    if match.empty:
        raise ValueError(f"Missing group in locus summary: {group_name}")

    return int(match.iloc[0])


def plot_locus_venn(
    male_only: int,
    female_only: int,
    shared: int,
    phenotype: str,
    output_file: Path,
) -> None:
    """
    Plot a Venn diagram showing male-only, female-only, and shared loci.
    """

    plt.figure(figsize=(7, 6))

    venn2(
        subsets=(male_only, female_only, shared),
        set_labels=("Male loci", "Female loci"),
    )

    plt.title(f"{phenotype} Sex-Specific Locus Overlap")
    plt.tight_layout()

    save_figure(output_file)


def prepare_summary_for_plot(summary_df: pd.DataFrame) -> pd.DataFrame:
    """
    Prepare locus summary table for barplot visualization.
    """

    plot_df = summary_df.copy()
    plot_df["label"] = plot_df["group"].map(LOCUS_GROUPS)
    plot_df = plot_df.dropna(subset=["label"]).copy()

    plot_df["label"] = pd.Categorical(
        plot_df["label"],
        categories=BARPLOT_ORDER,
        ordered=True,
    )

    return plot_df.sort_values("label")


def plot_locus_barplot(
    summary_df: pd.DataFrame,
    phenotype: str,
    output_file: Path,
) -> None:
    """
    Plot a barplot summarizing total, shared, and sex-specific loci.
    """

    plot_df = prepare_summary_for_plot(summary_df)

    plt.figure(figsize=(9, 5))
    bars = plt.bar(plot_df["label"], plot_df["count"])

    plt.ylabel("Number of loci")
    plt.title(f"{phenotype} Locus-Level Comparison Summary")
    plt.xticks(rotation=25, ha="right")

    for bar, count in zip(bars, plot_df["count"]):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{int(count):,}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    plt.tight_layout()

    save_figure(output_file)


def main() -> None:
    phenotype = parse_phenotype()
    paths = get_paths(phenotype)

    create_output_dirs(paths)

    summary_file = paths["tables_dir"] / f"{phenotype}_locus_summary.tsv"
    summary_df = load_locus_summary(summary_file)

    male_only = get_count(summary_df, "male_only_loci")
    female_only = get_count(summary_df, "female_only_loci")
    shared = get_count(summary_df, "shared_loci")

    venn_output = paths["figures_dir"] / f"{phenotype}_locus_overlap_venn.png"
    bar_output = paths["figures_dir"] / f"{phenotype}_locus_summary_barplot.png"

    plot_locus_venn(
        male_only=male_only,
        female_only=female_only,
        shared=shared,
        phenotype=phenotype,
        output_file=venn_output,
    )

    plot_locus_barplot(
        summary_df=summary_df,
        phenotype=phenotype,
        output_file=bar_output,
    )

    print("-" * 80)
    print(f"Locus-level visualization completed successfully for phenotype: {phenotype}")


if __name__ == "__main__":
    main()
