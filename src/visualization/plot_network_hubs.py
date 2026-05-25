from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

from src.utils.cli import parse_phenotype
from src.utils.config import get_paths, create_output_dirs


TOP_N = 15
NETWORK_NAMES = ["female", "male", "shared"]


def load_hub_genes(file_path: Path) -> pd.DataFrame:
    if not file_path.exists():
        raise FileNotFoundError(f"Missing hub gene file: {file_path}")

    df = pd.read_csv(file_path, sep="\t", low_memory=False)

    required_columns = [
        "gene",
        "degree_count",
        "weighted_degree",
        "betweenness_centrality",
    ]

    missing_columns = [col for col in required_columns if col not in df.columns]

    if missing_columns:
        raise ValueError(
            f"Missing columns in {file_path.name}: {missing_columns}"
        )

    df["gene"] = df["gene"].astype(str).str.strip()

    for col in ["degree_count", "weighted_degree", "betweenness_centrality"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(
        subset=["gene", "degree_count", "weighted_degree", "betweenness_centrality"]
    )

    df = df[
        (df["gene"] != "")
        & (df["gene"].str.lower() != "nan")
    ].copy()

    return df


def plot_top_hubs(
    df: pd.DataFrame,
    phenotype: str,
    network_name: str,
    metric: str,
    output_file: Path,
    top_n: int = TOP_N,
) -> None:
    top_df = df.sort_values(by=metric, ascending=False).head(top_n).copy()

    if top_df.empty:
        print(f"Skipped empty plot: {network_name} - {metric}")
        return

    top_df = top_df.sort_values(by=metric, ascending=True)

    plt.figure(figsize=(9, 6))
    plt.barh(top_df["gene"], top_df[metric])
    plt.xlabel(metric.replace("_", " ").title())
    plt.ylabel("Gene")
    plt.title(
        f"{phenotype} Top {top_n} {network_name.capitalize()} Hub Genes by {metric}"
    )
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Saved: {output_file}")


def plot_hub_venn(
    female_df: pd.DataFrame,
    male_df: pd.DataFrame,
    shared_df: pd.DataFrame,
    phenotype: str,
    output_file: Path,
    top_n: int = TOP_N,
) -> None:
    female_genes = set(
        female_df.sort_values(by="degree_count", ascending=False)
        .head(top_n)["gene"]
    )

    male_genes = set(
        male_df.sort_values(by="degree_count", ascending=False)
        .head(top_n)["gene"]
    )

    shared_genes = set(
        shared_df.sort_values(by="degree_count", ascending=False)
        .head(top_n)["gene"]
    )

    plt.figure(figsize=(7, 7))
    venn3(
        [female_genes, male_genes, shared_genes],
        set_labels=("Female hubs", "Male hubs", "Shared hubs"),
    )
    plt.title(f"{phenotype} Top {top_n} Hub Gene Overlap")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Saved: {output_file}")


def main() -> None:
    phenotype = parse_phenotype()
    paths = get_paths(phenotype)
    create_output_dirs(paths)

    hub_files = {
        name: paths["tables_dir"] / f"{phenotype}_{name}_hub_genes.tsv"
        for name in NETWORK_NAMES
    }

    hub_data = {
        name: load_hub_genes(path)
        for name, path in hub_files.items()
    }

    for network_name, df in hub_data.items():
        plot_top_hubs(
            df,
            phenotype,
            network_name,
            "degree_count",
            paths["figures_dir"] / f"{phenotype}_{network_name}_top_hubs_degree.png",
        )

        plot_top_hubs(
            df,
            phenotype,
            network_name,
            "weighted_degree",
            paths["figures_dir"] / f"{phenotype}_{network_name}_top_hubs_weighted_degree.png",
        )

        plot_top_hubs(
            df,
            phenotype,
            network_name,
            "betweenness_centrality",
            paths["figures_dir"] / f"{phenotype}_{network_name}_top_hubs_betweenness.png",
        )

    plot_hub_venn(
        hub_data["female"],
        hub_data["male"],
        hub_data["shared"],
        phenotype,
        paths["figures_dir"] / f"{phenotype}_hub_genes_overlap_venn.png",
    )

    print("-" * 60)
    print(f"Hub gene visualization completed successfully for: {phenotype}")


if __name__ == "__main__":
    main()
