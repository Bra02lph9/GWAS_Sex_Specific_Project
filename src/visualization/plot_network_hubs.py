from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib_venn import venn3

from src.utils.cli import parse_phenotype
from src.utils.config import get_paths, create_output_dirs


TOP_N = 15
NETWORK_NAMES = ["female", "male", "shared"]

METRICS_TO_PLOT = [
    "weighted_degree",
    "betweenness_centrality",
]

METRIC_LABELS = {
    "degree_count": "Degree Count",
    "weighted_degree": "Weighted Degree",
    "betweenness_centrality": "Betweenness Centrality",
}

CATEGORY_COLORS = {
    "Lipid metabolism": "#1f77b4",
    "Inflammation / immune": "#d62728",
    "Mitochondrial / metabolic": "#2ca02c",
    "ECM / vascular remodeling": "#ff7f0e",
    "Intracellular signaling": "#9467bd",
    "Regulatory / epigenetic": "#8c564b",
    "RNA / ribosomal / translation": "#7f7f7f",
    "Other": "#bdbdbd",
}


GENE_CATEGORIES = {
    "shared": {
        "APOE": "Lipid metabolism",
        "APOB": "Lipid metabolism",
        "LPL": "Lipid metabolism",
        "CETP": "Lipid metabolism",
        "LPA": "Lipid metabolism",
        "APOA5": "Lipid metabolism",

        "FN1": "ECM / vascular remodeling",
        "SMAD3": "ECM / vascular remodeling",

        "RELA": "Inflammation / immune",
        "RHOA": "Inflammation / immune",

        "AKT1": "Intracellular signaling",
        "CDKN2A": "Intracellular signaling",

        "KAT5": "Regulatory / epigenetic",
        "SMARCA4": "Regulatory / epigenetic",
        "UBC": "Regulatory / epigenetic",
    },

    "male": {
        "STAT3": "Inflammation / immune",
        "IL1B": "Inflammation / immune",
        "CD4": "Inflammation / immune",
        "FLT3LG": "Inflammation / immune",
        "IFI35": "Inflammation / immune",

        "PLCG1": "Intracellular signaling",
        "NRAS": "Intracellular signaling",
        "GNB1": "Intracellular signaling",
        "GNAI3": "Intracellular signaling",

        "RPS3": "RNA / ribosomal / translation",
        "RPS28": "RNA / ribosomal / translation",
        "RPS20": "RNA / ribosomal / translation",
        "RPLP2": "RNA / ribosomal / translation",
        "RPL6": "RNA / ribosomal / translation",
        "RPL5": "RNA / ribosomal / translation",
        "RPL37A": "RNA / ribosomal / translation",
        "RPL37": "RNA / ribosomal / translation",
        "RPL13A": "RNA / ribosomal / translation",
        "RPL10L": "RNA / ribosomal / translation",
        "MRPL24": "RNA / ribosomal / translation",
        "EEF2": "RNA / ribosomal / translation",
        "POLR2B": "RNA / ribosomal / translation",
        "H3-3B": "RNA / ribosomal / translation",

        "BRCA1": "Regulatory / epigenetic",
        "DOCK4": "Other",
        "CCDC105": "Other",
    },

    "female": {
        "PPARGC1A": "Mitochondrial / metabolic",
        "CS": "Mitochondrial / metabolic",
        "NDUFS3": "Mitochondrial / metabolic",
        "MRPL17": "Mitochondrial / metabolic",
        "HSPA9": "Mitochondrial / metabolic",
        "LEP": "Mitochondrial / metabolic",
        "POMC": "Mitochondrial / metabolic",

        "JAK2": "Inflammation / immune",
        "STAT2": "Inflammation / immune",

        "ERBB2": "Intracellular signaling",
        "GNAI2": "Intracellular signaling",
        "GNAQ": "Intracellular signaling",
        "KDM1A": "Intracellular signaling",
        "CCNB2": "Intracellular signaling",

        "MTREX": "RNA / ribosomal / translation",
        "UTP4": "RNA / ribosomal / translation",
        "SNRPE": "RNA / ribosomal / translation",
    },
}


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


def get_gene_category(gene: str, network_name: str) -> str:
    return GENE_CATEGORIES.get(network_name, {}).get(gene, "Other")


def add_categories(df: pd.DataFrame, network_name: str) -> pd.DataFrame:
    df = df.copy()
    df["category"] = df["gene"].apply(
        lambda gene: get_gene_category(gene, network_name)
    )
    df["color"] = df["category"].map(CATEGORY_COLORS).fillna(CATEGORY_COLORS["Other"])
    return df


def save_figure(output_file: Path) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)

    allowed_formats = [".png", ".jpg", ".jpeg"]

    if output_file.suffix.lower() not in allowed_formats:
        output_file = output_file.with_suffix(".png")

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(f"Saved: {output_file}")


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

    top_df = add_categories(top_df, network_name)
    top_df = top_df.sort_values(by=metric, ascending=True)

    metric_label = METRIC_LABELS.get(metric, metric.replace("_", " ").title())
    network_label = network_name.capitalize()

    fig, ax = plt.subplots(figsize=(10.5, 6.2))

    ax.barh(
        top_df["gene"],
        top_df[metric],
        color=top_df["color"],
        edgecolor="black",
        linewidth=0.4,
    )

    ax.set_xlabel(metric_label, fontsize=11)
    ax.set_ylabel("Gene", fontsize=11)
    ax.set_title(
        f"{phenotype} {network_label} Hub Genes Ranked by {metric_label}",
        fontsize=13,
        fontweight="bold",
    )

    ax.tick_params(axis="both", labelsize=10)
    ax.grid(axis="x", linestyle="--", alpha=0.25)
    ax.set_axisbelow(True)

    categories_present = list(dict.fromkeys(top_df["category"].tolist()))

    legend_handles = [
        Patch(
            facecolor=CATEGORY_COLORS[category],
            edgecolor="black",
            label=category,
        )
        for category in categories_present
    ]

    ax.legend(
        handles=legend_handles,
        title="Functional category",
        fontsize=8,
        title_fontsize=9,
        loc="lower right",
        frameon=True,
    )

    plt.tight_layout()
    save_figure(output_file)


def plot_hub_venn(
    female_df: pd.DataFrame,
    male_df: pd.DataFrame,
    shared_df: pd.DataFrame,
    phenotype: str,
    output_file: Path,
    top_n: int = TOP_N,
) -> None:
    female_genes = set(
        female_df.sort_values(by="weighted_degree", ascending=False)
        .head(top_n)["gene"]
    )

    male_genes = set(
        male_df.sort_values(by="weighted_degree", ascending=False)
        .head(top_n)["gene"]
    )

    shared_genes = set(
        shared_df.sort_values(by="weighted_degree", ascending=False)
        .head(top_n)["gene"]
    )

    plt.figure(figsize=(7, 7))

    venn3(
        [female_genes, male_genes, shared_genes],
        set_labels=("Female hubs", "Male hubs", "Shared hubs"),
    )

    plt.title(
        f"{phenotype} Top {top_n} Hub Gene Overlap",
        fontsize=13,
        fontweight="bold",
    )

    plt.tight_layout()
    save_figure(output_file)


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
        for metric in METRICS_TO_PLOT:
            output_file = (
                paths["figures_dir"]
                / f"{phenotype}_{network_name}_top_hubs_{metric}.png"
            )

            plot_top_hubs(
                df=df,
                phenotype=phenotype,
                network_name=network_name,
                metric=metric,
                output_file=output_file,
            )

    plot_hub_venn(
        female_df=hub_data["female"],
        male_df=hub_data["male"],
        shared_df=hub_data["shared"],
        phenotype=phenotype,
        output_file=paths["figures_dir"] / f"{phenotype}_hub_genes_overlap_venn.png",
    )

    print("-" * 60)
    print(f"Hub gene visualization completed successfully for: {phenotype}")


if __name__ == "__main__":
    main()
