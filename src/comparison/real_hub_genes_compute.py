from pathlib import Path
import pandas as pd
import networkx as nx

from src.utils.cli import parse_phenotype
from src.utils.config import get_paths, create_output_dirs


NETWORK_NAMES = ["female", "male", "shared"]


def load_network(file_path: Path) -> pd.DataFrame:
    if not file_path.exists():
        raise FileNotFoundError(f"Missing network file: {file_path}")

    df = pd.read_csv(file_path, sep="\t", low_memory=False)

    required_columns = ["source", "target"]
    missing = [col for col in required_columns if col not in df.columns]

    if missing:
        raise ValueError(f"Missing columns in {file_path.name}: {missing}")

    df["source"] = df["source"].astype(str).str.strip()
    df["target"] = df["target"].astype(str).str.strip()

    df = df[
        (df["source"] != "")
        & (df["target"] != "")
        & (df["source"].str.lower() != "nan")
        & (df["target"].str.lower() != "nan")
        & (df["source"] != df["target"])
    ].copy()

    if "combined_score" in df.columns:
        df["combined_score"] = pd.to_numeric(df["combined_score"], errors="coerce")
        df["combined_score"] = df["combined_score"].fillna(1.0)
    else:
        df["combined_score"] = 1.0

    df = df.drop_duplicates(subset=["source", "target"]).copy()

    return df


def calculate_hub_genes(df: pd.DataFrame) -> pd.DataFrame:
    graph = nx.from_pandas_edgelist(
        df,
        source="source",
        target="target",
        edge_attr="combined_score",
        create_using=nx.Graph(),
    )

    degree_count = dict(graph.degree())
    degree_centrality = nx.degree_centrality(graph)
    betweenness = nx.betweenness_centrality(graph)
    weighted_degree = dict(graph.degree(weight="combined_score"))

    hub_df = pd.DataFrame({
        "gene": list(graph.nodes()),
        "degree_count": [degree_count[gene] for gene in graph.nodes()],
        "degree_centrality": [degree_centrality[gene] for gene in graph.nodes()],
        "weighted_degree": [weighted_degree[gene] for gene in graph.nodes()],
        "betweenness_centrality": [betweenness[gene] for gene in graph.nodes()],
    })

    return hub_df.sort_values(
        by=["degree_count", "weighted_degree", "betweenness_centrality"],
        ascending=False,
    )


def main() -> None:
    phenotype = parse_phenotype()
    paths = get_paths(phenotype)
    create_output_dirs(paths)

    for network_name in NETWORK_NAMES:
        input_file = paths["cytoscape_dir"] / f"{network_name}_network_clean.tsv"

        print(f"Processing {phenotype} {network_name} network...")

        df = load_network(input_file)
        hub_df = calculate_hub_genes(df)

        output_file = paths["tables_dir"] / f"{phenotype}_{network_name}_hub_genes.tsv"

        hub_df.to_csv(output_file, sep="\t", index=False)

        print(f"Nodes: {len(hub_df):,}")
        print(f"Edges: {len(df):,}")
        print(f"Saved to: {output_file}")
        print("-" * 60)


if __name__ == "__main__":
    main()
