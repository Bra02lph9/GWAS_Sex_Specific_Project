from pathlib import Path
import pandas as pd
import networkx as nx


PROJECT_ROOT = Path(__file__).resolve().parents[2]

CYTOSCAPE_DIR = PROJECT_ROOT / "3_tools_results" / "cytoscape"
TABLES_DIR = PROJECT_ROOT / "4_results" / "tables"

TABLES_DIR.mkdir(parents=True, exist_ok=True)


NETWORKS = {
    "female": CYTOSCAPE_DIR / "female_network_clean.tsv",
    "male": CYTOSCAPE_DIR / "male_network_clean.tsv",
    "shared": CYTOSCAPE_DIR / "shared_network_clean.tsv",
}


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
    G = nx.from_pandas_edgelist(
        df,
        source="source",
        target="target",
        edge_attr="combined_score",
        create_using=nx.Graph()
    )

    degree_count = dict(G.degree())
    degree_centrality = nx.degree_centrality(G)
    betweenness = nx.betweenness_centrality(G)

    weighted_degree = dict(G.degree(weight="combined_score"))

    hub_df = pd.DataFrame({
        "gene": list(G.nodes()),
        "degree_count": [degree_count[gene] for gene in G.nodes()],
        "degree_centrality": [degree_centrality[gene] for gene in G.nodes()],
        "weighted_degree": [weighted_degree[gene] for gene in G.nodes()],
        "betweenness_centrality": [betweenness[gene] for gene in G.nodes()],
    })

    hub_df = hub_df.sort_values(
        by=["degree_count", "weighted_degree", "betweenness_centrality"],
        ascending=False
    )

    return hub_df


def main() -> None:
    for network_name, input_file in NETWORKS.items():
        print(f"Processing {network_name} network...")

        df = load_network(input_file)
        hub_df = calculate_hub_genes(df)

        output_file = TABLES_DIR / f"{network_name}_hub_genes.tsv"

        hub_df.to_csv(
            output_file,
            sep="\t",
            index=False
        )

        print(f"Nodes: {len(hub_df):,}")
        print(f"Edges: {len(df):,}")
        print(f"Saved to: {output_file}")
        print("-" * 60)


if __name__ == "__main__":
    main()
