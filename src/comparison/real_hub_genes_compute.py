import networkx as nx
import pandas as pd

base_path = "../GWAS_Sex_Specific_Project/3_tools_results/cytoscape"

networks = {
    "female": f"{base_path}/female_network_clean.tsv",
    "male":   f"{base_path}/male_network_clean.tsv",
    "shared": f"{base_path}/shared_network_clean.tsv",
}

for name, input_path in networks.items():
    df = pd.read_csv(input_path, sep="\t")

    G = nx.from_pandas_edgelist(df, source="source", target="target")

    degree = nx.degree_centrality(G)
    betweenness = nx.betweenness_centrality(G)

    hub_df = pd.DataFrame({
        "gene": list(degree.keys()),
        "degree_centrality": list(degree.values()),
        "betweenness": list(betweenness.values())
    })

    output_path = f"{base_path}/{name}_real_hub_genes.csv"
    hub_df.sort_values("degree_centrality", ascending=False).to_csv(output_path, index=False)