import pandas as pd
from pathlib import Path


INPUT = Path("3_tools_results/string/female_string_network.tsv")

OUTPUT = Path(
    "3_tools_results/cytoscape/female_network_clean.tsv"
)

THRESHOLD = 0.7


df = pd.read_csv(INPUT, sep="\t")

df = df[["#node1", "node2", "combined_score"]].copy()

df.columns = ["source", "target", "combined_score"]

df["combined_score"] = pd.to_numeric(
    df["combined_score"],
    errors="coerce"
)

df = df[df["combined_score"] >= THRESHOLD].copy()

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

df.to_csv(OUTPUT, sep="\t", index=False)

print(df.head())
print(f"Edges kept: {len(df)}")
print(f"Saved to: {OUTPUT}")
