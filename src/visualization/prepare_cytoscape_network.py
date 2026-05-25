from pathlib import Path
import pandas as pd

from src.utils.cli import parse_phenotype
from src.utils.config import get_paths, create_output_dirs


THRESHOLD = 0.7
NETWORK_GROUPS = ["female", "male", "shared"]


def prepare_network(input_file: Path, output_file: Path) -> None:
    if not input_file.exists():
        raise FileNotFoundError(f"Missing STRING network file: {input_file}")

    df = pd.read_csv(input_file, sep="\t", low_memory=False)

    required_columns = ["#node1", "node2", "combined_score"]
    missing_columns = [col for col in required_columns if col not in df.columns]

    if missing_columns:
        raise ValueError(
            f"Missing columns in {input_file.name}: {missing_columns}"
        )

    df = df[required_columns].copy()
    df.columns = ["source", "target", "combined_score"]

    df["source"] = df["source"].astype(str).str.strip()
    df["target"] = df["target"].astype(str).str.strip()
    df["combined_score"] = pd.to_numeric(df["combined_score"], errors="coerce")

    df = df.dropna(subset=["source", "target", "combined_score"])

    df = df[
        (df["source"] != "")
        & (df["target"] != "")
        & (df["source"].str.lower() != "nan")
        & (df["target"].str.lower() != "nan")
        & (df["source"] != df["target"])
        & (df["combined_score"] >= THRESHOLD)
    ].copy()

    df = df.drop_duplicates(subset=["source", "target"])

    output_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, sep="\t", index=False)

    print(df.head())
    print(f"Edges kept: {len(df):,}")
    print(f"Saved to: {output_file}")
    print("-" * 60)


def main() -> None:
    phenotype = parse_phenotype()
    paths = get_paths(phenotype)
    create_output_dirs(paths)

    for group in NETWORK_GROUPS:
        input_file = paths["string_dir"] / f"{group}_string_network.tsv"
        output_file = paths["cytoscape_dir"] / f"{group}_network_clean.tsv"

        print(f"Preparing {phenotype} {group} Cytoscape network...")
        prepare_network(input_file, output_file)

    print(f"Cytoscape network preparation completed for: {phenotype}")


if __name__ == "__main__":
    main()
