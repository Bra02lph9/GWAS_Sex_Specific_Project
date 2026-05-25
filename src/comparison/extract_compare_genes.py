from pathlib import Path
import pandas as pd

from src.utils.cli import parse_phenotype
from src.utils.config import (
    get_paths,
    create_output_dirs,
    MAGMA_P_THRESHOLD,
)


def load_magma_genes(file_path: Path) -> pd.DataFrame:
    if not file_path.exists():
        raise FileNotFoundError(f"Missing MAGMA file: {file_path}")

    df = pd.read_csv(file_path, sep="\t", low_memory=False)

    required_columns = ["GENE", "P"]
    missing_columns = [col for col in required_columns if col not in df.columns]

    if missing_columns:
        raise ValueError(
            f"Invalid MAGMA file: {file_path.name}. Missing columns: {missing_columns}"
        )

    df["GENE"] = df["GENE"].astype(str).str.strip()
    df["P"] = pd.to_numeric(df["P"], errors="coerce")

    df = df.dropna(subset=["GENE", "P"])

    df = df[
        (df["GENE"] != "")
        & (df["GENE"].str.lower() != "nan")
        & (df["P"] > 0)
        & (df["P"] <= 1)
    ].copy()

    before_duplicates = len(df)

    df = df.sort_values(by="P", ascending=True)
    df = df.drop_duplicates(subset=["GENE"], keep="first").copy()

    removed_duplicates = before_duplicates - len(df)

    print(f"Loaded: {file_path}")
    print(f"Valid genes: {len(df):,}")
    print(f"Duplicate genes removed: {removed_duplicates:,}")
    print("-" * 60)

    return df


def filter_significant_genes(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["P"] < MAGMA_P_THRESHOLD].sort_values(by="P").copy()


def compare_gene_sets(male_df: pd.DataFrame, female_df: pd.DataFrame):
    male_genes = set(male_df["GENE"])
    female_genes = set(female_df["GENE"])

    shared = male_genes & female_genes
    male_only = male_genes - female_genes
    female_only = female_genes - male_genes

    return shared, male_only, female_only


def save_gene_list(genes: set, output_file: Path) -> None:
    result = pd.DataFrame({"GENE": sorted(genes)})
    result.to_csv(output_file, sep="\t", index=False)

    print(f"Saved: {output_file.name}")
    print(f"Genes: {len(result):,}")
    print("-" * 60)


def save_gene_table(df: pd.DataFrame, genes: set, output_file: Path) -> pd.DataFrame:
    result = df[df["GENE"].isin(genes)].copy()

    if "P" in result.columns:
        result = result.sort_values(by="P", ascending=True)
    else:
        result = result.sort_values(by="GENE")

    result.to_csv(output_file, sep="\t", index=False)
    return result


def main() -> None:
    phenotype = parse_phenotype()
    paths = get_paths(phenotype)
    create_output_dirs(paths)

    male_file = paths["fuma_male_dir"] / "magma.genes.out"
    female_file = paths["fuma_female_dir"] / "magma.genes.out"

    male_df = load_magma_genes(male_file)
    female_df = load_magma_genes(female_file)

    male_sig = filter_significant_genes(male_df)
    female_sig = filter_significant_genes(female_df)

    shared, male_only, female_only = compare_gene_sets(male_sig, female_sig)

    save_gene_list(
        shared,
        paths["gene_lists_dir"] / "shared_genes.txt",
    )

    save_gene_list(
        male_only,
        paths["gene_lists_dir"] / "male_specific_genes.txt",
    )

    save_gene_list(
        female_only,
        paths["gene_lists_dir"] / "female_specific_genes.txt",
    )

    save_gene_table(
        male_sig,
        shared,
        paths["tables_dir"] / f"{phenotype}_shared_gene_table.tsv",
    )

    save_gene_table(
        male_sig,
        male_only,
        paths["tables_dir"] / f"{phenotype}_male_specific_gene_table.tsv",
    )

    save_gene_table(
        female_sig,
        female_only,
        paths["tables_dir"] / f"{phenotype}_female_specific_gene_table.tsv",
    )

    summary = pd.DataFrame({
        "group": [
            "male_significant_genes",
            "female_significant_genes",
            "shared_genes",
            "male_specific_genes",
            "female_specific_genes",
        ],
        "count": [
            len(male_sig),
            len(female_sig),
            len(shared),
            len(male_only),
            len(female_only),
        ],
    })

    summary_file = paths["tables_dir"] / f"{phenotype}_gene_summary.tsv"
    summary.to_csv(summary_file, sep="\t", index=False)

    print(summary)
    print("-" * 60)
    print(f"Summary saved to: {summary_file}")
    print(f"Gene extraction and comparison completed successfully for: {phenotype}")


if __name__ == "__main__":
    main()
