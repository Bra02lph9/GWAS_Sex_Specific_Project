from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

FUMA_DIR = PROJECT_ROOT / "3_tools_results" / "fuma"
GENE_LISTS_DIR = PROJECT_ROOT / "1_data" / "gene_lists"
TABLES_DIR = PROJECT_ROOT / "4_results" / "tables"

PHENOTYPE = "CAD"
P_THRESHOLD = 0.05


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
    result = df[df["P"] < P_THRESHOLD].copy()
    result = result.sort_values(by="P", ascending=True)
    return result


def compare_gene_sets(male_df: pd.DataFrame, female_df: pd.DataFrame):
    male_genes = set(male_df["GENE"])
    female_genes = set(female_df["GENE"])

    shared = male_genes & female_genes
    male_only = male_genes - female_genes
    female_only = female_genes - male_genes

    return shared, male_only, female_only


def save_gene_list(genes: set, output_file: Path) -> None:
    result = pd.DataFrame({
        "GENE": sorted(list(genes))
    })

    result.to_csv(
        output_file,
        sep="\t",
        index=False
    )

    print(f"Saved: {output_file.name}")
    print(f"Genes: {len(result):,}")
    print("-" * 60)


def save_gene_table(df: pd.DataFrame, genes: set, output_file: Path) -> pd.DataFrame:
    result = df[df["GENE"].isin(genes)].copy()

    if "P" in result.columns:
        result = result.sort_values(by="P", ascending=True)
    else:
        result = result.sort_values(by="GENE")

    result.to_csv(
        output_file,
        sep="\t",
        index=False
    )

    return result


def main() -> None:
    GENE_LISTS_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    male_file = FUMA_DIR / "male" / "magma.genes.out"
    female_file = FUMA_DIR / "female" / "magma.genes.out"

    male_df = load_magma_genes(male_file)
    female_df = load_magma_genes(female_file)

    male_sig = filter_significant_genes(male_df)
    female_sig = filter_significant_genes(female_df)

    shared, male_only, female_only = compare_gene_sets(male_sig, female_sig)

    save_gene_list(
        shared,
        GENE_LISTS_DIR / "shared_genes.txt"
    )

    save_gene_list(
        male_only,
        GENE_LISTS_DIR / "male_specific_genes.txt"
    )

    save_gene_list(
        female_only,
        GENE_LISTS_DIR / "female_specific_genes.txt"
    )

    save_gene_table(
        male_sig,
        shared,
        TABLES_DIR / f"{PHENOTYPE}_shared_gene_table.tsv"
    )

    save_gene_table(
        male_sig,
        male_only,
        TABLES_DIR / f"{PHENOTYPE}_male_specific_gene_table.tsv"
    )

    save_gene_table(
        female_sig,
        female_only,
        TABLES_DIR / f"{PHENOTYPE}_female_specific_gene_table.tsv"
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
        ]
    })

    summary_file = TABLES_DIR / "gene_summary.tsv"

    summary.to_csv(
        summary_file,
        sep="\t",
        index=False
    )

    print(summary)
    print("-" * 60)
    print(f"Summary saved to: {summary_file}")
    print("Gene extraction and comparison completed successfully.")


if __name__ == "__main__":
    main()
