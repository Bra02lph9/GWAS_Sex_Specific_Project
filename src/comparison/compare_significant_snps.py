from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

FILTERED_DIR = PROJECT_ROOT / "1_data" / "filtered"
TABLES_DIR = PROJECT_ROOT / "4_results" / "tables"

PHENOTYPE = "CAD"


def load_snps(file_path: Path) -> pd.DataFrame:
    if not file_path.exists():
        raise FileNotFoundError(f"Missing file: {file_path}")

    df = pd.read_csv(file_path, sep="\t", low_memory=False)

    if "SNP" not in df.columns:
        raise ValueError(f"SNP column missing in {file_path.name}")

    df["SNP"] = (
        df["SNP"]
        .astype(str)
        .str.strip()
    )

    df = df[
        (df["SNP"] != "")
        & (df["SNP"].str.lower() != "nan")
    ].copy()

    before_duplicates = len(df)
    df = df.drop_duplicates(subset=["SNP"]).copy()
    removed_duplicates = before_duplicates - len(df)

    print(f"Loaded: {file_path.name}")
    print(f"Rows after cleaning: {len(df):,}")
    print(f"Duplicate SNPs removed: {removed_duplicates:,}")
    print("-" * 60)

    return df


def compare_snps(male_df: pd.DataFrame, female_df: pd.DataFrame):
    male_snps = set(male_df["SNP"])
    female_snps = set(female_df["SNP"])

    shared = male_snps & female_snps
    male_only = male_snps - female_snps
    female_only = female_snps - male_snps

    return shared, male_only, female_only


def save_group(df: pd.DataFrame, snps: set, output_file: Path) -> pd.DataFrame:
    result = df[df["SNP"].isin(snps)].copy()

    if "p_value" in result.columns:
        result["p_value"] = pd.to_numeric(result["p_value"], errors="coerce")
        result = result.sort_values(by="p_value", ascending=True)
    else:
        result = result.sort_values(by="SNP")

    result.to_csv(output_file, sep="\t", index=False)

    print(f"Saved: {output_file.name}")
    print(f"Rows: {len(result):,}")
    print("-" * 60)

    return result


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    male_file = FILTERED_DIR / f"{PHENOTYPE}_male_significant_snps.tsv"
    female_file = FILTERED_DIR / f"{PHENOTYPE}_female_significant_snps.tsv"

    male_df = load_snps(male_file)
    female_df = load_snps(female_file)

    shared, male_only, female_only = compare_snps(male_df, female_df)

    shared_df = save_group(
        male_df,
        shared,
        TABLES_DIR / f"{PHENOTYPE}_shared_snps.tsv"
    )

    male_only_df = save_group(
        male_df,
        male_only,
        TABLES_DIR / f"{PHENOTYPE}_male_only_snps.tsv"
    )

    female_only_df = save_group(
        female_df,
        female_only,
        TABLES_DIR / f"{PHENOTYPE}_female_only_snps.tsv"
    )

    summary = pd.DataFrame({
        "group": [
            "male_total",
            "female_total",
            "shared",
            "male_only",
            "female_only",
        ],
        "count": [
            len(male_df),
            len(female_df),
            len(shared_df),
            len(male_only_df),
            len(female_only_df),
        ],
    })

    summary_file = TABLES_DIR / f"{PHENOTYPE}_summary.tsv"

    summary.to_csv(
        summary_file,
        sep="\t",
        index=False
    )

    print(summary)
    print("-" * 60)
    print(f"Summary saved to: {summary_file}")
    print("SNP comparison completed successfully.")


if __name__ == "__main__":
    main()
