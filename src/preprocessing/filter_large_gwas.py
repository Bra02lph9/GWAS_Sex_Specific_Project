from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

RAW_DIR = PROJECT_ROOT / "1_data" / "raw"
FILTERED_DIR = PROJECT_ROOT / "1_data" / "filtered"

PHENOTYPE = "CAD"
P_VALUE_THRESHOLD = 5e-8
CHUNK_SIZE = 500_000


USEFUL_COLUMNS = [
    "SNP",
    "chromosome",
    "position",
    "effect_allele",
    "effect_allele_frequency",
    "beta",
    "standard_error",
    "p_value",
    "sample_size",
]


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = (
        df.columns
        .str.strip()
        .str.lower()
        .str.replace(" ", "_", regex=False)
        .str.replace("-", "_", regex=False)
    )

    rename_map = {
        "rsid": "SNP",
        "snp": "SNP",
        "chr": "chromosome",
        "chrom": "chromosome",
        "bp": "position",
        "pos": "position",
        "reference_allele": "effect_allele",
        "a1": "effect_allele",
        "effect_allele": "effect_allele",
        "eaf": "effect_allele_frequency",
        "effect_allele_frequency": "effect_allele_frequency",
        "beta": "beta",
        "se": "standard_error",
        "standard_error": "standard_error",
        "p": "p_value",
        "p_value": "p_value",
        "n": "sample_size",
        "sample_size": "sample_size",
    }

    return df.rename(columns=rename_map)


def validate_required_columns(df: pd.DataFrame, file_name: str) -> None:
    required_columns = ["SNP", "p_value"]

    missing_columns = [col for col in required_columns if col not in df.columns]

    if missing_columns:
        raise ValueError(
            f"Missing required columns in {file_name}: {missing_columns}"
        )


def keep_available_columns(df: pd.DataFrame) -> pd.DataFrame:
    available_columns = [col for col in USEFUL_COLUMNS if col in df.columns]
    return df[available_columns].copy()


def clean_chunk(chunk: pd.DataFrame, file_name: str) -> pd.DataFrame:
    chunk = normalize_columns(chunk)

    validate_required_columns(chunk, file_name)

    chunk = keep_available_columns(chunk)

    chunk["SNP"] = (
        chunk["SNP"]
        .astype(str)
        .str.strip()
    )

    chunk["p_value"] = pd.to_numeric(
        chunk["p_value"],
        errors="coerce"
    )

    if "chromosome" in chunk.columns:
        chunk["chromosome"] = pd.to_numeric(
            chunk["chromosome"],
            errors="coerce"
        )

    if "position" in chunk.columns:
        chunk["position"] = pd.to_numeric(
            chunk["position"],
            errors="coerce"
        )

    if "beta" in chunk.columns:
        chunk["beta"] = pd.to_numeric(
            chunk["beta"],
            errors="coerce"
        )

    if "standard_error" in chunk.columns:
        chunk["standard_error"] = pd.to_numeric(
            chunk["standard_error"],
            errors="coerce"
        )

    if "effect_allele_frequency" in chunk.columns:
        chunk["effect_allele_frequency"] = pd.to_numeric(
            chunk["effect_allele_frequency"],
            errors="coerce"
        )

    if "sample_size" in chunk.columns:
        chunk["sample_size"] = pd.to_numeric(
            chunk["sample_size"],
            errors="coerce"
        )

    chunk = chunk.dropna(subset=["SNP", "p_value"])

    chunk = chunk[
        (chunk["SNP"] != "")
        & (chunk["SNP"].str.lower() != "nan")
        & (chunk["p_value"] > 0)
        & (chunk["p_value"] <= 1)
    ]

    return chunk


def filter_significant_snps(input_file: Path, output_file: Path) -> None:
    FILTERED_DIR.mkdir(parents=True, exist_ok=True)

    first_chunk = True
    total_rows = 0
    valid_rows = 0
    total_significant = 0

    print(f"Processing: {input_file.name}")
    print(f"P-value threshold: {P_VALUE_THRESHOLD}")

    for chunk in pd.read_csv(
        input_file,
        sep=",",
        chunksize=CHUNK_SIZE,
        low_memory=False,
    ):
        total_rows += len(chunk)

        chunk = clean_chunk(chunk, input_file.name)
        valid_rows += len(chunk)

        significant = chunk[chunk["p_value"] <= P_VALUE_THRESHOLD].copy()

        if significant.empty:
            continue

        total_significant += len(significant)

        significant.to_csv(
            output_file,
            sep="\t",
            index=False,
            mode="w" if first_chunk else "a",
            header=first_chunk,
        )

        first_chunk = False

    if total_significant == 0:
        pd.DataFrame(columns=USEFUL_COLUMNS).to_csv(
            output_file,
            sep="\t",
            index=False,
        )

    print(f"Total rows processed: {total_rows:,}")
    print(f"Valid rows after cleaning: {valid_rows:,}")
    print(f"Significant SNPs kept: {total_significant:,}")
    print(f"Saved to: {output_file}")
    print("-" * 60)


def main() -> None:
    files = {
        "male": RAW_DIR / f"{PHENOTYPE}_male.tsv",
        "female": RAW_DIR / f"{PHENOTYPE}_female.tsv",
    }

    for sex, file_path in files.items():
        if not file_path.exists():
            raise FileNotFoundError(f"Missing file: {file_path}")

        output_file = FILTERED_DIR / f"{PHENOTYPE}_{sex}_significant_snps.tsv"

        filter_significant_snps(file_path, output_file)


if __name__ == "__main__":
    main()
