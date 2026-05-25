from pathlib import Path
import pandas as pd

from src.utils.cli import parse_phenotype
from src.utils.config import get_paths, create_output_dirs, CHUNK_SIZE


FUMA_COLUMNS = [
    "SNP",
    "CHR",
    "BP",
    "A1",
    "BETA",
    "SE",
    "P",
    "N",
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
        "chr": "CHR",
        "chromosome": "CHR",
        "chrom": "CHR",
        "bp": "BP",
        "position": "BP",
        "pos": "BP",
        "reference_allele": "A1",
        "effect_allele": "A1",
        "a1": "A1",
        "beta": "BETA",
        "se": "SE",
        "standard_error": "SE",
        "p_value": "P",
        "p": "P",
        "n": "N",
        "sample_size": "N",
    }

    return df.rename(columns=rename_map)


def validate_fuma_columns(df: pd.DataFrame, file_name: str) -> None:
    missing_columns = [col for col in FUMA_COLUMNS if col not in df.columns]

    if missing_columns:
        raise ValueError(
            f"Missing required FUMA columns in {file_name}: {missing_columns}"
        )


def clean_fuma_chunk(chunk: pd.DataFrame, file_name: str) -> pd.DataFrame:
    chunk = normalize_columns(chunk)
    validate_fuma_columns(chunk, file_name)

    chunk = chunk[FUMA_COLUMNS].copy()

    chunk["SNP"] = chunk["SNP"].astype(str).str.strip()
    chunk["A1"] = chunk["A1"].astype(str).str.upper().str.strip()

    numeric_columns = ["CHR", "BP", "BETA", "SE", "P", "N"]

    for col in numeric_columns:
        chunk[col] = pd.to_numeric(chunk[col], errors="coerce")

    chunk = chunk.dropna(
        subset=["SNP", "CHR", "BP", "A1", "BETA", "SE", "P"]
    )

    chunk = chunk[
        (chunk["SNP"] != "")
        & (chunk["SNP"].str.lower() != "nan")
        & (chunk["A1"] != "")
        & (chunk["A1"].str.lower() != "nan")
        & (chunk["P"] > 0)
        & (chunk["P"] <= 1)
        & (chunk["CHR"] >= 1)
        & (chunk["CHR"] <= 22)
        & (chunk["BP"] > 0)
    ]

    chunk["CHR"] = chunk["CHR"].astype(int)
    chunk["BP"] = chunk["BP"].astype(int)

    return chunk


def prepare_fuma_file(input_file: Path, output_file: Path) -> None:
    if not input_file.exists():
        raise FileNotFoundError(f"Missing raw GWAS file: {input_file}")

    output_file.parent.mkdir(parents=True, exist_ok=True)

    first_chunk = True
    total_rows = 0
    kept_rows = 0

    print(f"Preparing FUMA input: {input_file}")

    for chunk in pd.read_csv(
        input_file,
        sep=",",
        chunksize=CHUNK_SIZE,
        low_memory=False,
    ):
        total_rows += len(chunk)

        chunk = clean_fuma_chunk(chunk, input_file.name)
        kept_rows += len(chunk)

        if chunk.empty:
            continue

        chunk.to_csv(
            output_file,
            sep="\t",
            index=False,
            mode="w" if first_chunk else "a",
            header=first_chunk,
        )

        first_chunk = False

    if kept_rows == 0:
        pd.DataFrame(columns=FUMA_COLUMNS).to_csv(
            output_file,
            sep="\t",
            index=False,
        )

    print(f"Total rows processed: {total_rows:,}")
    print(f"Rows kept for FUMA: {kept_rows:,}")
    print(f"Saved to: {output_file}")
    print("-" * 60)


def main() -> None:
    phenotype = parse_phenotype()
    paths = get_paths(phenotype)
    create_output_dirs(paths)

    files = {
        "male": paths["male_raw"],
        "female": paths["female_raw"],
    }

    outputs = {
        "male": paths["fuma_male_dir"] / f"{phenotype}_male_fuma_input.tsv",
        "female": paths["fuma_female_dir"] / f"{phenotype}_female_fuma_input.tsv",
    }

    for sex, input_file in files.items():
        prepare_fuma_file(input_file, outputs[sex])

    print(f"FUMA input preparation completed successfully for: {phenotype}")


if __name__ == "__main__":
    main()
