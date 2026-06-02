from pathlib import Path
from typing import Optional

import pandas as pd

from src.utils.cli import parse_phenotype
from src.utils.config import get_paths, create_output_dirs


SEXES = ("male", "female")

EQTL_REQUIRED_COLUMNS = [
    "gene",
    "group",
    "variant_id",
    "snp",
    "p_value",
    "nes",
    "tissue",
]

LOCI_REQUIRED_COLUMNS = [
    "GenomicLocus",
    "chr",
    "start",
    "end",
]


def read_tsv(file_path: Path) -> pd.DataFrame:
    if not file_path.exists():
        raise FileNotFoundError(f"Missing file: {file_path}")

    try:
        df = pd.read_csv(file_path, sep="\t", low_memory=False)

        if len(df.columns) > 1:
            df.columns = df.columns.str.strip()
            return df

    except Exception:
        pass

    df = pd.read_csv(file_path, sep=r"\s+", engine="python")
    df.columns = df.columns.str.strip()

    return df


def validate_required_columns(
    df: pd.DataFrame,
    required_columns: list[str],
    file_name: str,
) -> None:
    """
    Ensure that all required columns are present in a DataFrame.
    """

    missing_columns = [col for col in required_columns if col not in df.columns]

    if missing_columns:
        raise ValueError(
            f"Missing required columns in {file_name}: {missing_columns}"
        )


def load_eqtl_results(file_path: Path) -> pd.DataFrame:
    """
    Load and clean GTEx eQTL results.
    """

    df = read_tsv(file_path)
    validate_required_columns(df, EQTL_REQUIRED_COLUMNS, file_path.name)

    df = df[EQTL_REQUIRED_COLUMNS].copy()

    text_columns = ["gene", "group", "variant_id", "snp", "tissue"]
    numeric_columns = ["p_value", "nes"]

    for column in text_columns:
        df[column] = df[column].astype(str).str.strip()

    df["group"] = df["group"].str.lower()

    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    df = df.dropna(
        subset=[
            "gene",
            "group",
            "variant_id",
            "snp",
            "p_value",
            "nes",
        ]
    )

    if df.empty:
        raise ValueError(
           "No valid eQTL rows remained after cleaning. "
           "Check separators, missing values, or numeric columns p_value/nes."
        )

    return df.reset_index(drop=True)


def load_significant_snps(
    paths: dict[str, Path],
    phenotype: str,
) -> dict[str, set[str]]:
    """
    Load male and female significant SNP sets.
    """

    snp_sets: dict[str, set[str]] = {}

    for sex in SEXES:
        file_path = paths["filtered_dir"] / f"{phenotype}_{sex}_significant_snps.tsv"
        df = read_tsv(file_path)

        validate_required_columns(df, ["SNP"], file_path.name)

        df["SNP"] = df["SNP"].astype(str).str.strip()
        df = df[(df["SNP"] != "") & (df["SNP"].str.lower() != "nan")]

        snp_sets[sex] = set(df["SNP"])

    return snp_sets


def load_loci(
    paths: dict[str, Path],
    phenotype: str,
) -> dict[str, pd.DataFrame]:

    loci_sets: dict[str, pd.DataFrame] = {}

    for sex in SEXES:
        file_path = paths["loci_dir"] / f"{phenotype}_{sex}_loci.tsv"
        df = read_tsv(file_path)

        validate_required_columns(df, LOCI_REQUIRED_COLUMNS, file_path.name)

        df = df[LOCI_REQUIRED_COLUMNS].copy()

        df = df.rename(
            columns={
                "GenomicLocus": "locus_id",
                "chr": "chromosome",
            }
        )

        numeric_columns = ["locus_id", "chromosome", "start", "end"]

        for column in numeric_columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

        df = df.dropna(subset=numeric_columns)

        df = df[
            (df["chromosome"].between(1, 22))
            & (df["start"] > 0)
            & (df["end"] >= df["start"])
        ].copy()

        df["locus_id"] = df["locus_id"].astype(int)
        df["chromosome"] = df["chromosome"].astype(int)
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)

        df["locus_id"] = df["locus_id"].apply(
            lambda x: f"{phenotype}_{sex}_LOCUS_{x:04d}"
        )

        loci_sets[sex] = df.reset_index(drop=True)

    return loci_sets


def parse_gtex_variant_id(
    variant_id: str,
) -> tuple[Optional[int], Optional[int]]:
    """
    Parse a GTEx variant ID.

    Example:
        chr14_105499151_G_C_b38

    Returns:
        chromosome, position
    """

    try:
        parts = variant_id.split("_")

        if len(parts) < 2:
            return None, None

        chromosome = parts[0].replace("chr", "")
        position = int(parts[1])

        if chromosome in {"X", "Y", "MT", "M"}:
            return None, None

        return int(chromosome), position

    except (ValueError, IndexError):
        return None, None


def find_overlapping_loci(
    chromosome: Optional[int],
    position: Optional[int],
    loci_df: pd.DataFrame,
) -> str:
    """
    Find loci overlapping a given genomic position.
    """

    if chromosome is None or position is None:
        return ""

    overlaps = loci_df[
        (loci_df["chromosome"] == chromosome)
        & (loci_df["start"] <= position)
        & (loci_df["end"] >= position)
    ]

    if overlaps.empty:
        return ""

    return ";".join(overlaps["locus_id"].astype(str).tolist())


def assign_validation_level(row: pd.Series) -> str:
    """
    Assign validation category based on exact SNP and locus-level overlap.
    """

    exact_match = row["exact_snp_match_male"] or row["exact_snp_match_female"]
    locus_overlap = row["locus_overlap_male"] or row["locus_overlap_female"]

    if exact_match:
        return "strong_exact_snp_match"

    if locus_overlap:
        return "good_locus_overlap"

    return "supportive_only_no_gwas_overlap"


def validate_eqtl_results(
    eqtl_df: pd.DataFrame,
    snp_sets: dict[str, set[str]],
    loci_sets: dict[str, pd.DataFrame],
) -> pd.DataFrame:

    if eqtl_df.empty:
        return pd.DataFrame()

    records = []

    for _, row in eqtl_df.iterrows():
        chromosome, position = parse_gtex_variant_id(row["variant_id"])

        exact_male = row["snp"] in snp_sets["male"]
        exact_female = row["snp"] in snp_sets["female"]

        male_loci = find_overlapping_loci(chromosome, position, loci_sets["male"])
        female_loci = find_overlapping_loci(chromosome, position, loci_sets["female"])

        record = row.to_dict()

        record.update(
            {
                "chromosome": chromosome,
                "position": position,
                "exact_snp_match_male": exact_male,
                "exact_snp_match_female": exact_female,
                "locus_overlap_male": bool(male_loci),
                "locus_overlap_female": bool(female_loci),
                "matched_male_loci": male_loci,
                "matched_female_loci": female_loci,
            }
        )

        records.append(record)

    validated_df = pd.DataFrame(records)

    if validated_df.empty:
        return validated_df

    validated_df["validation_level"] = validated_df.apply(
        assign_validation_level,
        axis=1,
    )

    return validated_df


def get_eqtl_file_path(
    paths: dict[str, Path],
    phenotype: str,
) -> Path:
    """
    Return GTEx eQTL file path.
    """

    if "gtex_dir" in paths:
        return paths["gtex_dir"] / f"{phenotype}_eqtl_results.tsv"

    return (
        Path("3_tools_results")
        / "gtex"
        / phenotype
        / f"{phenotype}_eqtl_results.tsv"
    )


def save_validation_outputs(
    validated_df: pd.DataFrame,
    paths: dict[str, Path],
    phenotype: str,
) -> tuple[Path, Path]:
    """
    Save validated eQTL table and summary table.
    """

    output_file = paths["tables_dir"] / f"{phenotype}_eqtl_validated.tsv"
    summary_file = paths["tables_dir"] / f"{phenotype}_eqtl_validation_summary.tsv"

    validated_df.to_csv(output_file, sep="\t", index=False)

    summary_df = (
        validated_df["validation_level"]
        .value_counts()
        .reset_index()
    )

    summary_df.columns = ["validation_level", "count"]
    summary_df.to_csv(summary_file, sep="\t", index=False)

    return output_file, summary_file


def main() -> None:
    phenotype = parse_phenotype()
    paths = get_paths(phenotype)

    create_output_dirs(paths)

    eqtl_file = get_eqtl_file_path(paths, phenotype)

    print("=" * 80)
    print(f"Validating GTEx eQTL results for phenotype: {phenotype}")
    print(f"Input eQTL file: {eqtl_file}")
    print("=" * 80)

    eqtl_df = load_eqtl_results(eqtl_file)
    snp_sets = load_significant_snps(paths, phenotype)
    loci_sets = load_loci(paths, phenotype)

    validated_df = validate_eqtl_results(
        eqtl_df=eqtl_df,
        snp_sets=snp_sets,
        loci_sets=loci_sets,
    )

    output_file, summary_file = save_validation_outputs(
        validated_df=validated_df,
        paths=paths,
        phenotype=phenotype,
    )

    summary_df = read_tsv(summary_file)

    print(summary_df)
    print("-" * 80)
    print(f"Saved validated eQTL table: {output_file}")
    print(f"Saved summary table       : {summary_file}")


if __name__ == "__main__":
    main()