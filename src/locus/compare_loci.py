from pathlib import Path
from typing import Any

import pandas as pd

from src.utils.cli import parse_phenotype
from src.utils.config import get_paths, create_output_dirs


REQUIRED_FUMA_COLUMNS = [
    "GenomicLocus",
    "rsID",
    "chr",
    "pos",
    "p",
    "start",
    "end",
    "nSNPs",
    "LeadSNPs",
]

STANDARD_COLUMNS = [
    "locus_id",
    "chromosome",
    "start",
    "end",
    "lead_snp",
    "lead_snp_position",
    "lead_p_value",
    "n_snps",
    "lead_snps",
]

SHARED_COLUMNS = [
    "chromosome",
    "shared_start",
    "shared_end",
    "male_locus_id",
    "female_locus_id",
    "male_lead_snp",
    "female_lead_snp",
    "male_lead_p_value",
    "female_lead_p_value",
    "male_start",
    "male_end",
    "female_start",
    "female_end",
]


def load_loci(file_path: Path, phenotype: str, sex: str) -> pd.DataFrame:
    """
    Load FUMA-like loci file and convert it to a standard internal format.
    """

    if not file_path.exists():
        raise FileNotFoundError(f"Missing loci file: {file_path}")

    df = pd.read_csv(file_path, sep="\t", low_memory=False)

    missing_columns = [
        col for col in REQUIRED_FUMA_COLUMNS if col not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Missing required FUMA columns in {file_path.name}: {missing_columns}"
        )

    df = df[REQUIRED_FUMA_COLUMNS].copy()

    df = df.rename(
        columns={
            "GenomicLocus": "original_locus_id",
            "chr": "chromosome",
            "rsID": "lead_snp",
            "pos": "lead_snp_position",
            "p": "lead_p_value",
            "nSNPs": "n_snps",
            "LeadSNPs": "lead_snps",
        }
    )

    numeric_columns = [
        "original_locus_id",
        "chromosome",
        "start",
        "end",
        "lead_snp_position",
        "lead_p_value",
        "n_snps",
    ]

    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    df = df.dropna(
        subset=[
            "original_locus_id",
            "chromosome",
            "start",
            "end",
            "lead_snp_position",
            "lead_p_value",
            "n_snps",
        ]
    )

    df = df[
        (df["chromosome"].between(1, 22))
        & (df["start"] > 0)
        & (df["end"] >= df["start"])
        & (df["lead_snp_position"] > 0)
        & (df["lead_p_value"].between(0, 1, inclusive="right"))
    ].copy()

    df["original_locus_id"] = df["original_locus_id"].astype(int)
    df["chromosome"] = df["chromosome"].astype(int)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["lead_snp_position"] = df["lead_snp_position"].astype(int)
    df["n_snps"] = df["n_snps"].astype(int)

    df["lead_snp"] = df["lead_snp"].astype(str).str.strip()
    df["lead_snps"] = df["lead_snps"].astype(str).str.strip()

    df["locus_id"] = df["original_locus_id"].apply(
        lambda x: f"{phenotype}_{sex}_LOCUS_{x:04d}"
    )

    return df[STANDARD_COLUMNS].sort_values(
        ["chromosome", "start", "end"]
    ).reset_index(drop=True)


def loci_overlap(locus_a: pd.Series, locus_b: pd.Series) -> bool:
    if locus_a["chromosome"] != locus_b["chromosome"]:
        return False

    return locus_a["start"] <= locus_b["end"] and locus_b["start"] <= locus_a["end"]


def build_shared_locus_record(
    male_locus: pd.Series,
    female_locus: pd.Series,
) -> dict[str, Any]:

    return {
        "chromosome": int(male_locus["chromosome"]),
        "shared_start": int(max(male_locus["start"], female_locus["start"])),
        "shared_end": int(min(male_locus["end"], female_locus["end"])),
        "male_locus_id": male_locus["locus_id"],
        "female_locus_id": female_locus["locus_id"],
        "male_lead_snp": male_locus["lead_snp"],
        "female_lead_snp": female_locus["lead_snp"],
        "male_lead_p_value": male_locus["lead_p_value"],
        "female_lead_p_value": female_locus["lead_p_value"],
        "male_start": int(male_locus["start"]),
        "male_end": int(male_locus["end"]),
        "female_start": int(female_locus["start"]),
        "female_end": int(female_locus["end"]),
    }


def compare_loci(
    male_df: pd.DataFrame,
    female_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:

    shared_male_indices: set[int] = set()
    shared_female_indices: set[int] = set()
    shared_records: list[dict[str, Any]] = []

    female_by_chromosome = {
        chromosome: chrom_df
        for chromosome, chrom_df in female_df.groupby("chromosome")
    }

    for male_idx, male_locus in male_df.iterrows():
        same_chrom_female = female_by_chromosome.get(
            male_locus["chromosome"],
            pd.DataFrame(),
        )

        for female_idx, female_locus in same_chrom_female.iterrows():
            if not loci_overlap(male_locus, female_locus):
                continue

            shared_male_indices.add(male_idx)
            shared_female_indices.add(female_idx)

            shared_records.append(
                build_shared_locus_record(male_locus, female_locus)
            )

    shared_df = pd.DataFrame(shared_records, columns=SHARED_COLUMNS)

    male_only_df = male_df.loc[
        ~male_df.index.isin(shared_male_indices)
    ].copy()

    female_only_df = female_df.loc[
        ~female_df.index.isin(shared_female_indices)
    ].copy()

    return shared_df, male_only_df, female_only_df


def build_summary(
    male_df: pd.DataFrame,
    female_df: pd.DataFrame,
    shared_df: pd.DataFrame,
    male_only_df: pd.DataFrame,
    female_only_df: pd.DataFrame,
) -> pd.DataFrame:

    return pd.DataFrame(
        {
            "group": [
                "male_total_loci",
                "female_total_loci",
                "shared_loci",
                "male_only_loci",
                "female_only_loci",
            ],
            "count": [
                len(male_df),
                len(female_df),
                len(shared_df),
                len(male_only_df),
                len(female_only_df),
            ],
        }
    )


def save_results(
    phenotype: str,
    paths: dict[str, Path],
    shared_df: pd.DataFrame,
    male_only_df: pd.DataFrame,
    female_only_df: pd.DataFrame,
    summary_df: pd.DataFrame,
) -> None:

    output_files = {
        "shared_loci": paths["tables_dir"] / f"{phenotype}_shared_loci.tsv",
        "male_only_loci": paths["tables_dir"] / f"{phenotype}_male_only_loci.tsv",
        "female_only_loci": paths["tables_dir"] / f"{phenotype}_female_only_loci.tsv",
        "summary": paths["tables_dir"] / f"{phenotype}_locus_summary.tsv",
    }

    shared_df.to_csv(output_files["shared_loci"], sep="\t", index=False)
    male_only_df.to_csv(output_files["male_only_loci"], sep="\t", index=False)
    female_only_df.to_csv(output_files["female_only_loci"], sep="\t", index=False)
    summary_df.to_csv(output_files["summary"], sep="\t", index=False)

    for output_file in output_files.values():
        print(f"Saved: {output_file}")


def main() -> None:
    phenotype = parse_phenotype()
    paths = get_paths(phenotype)

    create_output_dirs(paths)

    male_file = paths["loci_dir"] / f"{phenotype}_male_loci.tsv"
    female_file = paths["loci_dir"] / f"{phenotype}_female_loci.tsv"

    male_df = load_loci(male_file, phenotype, "male")
    female_df = load_loci(female_file, phenotype, "female")

    shared_df, male_only_df, female_only_df = compare_loci(male_df, female_df)

    summary_df = build_summary(
        male_df=male_df,
        female_df=female_df,
        shared_df=shared_df,
        male_only_df=male_only_df,
        female_only_df=female_only_df,
    )

    save_results(
        phenotype=phenotype,
        paths=paths,
        shared_df=shared_df,
        male_only_df=male_only_df,
        female_only_df=female_only_df,
        summary_df=summary_df,
    )

    print("=" * 80)
    print(f"FUMA-like locus comparison completed successfully for phenotype: {phenotype}")
    print("=" * 80)
    print(summary_df)


if __name__ == "__main__":
    main()
