from pathlib import Path
from typing import Any

import pandas as pd

from src.utils.cli import parse_phenotype
from src.utils.config import get_paths, create_output_dirs


MERGE_DISTANCE = 250_000
SEXES = ("male", "female")

REQUIRED_COLUMNS = ["SNP", "chromosome", "position", "p_value"]

OUTPUT_COLUMNS = [
    "GenomicLocus",
    "uniqID",
    "rsID",
    "chr",
    "pos",
    "p",
    "start",
    "end",
    "nSNPs",
    "nGWASSNPs",
    "nIndSigSNPs",
    "IndSigSNPs",
    "nLeadSNPs",
    "LeadSNPs",
]


def load_significant_snps(file_path: Path) -> pd.DataFrame:
    if not file_path.exists():
        raise FileNotFoundError(f"Missing SNP file: {file_path}")

    df = pd.read_csv(file_path, sep="\t", low_memory=False)

    missing_columns = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns in {file_path.name}: {missing_columns}")

    df = df[REQUIRED_COLUMNS].copy()

    df["SNP"] = df["SNP"].astype(str).str.strip()
    df["chromosome"] = pd.to_numeric(df["chromosome"], errors="coerce")
    df["position"] = pd.to_numeric(df["position"], errors="coerce")
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")

    df = df.dropna(subset=REQUIRED_COLUMNS)

    df = df[
        (df["SNP"] != "")
        & (df["SNP"].str.lower() != "nan")
        & (df["chromosome"].between(1, 22))
        & (df["position"] > 0)
        & (df["p_value"].between(0, 1, inclusive="right"))
    ].copy()

    df["chromosome"] = df["chromosome"].astype(int)
    df["position"] = df["position"].astype(int)

    return df.sort_values(["chromosome", "position"]).reset_index(drop=True)


def build_fuma_like_locus(
    locus_number: int,
    locus_snps: list[pd.Series],
) -> dict[str, Any]:
    locus_df = pd.DataFrame(locus_snps).sort_values("p_value")

    lead = locus_df.iloc[0]

    chromosome = int(lead["chromosome"])
    lead_position = int(lead["position"])
    lead_snp = str(lead["SNP"])
    lead_p = lead["p_value"]

    start = int(locus_df["position"].min())
    end = int(locus_df["position"].max())

    all_snps = locus_df["SNP"].astype(str).tolist()

    return {
        "GenomicLocus": locus_number,
        "uniqID": f"{chromosome}:{lead_position}:NA:NA",
        "rsID": lead_snp,
        "chr": chromosome,
        "pos": lead_position,
        "p": lead_p,
        "start": start,
        "end": end,
        "nSNPs": len(locus_df),
        "nGWASSNPs": len(locus_df),
        "nIndSigSNPs": len(locus_df),
        "IndSigSNPs": ";".join(all_snps),
        "nLeadSNPs": 1,
        "LeadSNPs": lead_snp,
    }


def define_fuma_like_loci(snps_df: pd.DataFrame) -> pd.DataFrame:
    if snps_df.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)

    loci = []
    locus_number = 1

    for chromosome, chrom_df in snps_df.groupby("chromosome"):
        chrom_df = chrom_df.sort_values("position").reset_index(drop=True)

        current_locus_snps = []
        current_locus_end = None

        for _, row in chrom_df.iterrows():
            position = int(row["position"])

            if current_locus_end is None:
                current_locus_snps = [row]
                current_locus_end = position
                continue

            if position - current_locus_end <= MERGE_DISTANCE:
                current_locus_snps.append(row)
                current_locus_end = position
            else:
                loci.append(build_fuma_like_locus(locus_number, current_locus_snps))
                locus_number += 1

                current_locus_snps = [row]
                current_locus_end = position

        if current_locus_snps:
            loci.append(build_fuma_like_locus(locus_number, current_locus_snps))
            locus_number += 1

    return pd.DataFrame(loci, columns=OUTPUT_COLUMNS)


def process_sex(phenotype: str, sex: str, paths: dict[str, Path]) -> None:
    input_file = paths["filtered_dir"] / f"{phenotype}_{sex}_significant_snps.tsv"
    output_file = paths["loci_dir"] / f"{phenotype}_{sex}_fuma_like_loci.tsv"

    print("=" * 80)
    print(f"Defining FUMA-like loci for phenotype: {phenotype} | sex: {sex}")
    print(f"Input : {input_file}")
    print(f"Output: {output_file}")
    print("=" * 80)

    snps_df = load_significant_snps(input_file)
    loci_df = define_fuma_like_loci(snps_df)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    loci_df.to_csv(output_file, sep="\t", index=False)

    print(f"Significant SNPs loaded: {len(snps_df):,}")
    print(f"FUMA-like loci defined : {len(loci_df):,}")
    print(f"Saved successfully to  : {output_file}")


def main() -> None:
    phenotype = parse_phenotype()
    paths = get_paths(phenotype)

    create_output_dirs(paths)

    for sex in SEXES:
        process_sex(phenotype, sex, paths)

    print("-" * 80)
    print(f"FUMA-like locus definition completed successfully for phenotype: {phenotype}")


if __name__ == "__main__":
    main()