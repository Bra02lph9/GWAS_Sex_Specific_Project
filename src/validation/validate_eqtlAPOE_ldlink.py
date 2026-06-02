from pathlib import Path
from typing import Any, Optional

import os
import time
import pandas as pd
import requests
import re

from src.utils.cli import parse_phenotype
from src.utils.config import get_paths, create_output_dirs


POPULATION = "ALL"
GENOME_BUILD = "grch38"

TARGET_GENES = {"APOE"}

R2_STRONG = 0.80
R2_MODERATE = 0.60
REQUEST_SLEEP_SECONDS = 1.2
REQUEST_TIMEOUT_SECONDS = 30

VALIDATION_LEVELS_FOR_LD = {
    "good_locus_overlap",
    "strong_exact_snp_match",
}

EQTL_REQUIRED_COLUMNS = [
    "gene",
    "group",
    "variant_id",
    "snp",
    "p_value",
    "nes",
    "tissue",
    "validation_level",
    "matched_male_loci",
    "matched_female_loci",
]

LOCI_REQUIRED_COLUMNS = [
    "GenomicLocus",
    "chr",
    "start",
    "end",
    "IndSigSNPs",
    "LeadSNPs",
]


def read_tsv(file_path: Path) -> pd.DataFrame:
    if not file_path.exists():
        raise FileNotFoundError(f"Missing file: {file_path}")

    return pd.read_csv(file_path, sep="\t", low_memory=False)


def validate_required_columns(
    df: pd.DataFrame,
    required_columns: list[str],
    file_name: str,
) -> None:
    missing = [col for col in required_columns if col not in df.columns]

    if missing:
        raise ValueError(f"Missing required columns in {file_name}: {missing}")


def load_validated_eqtls(file_path: Path) -> pd.DataFrame:
    df = read_tsv(file_path)
    validate_required_columns(df, EQTL_REQUIRED_COLUMNS, file_path.name)

    df = df[EQTL_REQUIRED_COLUMNS].copy()

    text_columns = [
        "gene",
        "group",
        "variant_id",
        "snp",
        "tissue",
        "validation_level",
        "matched_male_loci",
        "matched_female_loci",
    ]

    for col in text_columns:
        df[col] = df[col].astype(str).str.strip()

    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    df["nes"] = pd.to_numeric(df["nes"], errors="coerce")

    df = df.dropna(subset=["gene", "snp", "variant_id", "p_value", "nes"])

    df = df[
        df["gene"].isin(TARGET_GENES)
        & df["validation_level"].isin(VALIDATION_LEVELS_FOR_LD)
    ].copy()

    return df.reset_index(drop=True)


def load_loci(file_path: Path, phenotype: str, sex: str) -> pd.DataFrame:
    df = read_tsv(file_path)
    validate_required_columns(df, LOCI_REQUIRED_COLUMNS, file_path.name)

    df = df[LOCI_REQUIRED_COLUMNS].copy()

    df = df.rename(
        columns={
            "GenomicLocus": "locus_id",
            "chr": "chromosome",
        }
    )

    for col in ["locus_id", "chromosome", "start", "end"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["locus_id", "chromosome", "start", "end"])

    df["locus_id"] = df["locus_id"].astype(int)
    df["chromosome"] = df["chromosome"].astype(int)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    df["locus_id"] = df["locus_id"].apply(
        lambda x: f"{phenotype}_{sex}_LOCUS_{x:04d}"
    )

    df["IndSigSNPs"] = df["IndSigSNPs"].astype(str).str.strip()
    df["LeadSNPs"] = df["LeadSNPs"].astype(str).str.strip()

    df["snp_list"] = df.apply(
        lambda row: merge_snp_lists(row["IndSigSNPs"], row["LeadSNPs"]),
        axis=1,
    )

    return df[["locus_id", "chromosome", "start", "end", "snp_list"]]


def merge_snp_lists(ind_sig_snps: str, lead_snps: str) -> str:
    snps = []

    for value in [ind_sig_snps, lead_snps]:
        if value.lower() in {"nan", "none", ""}:
            continue

        snps.extend([s.strip() for s in value.split(";") if s.strip()])

    return ";".join(sorted(set(snps)))


def parse_locus_ids(value: Any) -> list[str]:
    if pd.isna(value):
        return []

    value = str(value).strip()

    if value == "" or value.lower() in {"nan", "none"}:
        return []

    return [x.strip() for x in value.split(";") if x.strip()]


def get_locus_snps(loci_df: pd.DataFrame, locus_id: str) -> list[str]:
    row = loci_df[loci_df["locus_id"] == locus_id]

    if row.empty:
        return []

    snp_list = str(row["snp_list"].iloc[0])

    return [snp.strip() for snp in snp_list.split(";") if snp.strip()]


def parse_float(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def query_ldpair(var1: str, var2: str, token: str) -> dict:
    url = "https://ldlink.nih.gov/LDlinkRest/ldpair"

    params = {
        "var1": var1,
        "var2": var2,
        "pop": POPULATION,
        "genome_build": GENOME_BUILD,
        "json_out": "false",
        "token": token,
    }

    response = requests.get(url, params=params, timeout=30)
    text = response.text

    if response.status_code != 200:
        return {
            "ld_status": "api_error",
            "r2": None,
            "d_prime": None,
            "raw_response": text[:500],
        }

    if "error" in text.lower() or "warning" in text.lower():
        return {
            "ld_status": "ldlink_error_or_warning",
            "r2": None,
            "d_prime": None,
            "raw_response": text[:500],
        }

    r2_match = re.search(r"R\^2:\s*([0-9.]+)", text)
    d_match = re.search(r"D':\s*([0-9.]+)", text)
    d_prime = float(d_match.group(1)) if d_match else None
    r2 = float(r2_match.group(1)) if r2_match else None

    return {
        "ld_status": "ok",
        "r2": r2,
        "d_prime": d_prime,
        "raw_response": text[:500],
    }


def classify_ld(r2: float | None, d_prime: float | None) -> str:

    if r2 is not None and r2 >= 0.8:
        return "strong_ld_r2"

    if d_prime is not None and d_prime >= 0.8:
        return "strong_ld_dprime"

    if r2 is not None and r2 >= 0.6:
        return "moderate_ld"

    if r2 is not None and r2 >= 0.3:
        return "weak_ld"

    if r2 is None and d_prime is None:
        return "no_ld_result"

    return "weak_or_no_ld"


def build_pairs_for_eqtl(
    row: pd.Series,
    male_loci: pd.DataFrame,
    female_loci: pd.DataFrame,
) -> list[dict[str, Any]]:
    pairs = []
    eqtl_snp = row["snp"]

    sources = [
        ("male", row.get("matched_male_loci", ""), male_loci),
        ("female", row.get("matched_female_loci", ""), female_loci),
    ]

    for sex, locus_ids_raw, loci_df in sources:
        for locus_id in parse_locus_ids(locus_ids_raw):
            locus_snps = get_locus_snps(loci_df, locus_id)

            for cad_snp in locus_snps:
                if cad_snp == eqtl_snp:
                    continue

                pairs.append(
                    {
                        "gene": row["gene"],
                        "group": row["group"],
                        "eqtl_snp": eqtl_snp,
                        "variant_id": row["variant_id"],
                        "tissue": row["tissue"],
                        "eqtl_p_value": row["p_value"],
                        "eqtl_nes": row["nes"],
                        "sex": sex,
                        "locus_id": locus_id,
                        "cad_snp": cad_snp,
                    }
                )

    return pairs


def build_all_ld_pairs(
    eqtl_df: pd.DataFrame,
    male_loci: pd.DataFrame,
    female_loci: pd.DataFrame,
) -> pd.DataFrame:
    pairs = []

    for _, row in eqtl_df.iterrows():
        pairs.extend(build_pairs_for_eqtl(row, male_loci, female_loci))

    if not pairs:
        return pd.DataFrame()

    return (
        pd.DataFrame(pairs)
        .drop_duplicates(subset=["eqtl_snp", "cad_snp", "sex", "locus_id"])
        .reset_index(drop=True)
    )


def run_ld_validation(pairs_df: pd.DataFrame, token: str) -> pd.DataFrame:
    results = []
    total = len(pairs_df)

    for i, row in pairs_df.iterrows():
        print(
            f"[{i + 1}/{total}] LDpair: "
            f"{row['eqtl_snp']} vs {row['cad_snp']} "
            f"({row['sex']} | {row['locus_id']})"
        )

        ld_result = query_ldpair(row["eqtl_snp"], row["cad_snp"], token)

        record = row.to_dict()
        record.update(ld_result)
        record["ld_class"] = classify_ld(ld_result["r2"], ld_result["d_prime"])

        results.append(record)
        time.sleep(REQUEST_SLEEP_SECONDS)

    return (
        pd.DataFrame(results)
        .sort_values(
            by=["gene", "sex", "locus_id", "r2"],
            ascending=[True, True, True, False],
            na_position="last",
        )
        .reset_index(drop=True)
    )


def save_outputs(
    result_df: pd.DataFrame,
    phenotype: str,
    paths: dict[str, Path],
) -> tuple[Path, Path]:
    output_file = paths["tables_dir"] / f"{phenotype}_APOE_eqtl_ld_validation.tsv"
    summary_file = paths["tables_dir"] / f"{phenotype}_APOE_eqtl_ld_validation_summary.tsv"

    result_df.to_csv(output_file, sep="\t", index=False)

    summary_df = result_df["ld_class"].value_counts().reset_index()
    summary_df.columns = ["ld_class", "count"]
    summary_df.to_csv(summary_file, sep="\t", index=False)

    return output_file, summary_file


def get_ldlink_token() -> str:
    token = os.getenv("LDLINK_TOKEN")

    if not token:
        raise EnvironmentError(
            "Missing LDlink token. Set it in PowerShell with:\n"
            '$env:LDLINK_TOKEN="YOUR_TOKEN_HERE"'
        )

    return token


def main() -> None:
    phenotype = parse_phenotype()
    paths = get_paths(phenotype)
    create_output_dirs(paths)

    token = get_ldlink_token()

    validated_eqtl_file = paths["tables_dir"] / f"{phenotype}_eqtl_validated.tsv"
    male_loci_file = paths["loci_dir"] / f"{phenotype}_male_loci.tsv"
    female_loci_file = paths["loci_dir"] / f"{phenotype}_female_loci.tsv"

    print("=" * 80)
    print(f"Running APOE-only LD validation for phenotype: {phenotype}")
    print(f"Validated eQTL file: {validated_eqtl_file}")
    print("=" * 80)

    eqtl_df = load_validated_eqtls(validated_eqtl_file)
    male_loci = load_loci(male_loci_file, phenotype, "male")
    female_loci = load_loci(female_loci_file, phenotype, "female")

    print(f"APOE validated eQTL rows kept: {len(eqtl_df)}")

    pairs_df = build_all_ld_pairs(eqtl_df, male_loci, female_loci)

    output_file = paths["tables_dir"] / f"{phenotype}_APOE_eqtl_ld_validation.tsv"

    if pairs_df.empty:
        pd.DataFrame().to_csv(output_file, sep="\t", index=False)
        print("No APOE LD pairs found.")
        print(f"Saved empty output file: {output_file}")
        return

    print(f"APOE LD pairs to test: {len(pairs_df)}")

    result_df = run_ld_validation(pairs_df, token)

    output_file, summary_file = save_outputs(result_df, phenotype, paths)

    print("-" * 80)
    print(f"Saved APOE LD validation table  : {output_file}")
    print(f"Saved APOE LD validation summary: {summary_file}")
    print(read_tsv(summary_file))


if __name__ == "__main__":
    main()
