from pathlib import Path
from typing import Any, Optional

import os
import time

import pandas as pd
import requests

from src.utils.cli import parse_phenotype
from src.utils.config import get_paths, create_output_dirs


POPULATION = "EUR"
GENOME_BUILD = "grch38"

R2_STRONG = 0.80
R2_MODERATE = 0.60
REQUEST_SLEEP_SECONDS = 1.2
REQUEST_TIMEOUT_SECONDS = 30

SEXES = ("male", "female")

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
]

LOCI_REQUIRED_COLUMNS = [
    "GenomicLocus",
    "chr",
    "start",
    "end",
    "IndSigSNPs",
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
    missing_columns = [col for col in required_columns if col not in df.columns]

    if missing_columns:
        raise ValueError(
            f"Missing required columns in {file_name}: {missing_columns}"
        )


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
    ]

    for column in text_columns:
        df[column] = df[column].astype(str).str.strip()

    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    df["nes"] = pd.to_numeric(df["nes"], errors="coerce")

    df = df.dropna(subset=["gene", "snp", "variant_id", "p_value", "nes"])

    return df.reset_index(drop=True)


def load_loci(file_path: Path, phenotype: str, sex: str) -> pd.DataFrame:
    df = read_tsv(file_path)
    validate_required_columns(df, LOCI_REQUIRED_COLUMNS, file_path.name)

    df = df[LOCI_REQUIRED_COLUMNS].copy()

    df = df.rename(
        columns={
            "GenomicLocus": "locus_id",
            "chr": "chromosome",
            "IndSigSNPs": "snp_list",
        }
    )

    numeric_columns = ["locus_id", "chromosome", "start", "end"]

    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    df = df.dropna(subset=numeric_columns)

    df["locus_id"] = df["locus_id"].astype(int)
    df["chromosome"] = df["chromosome"].astype(int)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    df["locus_id"] = df["locus_id"].apply(
        lambda x: f"{phenotype}_{sex}_LOCUS_{x:04d}"
    )

    df["snp_list"] = df["snp_list"].astype(str).str.strip()

    return df[["locus_id", "chromosome", "start", "end", "snp_list"]]


def parse_locus_ids(value: Any) -> list[str]:
    if pd.isna(value) or str(value).strip() == "":
        return []

    return [item.strip() for item in str(value).split(";") if item.strip()]


def get_locus_snps(loci_df: pd.DataFrame, locus_id: str) -> list[str]:
    matched_locus = loci_df[loci_df["locus_id"] == locus_id]

    if matched_locus.empty:
        return []

    snp_list = str(matched_locus["snp_list"].iloc[0])

    return [snp.strip() for snp in snp_list.split(";") if snp.strip()]


def query_ldpair(var1: str, var2: str, token: str) -> dict[str, Any]:
    url = "https://ldlink.nih.gov/LDlinkRest/ldpair"

    params = {
        "var1": var1,
        "var2": var2,
        "pop": POPULATION,
        "genome_build": GENOME_BUILD,
        "json_out": "true",
        "token": token,
    }

    try:
        response = requests.get(
            url,
            params=params,
            timeout=REQUEST_TIMEOUT_SECONDS,
        )
    except requests.RequestException as error:
        return {
            "ld_status": "request_error",
            "r2": None,
            "d_prime": None,
            "raw_response": str(error)[:300],
        }

    if response.status_code != 200:
        return {
            "ld_status": "api_error",
            "r2": None,
            "d_prime": None,
            "raw_response": response.text[:300],
        }

    try:
        data = response.json()
    except ValueError:
        return {
            "ld_status": "parse_error",
            "r2": None,
            "d_prime": None,
            "raw_response": response.text[:300],
        }

    if "error" in data:
        return {
            "ld_status": "ldlink_error",
            "r2": None,
            "d_prime": None,
            "raw_response": str(data["error"])[:300],
        }

    statistics = data.get("statistics", {})

    r2 = parse_float(statistics.get("r2") or statistics.get("R2"))
    d_prime = parse_float(
        statistics.get("d_prime")
        or statistics.get("D_prime")
        or statistics.get("D'")
    )

    return {
        "ld_status": "ok",
        "r2": r2,
        "d_prime": d_prime,
        "raw_response": "",
    }


def parse_float(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def classify_ld(r2: Optional[float]) -> str:
    if r2 is None:
        return "no_ld_result"

    if r2 >= R2_STRONG:
        return "strong_ld"

    if r2 >= R2_MODERATE:
        return "moderate_ld"

    return "weak_or_no_ld"


def build_pairs_for_eqtl(
    row: pd.Series,
    male_loci: pd.DataFrame,
    female_loci: pd.DataFrame,
) -> list[dict[str, Any]]:
    pairs: list[dict[str, Any]] = []

    eqtl_snp = row["snp"]

    sources = [
        ("male", row.get("matched_male_loci", ""), male_loci),
        ("female", row.get("matched_female_loci", ""), female_loci),
    ]

    for sex, locus_ids_raw, loci_df in sources:
        locus_ids = parse_locus_ids(locus_ids_raw)

        for locus_id in locus_ids:
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
    all_pairs: list[dict[str, Any]] = []

    for _, row in eqtl_df.iterrows():
        all_pairs.extend(
            build_pairs_for_eqtl(
                row=row,
                male_loci=male_loci,
                female_loci=female_loci,
            )
        )

    if not all_pairs:
        return pd.DataFrame()

    return (
        pd.DataFrame(all_pairs)
        .drop_duplicates(subset=["eqtl_snp", "cad_snp", "sex", "locus_id"])
        .reset_index(drop=True)
    )


def run_ld_validation(
    pairs_df: pd.DataFrame,
    token: str,
) -> pd.DataFrame:
    results: list[dict[str, Any]] = []

    total_pairs = len(pairs_df)

    for index, row in pairs_df.iterrows():
        print(
            f"[{index + 1}/{total_pairs}] LDpair: "
            f"{row['eqtl_snp']} vs {row['cad_snp']} "
            f"({row['sex']} | {row['locus_id']})"
        )

        ld_result = query_ldpair(
            var1=row["eqtl_snp"],
            var2=row["cad_snp"],
            token=token,
        )

        record = row.to_dict()
        record.update(ld_result)
        record["ld_class"] = classify_ld(ld_result["r2"])

        results.append(record)

        time.sleep(REQUEST_SLEEP_SECONDS)

    result_df = pd.DataFrame(results)

    return result_df.sort_values(
        by=["gene", "sex", "locus_id", "r2"],
        ascending=[True, True, True, False],
        na_position="last",
    ).reset_index(drop=True)


def save_outputs(
    result_df: pd.DataFrame,
    phenotype: str,
    paths: dict[str, Path],
) -> tuple[Path, Path]:
    output_file = paths["tables_dir"] / f"{phenotype}_eqtl_ld_validation.tsv"
    summary_file = paths["tables_dir"] / f"{phenotype}_eqtl_ld_validation_summary.tsv"

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
            'setx LDLINK_TOKEN "YOUR_TOKEN_HERE"\n'
            "Then close and reopen your terminal."
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
    print(f"Running LD validation for phenotype: {phenotype}")
    print(f"Validated eQTL file: {validated_eqtl_file}")
    print("=" * 80)

    eqtl_df = load_validated_eqtls(validated_eqtl_file)

    eqtl_df = eqtl_df[
        eqtl_df["validation_level"].isin(VALIDATION_LEVELS_FOR_LD)
    ].copy()

    male_loci = load_loci(male_loci_file, phenotype, "male")
    female_loci = load_loci(female_loci_file, phenotype, "female")

    pairs_df = build_all_ld_pairs(
        eqtl_df=eqtl_df,
        male_loci=male_loci,
        female_loci=female_loci,
    )

    output_file = paths["tables_dir"] / f"{phenotype}_eqtl_ld_validation.tsv"

    if pairs_df.empty:
        pd.DataFrame().to_csv(output_file, sep="\t", index=False)
        print("No locus-overlapping eQTLs available for LD validation.")
        print(f"Saved empty output file: {output_file}")
        return

    result_df = run_ld_validation(
        pairs_df=pairs_df,
        token=token,
    )

    output_file, summary_file = save_outputs(
        result_df=result_df,
        phenotype=phenotype,
        paths=paths,
    )

    summary_df = read_tsv(summary_file)

    print("-" * 80)
    print(f"Saved LD validation table  : {output_file}")
    print(f"Saved LD validation summary: {summary_file}")
    print(summary_df)


if __name__ == "__main__":
    main()
