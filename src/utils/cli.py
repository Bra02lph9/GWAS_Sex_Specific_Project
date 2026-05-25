import argparse


def parse_phenotype() -> str:

    parser = argparse.ArgumentParser(
        description="Run sex-specific GWAS pipeline for one phenotype."
    )

    parser.add_argument(
        "--phenotype",
        type=str,
        required=True,
        help="Phenotype/disease name (e.g. CAD, MI, AP)"
    )

    phenotype = parser.parse_args().phenotype.upper().strip()

    if not phenotype:
        raise ValueError("Phenotype name cannot be empty.")

    return phenotype
