from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[2]

RAW_DIR = PROJECT_ROOT / "1_data" / "raw"
FILTERED_DIR = PROJECT_ROOT / "1_data" / "filtered"
TABLES_DIR = PROJECT_ROOT / "4_results" / "tables"
FIGURES_DIR = PROJECT_ROOT / "4_results" / "figures"

PHENOTYPE = "CAD"


def count_rows_csv(file_path: Path) -> int:
    return sum(1 for _ in open(file_path)) - 1


def count_rows_tsv(file_path: Path) -> int:
    return sum(1 for _ in open(file_path)) - 1


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    male_raw = RAW_DIR / f"{PHENOTYPE}_male.tsv"
    female_raw = RAW_DIR / f"{PHENOTYPE}_female.tsv"

    male_filtered = FILTERED_DIR / f"{PHENOTYPE}_male_significant_snps.tsv"
    female_filtered = FILTERED_DIR / f"{PHENOTYPE}_female_significant_snps.tsv"

    shared_snps = TABLES_DIR / f"{PHENOTYPE}_shared_snps.tsv"
    male_only = TABLES_DIR / f"{PHENOTYPE}_male_only_snps.tsv"
    female_only = TABLES_DIR / f"{PHENOTYPE}_female_only_snps.tsv"

    values = {
        "Male Raw SNPs": count_rows_csv(male_raw),
        "Female Raw SNPs": count_rows_csv(female_raw),

        "Male Significant SNPs": count_rows_tsv(male_filtered),
        "Female Significant SNPs": count_rows_tsv(female_filtered),

        "Shared SNPs": count_rows_tsv(shared_snps),
        "Male-specific SNPs": count_rows_tsv(male_only),
        "Female-specific SNPs": count_rows_tsv(female_only),
    }

    labels = list(values.keys())
    counts = list(values.values())

    plt.figure(figsize=(12, 6))

    bars = plt.bar(labels, counts)

    plt.ylabel("Number of SNPs")
    plt.title(f"{PHENOTYPE} SNP Filtering and Comparison Summary")

    plt.xticks(rotation=30, ha="right")

    for bar, count in zip(bars, counts):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{count:,}",
            ha="center",
            va="bottom",
            fontsize=9
        )

    plt.tight_layout()

    output_file = FIGURES_DIR / f"{PHENOTYPE}_snp_filtering_summary.png"

    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Saved figure: {output_file}")


if __name__ == "__main__":
    main()
