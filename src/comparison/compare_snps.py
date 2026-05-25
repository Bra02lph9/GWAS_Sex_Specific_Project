from pathlib import Path
import matplotlib.pyplot as plt

from src.utils.cli import parse_phenotype
from src.utils.config import get_paths, create_output_dirs


def count_rows(file_path: Path) -> int:
    if not file_path.exists():
        raise FileNotFoundError(f"Missing file: {file_path}")

    with open(file_path, "r", encoding="utf-8") as file:
        return max(sum(1 for _ in file) - 1, 0)


def main() -> None:
    phenotype = parse_phenotype()
    paths = get_paths(phenotype)
    create_output_dirs(paths)

    male_raw = paths["male_raw"]
    female_raw = paths["female_raw"]

    male_filtered = paths["filtered_dir"] / f"{phenotype}_male_significant_snps.tsv"
    female_filtered = paths["filtered_dir"] / f"{phenotype}_female_significant_snps.tsv"

    shared_snps = paths["tables_dir"] / f"{phenotype}_shared_snps.tsv"
    male_only = paths["tables_dir"] / f"{phenotype}_male_only_snps.tsv"
    female_only = paths["tables_dir"] / f"{phenotype}_female_only_snps.tsv"

    values = {
        "Male Raw SNPs": count_rows(male_raw),
        "Female Raw SNPs": count_rows(female_raw),
        "Male Significant SNPs": count_rows(male_filtered),
        "Female Significant SNPs": count_rows(female_filtered),
        "Shared SNPs": count_rows(shared_snps),
        "Male-specific SNPs": count_rows(male_only),
        "Female-specific SNPs": count_rows(female_only),
    }

    labels = list(values.keys())
    counts = list(values.values())

    plt.figure(figsize=(12, 6))
    bars = plt.bar(labels, counts)

    plt.ylabel("Number of SNPs")
    plt.title(f"{phenotype} SNP Filtering and Comparison Summary")
    plt.xticks(rotation=30, ha="right")

    for bar, count in zip(bars, counts):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{count:,}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    plt.tight_layout()

    output_file = paths["figures_dir"] / f"{phenotype}_snp_filtering_summary.png"
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Saved figure: {output_file}")


if __name__ == "__main__":
    main()
