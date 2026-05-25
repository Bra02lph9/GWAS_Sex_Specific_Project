import subprocess
import sys
from pathlib import Path

from src.utils.cli import parse_phenotype
from src.utils.config import get_paths, create_output_dirs


PROJECT_ROOT = Path(__file__).resolve().parent

STEPS = [
    "src/comparison/extract_compare_genes.py",
    "src/comparison/compare_pathways.py",
    "src/visualization/prepare_cytoscape_network.py",
    "src/comparison/real_hub_genes_compute.py",
    "src/visualization/gwas_plots.py",
    "src/visualization/plot_network_hubs.py",
    "src/visualization/plot_pathway.py",
]


def check_required_files(paths: dict) -> None:
    required_files = [
        paths["fuma_male_dir"] / "magma.genes.out",
        paths["fuma_female_dir"] / "magma.genes.out",

        paths["gprofiler_dir"] / "male_pathways.tsv",
        paths["gprofiler_dir"] / "female_pathways.tsv",
        paths["gprofiler_dir"] / "shared_pathways.tsv",

        paths["cytoscape_dir"] / "female_network_clean.tsv",
        paths["cytoscape_dir"] / "male_network_clean.tsv",
        paths["cytoscape_dir"] / "shared_network_clean.tsv",
    ]

    missing_files = [
        file_path for file_path in required_files
        if not file_path.exists()
    ]

    if missing_files:
        print("=" * 80)
        print("Missing required external tool result files:")
        print("=" * 80)

        for file_path in missing_files:
            print(f"- {file_path}")

        raise FileNotFoundError(
            "\nSome required files are missing.\n"
            "Please complete FUMA, g:Profiler, STRING, "
            "and Cytoscape steps before running this pipeline."
        )


def run_step(script_path: str, phenotype: str) -> None:
    full_path = PROJECT_ROOT / script_path

    if not full_path.exists():
        raise FileNotFoundError(f"Script not found: {full_path}")

    print("=" * 80)
    print(f"Running: {script_path}")
    print(f"Phenotype: {phenotype}")
    print("=" * 80)

    result = subprocess.run(
        [
            sys.executable,
            str(full_path),
            "--phenotype",
            phenotype,
        ],
        cwd=PROJECT_ROOT,
    )

    if result.returncode != 0:
        raise RuntimeError(f"Pipeline stopped at: {script_path}")


def main() -> None:
    phenotype = parse_phenotype()
    paths = get_paths(phenotype)

    create_output_dirs(paths)

    print("=" * 80)
    print(f"Checking required external tool results for: {phenotype}")
    print("=" * 80)

    check_required_files(paths)

    print("=" * 80)
    print(f"Starting post-tools pipeline for: {phenotype}")
    print("=" * 80)

    for step in STEPS:
        run_step(step, phenotype)

    print("=" * 80)
    print(f"Post-tools pipeline completed successfully for: {phenotype}")
    print(f"Final tables saved in: {paths['tables_dir']}")
    print(f"Final figures saved in: {paths['figures_dir']}")
    print("=" * 80)


if __name__ == "__main__":
    main()
