import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent

REQUIRED_FILES = [
    "3_tools_results/fuma/male/magma.genes.out",
    "3_tools_results/fuma/female/magma.genes.out",

    "3_tools_results/gprofiler/male_pathways.tsv",
    "3_tools_results/gprofiler/female_pathways.tsv",
    "3_tools_results/gprofiler/shared_pathways.tsv",

    "3_tools_results/cytoscape/female_network_clean.tsv",
    "3_tools_results/cytoscape/male_network_clean.tsv",
    "3_tools_results/cytoscape/shared_network_clean.tsv",
]

STEPS = [
    "src/comparison/extract_compare_genes.py",
    "src/comparison/compare_pathways.py",
    "src/visualization/prepare_cytoscape_network.py",
    "src/comparison/real_hub_genes_compute.py",
    "src/visualization/gwas_plots.py",
    "src/visualization/plot_network_hubs.py",
    "src/visualization/plot_pathway.py",
]


def check_required_files() -> None:
    missing_files = []

    for file_path in REQUIRED_FILES:
        full_path = PROJECT_ROOT / file_path
        if not full_path.exists():
            missing_files.append(file_path)

    if missing_files:
        print("Missing required external tool result files:")
        for file in missing_files:
            print(f"- {file}")

        raise FileNotFoundError(
            "Some required files are missing. "
            "Please complete FUMA, g:Profiler, STRING/Cytoscape steps first."
        )


def run_step(script_path: str) -> None:
    full_path = PROJECT_ROOT / script_path

    if not full_path.exists():
        raise FileNotFoundError(f"Script not found: {full_path}")

    print("=" * 80)
    print(f"Running: {script_path}")
    print("=" * 80)

    result = subprocess.run(
        [sys.executable, str(full_path)],
        cwd=PROJECT_ROOT
    )

    if result.returncode != 0:
        raise RuntimeError(f"Pipeline stopped at: {script_path}")


def main() -> None:
    print("Checking required external tool results...")
    check_required_files()

    print("Starting post-tools pipeline...")

    for step in STEPS:
        run_step(step)

    print("=" * 80)
    print("Post-tools pipeline completed successfully.")
    print("Final tables saved in: 4_results/tables/")
    print("Final figures saved in: 4_results/figures/")
    print("=" * 80)


if __name__ == "__main__":
    main()
