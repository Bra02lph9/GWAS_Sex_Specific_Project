import os
import subprocess
import sys
from pathlib import Path

from src.utils.cli import parse_phenotype
from src.utils.config import get_paths, create_output_dirs


PROJECT_ROOT = Path(__file__).resolve().parent

STEPS = [
    "src/preprocessing/filter_large_gwas.py",
    "src/comparison/compare_significant_snps.py",
    "src/comparison/compare_snps.py",
    "src/locus/define_loci.py",
    "src/locus/compare_loci.py",
    "src/preprocessing/prepare_fuma_input.py",
]


def run_step(script_path: str, phenotype: str) -> None:
    full_path = PROJECT_ROOT / script_path

    if not full_path.exists():
        raise FileNotFoundError(f"Script not found: {full_path}")

    print("=" * 80)
    print(f"Running: {script_path}")
    print(f"Phenotype: {phenotype}")
    print("=" * 80)

    env = os.environ.copy()
    env["PYTHONPATH"] = str(PROJECT_ROOT)

    result = subprocess.run(
        [
            sys.executable,
            str(full_path),
            "--phenotype",
            phenotype,
        ],
        cwd=PROJECT_ROOT,
        env=env,
    )

    if result.returncode != 0:
        raise RuntimeError(f"Pipeline stopped at: {script_path}")


def main() -> None:
    phenotype = parse_phenotype()
    paths = get_paths(phenotype)

    create_output_dirs(paths)

    print("=" * 80)
    print(f"Starting pre-FUMA pipeline for phenotype: {phenotype}")
    print("=" * 80)

    for step in STEPS:
        run_step(step, phenotype)

    print("=" * 80)
    print(f"Pre-FUMA pipeline completed successfully for: {phenotype}")
    print(f"Next step: upload FUMA input files from:")
    print(paths["fuma_dir"])
    print("=" * 80)


if __name__ == "__main__":
    main()
