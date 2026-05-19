import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent

STEPS = [
    "src/preprocessing/filter_large_gwas.py",
    "src/comparison/compare_significant_snps.py",
    "src/comparison/compare_snps.py",
    "src/preprocessing/prepare_fuma_input.py",
]


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
    print("Starting pre-FUMA pipeline...")

    for step in STEPS:
        run_step(step)

    print("=" * 80)
    print("Pre-FUMA pipeline completed successfully.")
    print("Next step: upload FUMA input files from 3_tools_results/fuma/")
    print("=" * 80)


if __name__ == "__main__":
    main()
