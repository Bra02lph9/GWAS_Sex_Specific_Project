from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

RAW_DIR = PROJECT_ROOT / "1_data" / "raw"

FILTERED_ROOT = PROJECT_ROOT / "1_data" / "filtered"
GENE_LISTS_ROOT = PROJECT_ROOT / "1_data" / "gene_lists"

FUMA_ROOT = PROJECT_ROOT / "3_tools_results" / "fuma"
GPROFILER_ROOT = PROJECT_ROOT / "3_tools_results" / "gprofiler"
STRING_ROOT = PROJECT_ROOT / "3_tools_results" / "string"
CYTOSCAPE_ROOT = PROJECT_ROOT / "3_tools_results" / "cytoscape"

RESULTS_ROOT = PROJECT_ROOT / "4_results"
REPORTS_ROOT = PROJECT_ROOT / "reports"
LOCI_ROOT = PROJECT_ROOT / "1_data" / "loci"

GWAS_P_THRESHOLD = 5e-8
MAGMA_P_THRESHOLD = 0.05
CHUNK_SIZE = 500_000


def get_paths(phenotype: str) -> dict:
    """
    Return all project paths for one phenotype/disease.

    Expected raw input:
        1_data/raw/{PHENOTYPE}/male.tsv
        1_data/raw/{PHENOTYPE}/female.tsv
    """
    phenotype = phenotype.upper().strip()

    return {
        "phenotype": phenotype,

        "raw_dir": RAW_DIR / phenotype,
        "male_raw": RAW_DIR / phenotype / "male.tsv",
        "female_raw": RAW_DIR / phenotype / "female.tsv",

        "filtered_dir": FILTERED_ROOT / phenotype,
        "gene_lists_dir": GENE_LISTS_ROOT / phenotype,

        "fuma_dir": FUMA_ROOT / phenotype,
        "fuma_male_dir": FUMA_ROOT / phenotype / "male",
        "fuma_female_dir": FUMA_ROOT / phenotype / "female",

        "gprofiler_dir": GPROFILER_ROOT / phenotype,
        "string_dir": STRING_ROOT / phenotype,
        "cytoscape_dir": CYTOSCAPE_ROOT / phenotype,

        "results_dir": RESULTS_ROOT / phenotype,
        "tables_dir": RESULTS_ROOT / phenotype / "tables",
        "figures_dir": RESULTS_ROOT / phenotype / "figures",

        "report_dir": REPORTS_ROOT / phenotype,
        "loci_dir": LOCI_ROOT / phenotype,
    }


def create_output_dirs(paths: dict) -> None:
    """Create all phenotype-specific output directories."""
    output_keys = [
        "filtered_dir",
        "gene_lists_dir",
        "fuma_dir",
        "fuma_male_dir",
        "fuma_female_dir",
        "gprofiler_dir",
        "string_dir",
        "cytoscape_dir",
        "results_dir",
        "tables_dir",
        "figures_dir",
        "report_dir",
        "loci_dir",
    ]

    for key in output_keys:
        paths[key].mkdir(parents=True, exist_ok=True)
