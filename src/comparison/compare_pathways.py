import pandas as pd

from src.utils.cli import parse_phenotype
from src.utils.config import get_paths, create_output_dirs


COLUMNS_TO_KEEP = [
    "term_name",
    "source",
    "adjusted_p_value",
    "intersection_size",
]


def load_pathway_file(file_path):
    if not file_path.exists():
        raise FileNotFoundError(f"Pathway file not found: {file_path}")

    df = pd.read_csv(file_path, sep="\t")

    missing_columns = [col for col in COLUMNS_TO_KEEP if col not in df.columns]
    if missing_columns:
        raise ValueError(
            f"Missing columns in {file_path}: {missing_columns}"
        )

    return df[COLUMNS_TO_KEEP].copy()


def main():
    phenotype = parse_phenotype()
    paths = get_paths(phenotype)
    create_output_dirs(paths)

    female_file = paths["gprofiler_dir"] / "female_pathways.tsv"
    male_file = paths["gprofiler_dir"] / "male_pathways.tsv"

    female_df = load_pathway_file(female_file)
    male_df = load_pathway_file(male_file)

    female_pathways = set(female_df["term_name"])
    male_pathways = set(male_df["term_name"])

    shared_pathways = female_pathways.intersection(male_pathways)
    female_specific = female_pathways - male_pathways
    male_specific = male_pathways - female_pathways

    shared_df = pd.concat([
        female_df[female_df["term_name"].isin(shared_pathways)],
        male_df[male_df["term_name"].isin(shared_pathways)],
    ]).drop_duplicates()

    female_specific_df = female_df[
        female_df["term_name"].isin(female_specific)
    ].sort_values("adjusted_p_value")

    male_specific_df = male_df[
        male_df["term_name"].isin(male_specific)
    ].sort_values("adjusted_p_value")

    shared_output = paths["tables_dir"] / "shared_pathways.csv"
    female_output = paths["tables_dir"] / "female_specific_pathways.csv"
    male_output = paths["tables_dir"] / "male_specific_pathways.csv"

    shared_df.to_csv(shared_output, index=False)
    female_specific_df.to_csv(female_output, index=False)
    male_specific_df.to_csv(male_output, index=False)

    print("===== PATHWAY COMPARISON =====")
    print(f"Phenotype: {phenotype}")
    print(f"Female pathways: {len(female_pathways)}")
    print(f"Male pathways: {len(male_pathways)}")
    print(f"Shared pathways: {len(shared_pathways)}")
    print(f"Female-specific pathways: {len(female_specific)}")
    print(f"Male-specific pathways: {len(male_specific)}")

    print("\nFiles generated:")
    print(f"- {shared_output}")
    print(f"- {female_output}")
    print(f"- {male_output}")


if __name__ == "__main__":
    main()
