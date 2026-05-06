from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


PROJECT_ROOT = Path(__file__).resolve().parents[2]

RAW_DIR = PROJECT_ROOT / "1_data" / "raw"
FILTERED_DIR = PROJECT_ROOT / "1_data" / "filtered"
GENE_LISTS_DIR = PROJECT_ROOT / "1_data" / "gene_lists"
TABLES_DIR = PROJECT_ROOT / "4_results" / "tables"
FIGURES_DIR = PROJECT_ROOT / "4_results" / "figures"

PHENOTYPE = "CAD"
CHUNK_SIZE = 500_000
GWAS_THRESHOLD = 5e-8


PLOT_COLUMNS = ["CHR", "BP", "P", "BETA", "SNP"]


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = (
        df.columns
        .str.strip()
        .str.lower()
        .str.replace(" ", "_", regex=False)
        .str.replace("-", "_", regex=False)
    )

    rename_map = {
        "rsid": "SNP",
        "snp": "SNP",
        "chr": "CHR",
        "chromosome": "CHR",
        "chrom": "CHR",
        "bp": "BP",
        "position": "BP",
        "pos": "BP",
        "p": "P",
        "p_value": "P",
        "beta": "BETA",
    }

    return df.rename(columns=rename_map)


def load_gwas_for_plot(file_path: Path) -> pd.DataFrame:
    if not file_path.exists():
        raise FileNotFoundError(f"Missing file: {file_path}")

    chunks = []

    print(f"Loading GWAS file for plotting: {file_path.name}")

    for chunk in pd.read_csv(
        file_path,
        sep=",",
        chunksize=CHUNK_SIZE,
        low_memory=False,
    ):
        chunk = normalize_columns(chunk)

        missing_columns = [col for col in PLOT_COLUMNS if col not in chunk.columns]
        if missing_columns:
            raise ValueError(
                f"Missing required columns in {file_path.name}: {missing_columns}"
            )

        chunk = chunk[PLOT_COLUMNS].copy()

        chunk["SNP"] = chunk["SNP"].astype(str).str.strip()
        chunk["CHR"] = pd.to_numeric(chunk["CHR"], errors="coerce")
        chunk["BP"] = pd.to_numeric(chunk["BP"], errors="coerce")
        chunk["P"] = pd.to_numeric(chunk["P"], errors="coerce")
        chunk["BETA"] = pd.to_numeric(chunk["BETA"], errors="coerce")

        chunk = chunk.dropna(subset=["CHR", "BP", "P", "BETA"])

        chunk = chunk[
            (chunk["P"] > 0)
            & (chunk["P"] <= 1)
            & (chunk["CHR"] >= 1)
            & (chunk["CHR"] <= 22)
            & (chunk["BP"] > 0)
        ].copy()

        if not chunk.empty:
            chunk["CHR"] = chunk["CHR"].astype(int)
            chunk["BP"] = chunk["BP"].astype(int)
            chunks.append(chunk)

    if not chunks:
        raise ValueError(f"No valid GWAS rows found in {file_path.name}")

    df = pd.concat(chunks, ignore_index=True)

    print(f"Valid rows loaded: {len(df):,}")
    print("-" * 60)

    return df


def plot_manhattan(df: pd.DataFrame, title: str, output_file: Path) -> None:
    df = df.copy()
    df["minus_log10_p"] = -np.log10(df["P"])
    df = df.sort_values(["CHR", "BP"]).reset_index(drop=True)

    current_offset = 0
    tick_positions = []
    tick_labels = []
    plot_frames = []

    for chrom in sorted(df["CHR"].unique()):
        chrom_df = df[df["CHR"] == chrom].copy()
        chrom_df["plot_position"] = chrom_df["BP"] + current_offset

        tick_positions.append(
            chrom_df["plot_position"].min()
            + (chrom_df["plot_position"].max() - chrom_df["plot_position"].min()) / 2
        )
        tick_labels.append(str(chrom))

        current_offset = chrom_df["plot_position"].max()
        plot_frames.append(chrom_df)

    plot_df = pd.concat(plot_frames, ignore_index=True)

    plt.figure(figsize=(14, 6))

    for chrom in sorted(plot_df["CHR"].unique()):
        chrom_df = plot_df[plot_df["CHR"] == chrom]
        plt.scatter(
            chrom_df["plot_position"],
            chrom_df["minus_log10_p"],
            s=4,
            alpha=0.7,
        )

    plt.axhline(
        -np.log10(GWAS_THRESHOLD),
        linestyle="--",
        linewidth=1,
        label=f"GWAS threshold: {GWAS_THRESHOLD}"
    )

    plt.xlabel("Chromosome")
    plt.ylabel("-log10(P-value)")
    plt.title(title)
    plt.xticks(tick_positions, tick_labels)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def plot_qq(df: pd.DataFrame, title: str, output_file: Path) -> None:
    p_values = df["P"].dropna()
    p_values = p_values[(p_values > 0) & (p_values <= 1)].sort_values()

    observed = -np.log10(p_values)
    expected = -np.log10(np.arange(1, len(p_values) + 1) / (len(p_values) + 1))

    plt.figure(figsize=(6, 6))
    plt.scatter(expected, observed, s=5, alpha=0.6)

    max_value = max(expected.max(), observed.max())
    plt.plot([0, max_value], [0, max_value], linestyle="--", linewidth=1)

    plt.xlabel("Expected -log10(P-value)")
    plt.ylabel("Observed -log10(P-value)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def plot_volcano(df: pd.DataFrame, title: str, output_file: Path) -> None:
    df = df.copy()
    df["minus_log10_p"] = -np.log10(df["P"])

    plt.figure(figsize=(8, 6))
    plt.scatter(df["BETA"], df["minus_log10_p"], s=5, alpha=0.6)

    plt.axhline(
        -np.log10(GWAS_THRESHOLD),
        linestyle="--",
        linewidth=1,
        label=f"GWAS threshold: {GWAS_THRESHOLD}"
    )

    plt.axvline(0, linestyle="--", linewidth=1)
    plt.xlabel("Beta effect size")
    plt.ylabel("-log10(P-value)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def plot_count_bar(
    table_file: Path,
    x_col: str,
    y_col: str,
    title: str,
    output_file: Path,
) -> None:
    if not table_file.exists():
        print(f"Skipped missing table: {table_file.name}")
        return

    df = pd.read_csv(table_file, sep="\t")

    if x_col not in df.columns or y_col not in df.columns:
        raise ValueError(f"Missing columns in {table_file.name}: {x_col}, {y_col}")

    df[y_col] = pd.to_numeric(df[y_col], errors="coerce")
    df = df.dropna(subset=[x_col, y_col])

    if df.empty:
        print(f"Skipped empty table: {table_file.name}")
        return

    plt.figure(figsize=(8, 5))
    plt.bar(df[x_col], df[y_col])
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title(title)
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def load_gene_list(file_path: Path) -> set:
    if not file_path.exists():
        print(f"Skipped missing gene list: {file_path.name}")
        return set()

    df = pd.read_csv(file_path, sep="\t")

    if "GENE" not in df.columns:
        raise ValueError(f"GENE column missing in {file_path.name}")

    return set(df["GENE"].astype(str).str.strip())


def plot_venn(set_a: set, set_b: set, label_a: str, label_b: str, title: str, output_file: Path) -> None:
    plt.figure(figsize=(6, 6))
    venn2([set_a, set_b], set_labels=(label_a, label_b))
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def plot_chromosome_distribution(table_file: Path, title: str, output_file: Path) -> None:
    if not table_file.exists():
        print(f"Skipped missing SNP table: {table_file.name}")
        return

    df = pd.read_csv(table_file, sep="\t")

    chr_col = "chromosome" if "chromosome" in df.columns else "CHR"

    if chr_col not in df.columns:
        print(f"Skipped chromosome plot. No chromosome column in {table_file.name}")
        return

    df[chr_col] = pd.to_numeric(df[chr_col], errors="coerce")
    df = df.dropna(subset=[chr_col])
    df = df[(df[chr_col] >= 1) & (df[chr_col] <= 22)].copy()
    df[chr_col] = df[chr_col].astype(int)

    counts = df[chr_col].value_counts().sort_index()

    plt.figure(figsize=(10, 5))
    plt.bar(counts.index.astype(str), counts.values)
    plt.xlabel("Chromosome")
    plt.ylabel("Number of significant SNPs")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def plot_top_genes(table_file: Path, title: str, output_file: Path, top_n: int = 10) -> None:
    if not table_file.exists():
        print(f"Skipped missing gene table: {table_file.name}")
        return

    df = pd.read_csv(table_file, sep="\t")

    if "GENE" not in df.columns or "P" not in df.columns:
        print(f"Skipped top genes plot. Missing GENE or P in {table_file.name}")
        return

    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df = df.dropna(subset=["GENE", "P"])
    df = df[(df["P"] > 0) & (df["P"] <= 1)].copy()

    if df.empty:
        print(f"Skipped empty gene table: {table_file.name}")
        return

    df = df.sort_values(by="P", ascending=True).head(top_n)
    df["minus_log10_p"] = -np.log10(df["P"])

    plt.figure(figsize=(9, 5))
    plt.bar(df["GENE"], df["minus_log10_p"])
    plt.xlabel("Gene")
    plt.ylabel("-log10(P-value)")
    plt.title(title)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    male_file = RAW_DIR / f"{PHENOTYPE}_male.tsv"
    female_file = RAW_DIR / f"{PHENOTYPE}_female.tsv"

    male_df = load_gwas_for_plot(male_file)
    female_df = load_gwas_for_plot(female_file)

    plot_manhattan(
        male_df,
        f"{PHENOTYPE} Male Manhattan Plot",
        FIGURES_DIR / f"{PHENOTYPE}_male_manhattan.png"
    )

    plot_manhattan(
        female_df,
        f"{PHENOTYPE} Female Manhattan Plot",
        FIGURES_DIR / f"{PHENOTYPE}_female_manhattan.png"
    )

    plot_qq(
        male_df,
        f"{PHENOTYPE} Male QQ Plot",
        FIGURES_DIR / f"{PHENOTYPE}_male_qq.png"
    )

    plot_qq(
        female_df,
        f"{PHENOTYPE} Female QQ Plot",
        FIGURES_DIR / f"{PHENOTYPE}_female_qq.png"
    )

    plot_volcano(
        male_df,
        f"{PHENOTYPE} Male Volcano Plot",
        FIGURES_DIR / f"{PHENOTYPE}_male_volcano.png"
    )

    plot_volcano(
        female_df,
        f"{PHENOTYPE} Female Volcano Plot",
        FIGURES_DIR / f"{PHENOTYPE}_female_volcano.png"
    )

    plot_count_bar(
        TABLES_DIR / f"{PHENOTYPE}_summary.tsv",
        "group",
        "count",
        f"{PHENOTYPE} SNP Comparison",
        FIGURES_DIR / f"{PHENOTYPE}_snp_summary_barplot.png"
    )

    plot_count_bar(
        TABLES_DIR / "gene_summary.tsv",
        "group",
        "count",
        f"{PHENOTYPE} Gene Comparison",
        FIGURES_DIR / f"{PHENOTYPE}_gene_summary_barplot.png"
    )

    male_snps = pd.read_csv(FILTERED_DIR / f"{PHENOTYPE}_male_significant_snps.tsv", sep="\t")
    female_snps = pd.read_csv(FILTERED_DIR / f"{PHENOTYPE}_female_significant_snps.tsv", sep="\t")

    plot_venn(
        set(male_snps["SNP"].astype(str)),
        set(female_snps["SNP"].astype(str)),
        "Male",
        "Female",
        f"{PHENOTYPE} Significant SNPs Venn Diagram",
        FIGURES_DIR / f"{PHENOTYPE}_snp_venn.png"
    )

    male_genes = load_gene_list(GENE_LISTS_DIR / "male_specific_genes.txt")
    female_genes = load_gene_list(GENE_LISTS_DIR / "female_specific_genes.txt")
    shared_genes = load_gene_list(GENE_LISTS_DIR / "shared_genes.txt")

    plot_venn(
        male_genes | shared_genes,
        female_genes | shared_genes,
        "Male",
        "Female",
        f"{PHENOTYPE} Significant Genes Venn Diagram",
        FIGURES_DIR / f"{PHENOTYPE}_gene_venn.png"
    )

    plot_chromosome_distribution(
        FILTERED_DIR / f"{PHENOTYPE}_male_significant_snps.tsv",
        f"{PHENOTYPE} Male Significant SNPs by Chromosome",
        FIGURES_DIR / f"{PHENOTYPE}_male_chr_distribution.png"
    )

    plot_chromosome_distribution(
        FILTERED_DIR / f"{PHENOTYPE}_female_significant_snps.tsv",
        f"{PHENOTYPE} Female Significant SNPs by Chromosome",
        FIGURES_DIR / f"{PHENOTYPE}_female_chr_distribution.png"
    )

    plot_top_genes(
        TABLES_DIR / f"{PHENOTYPE}_male_specific_gene_table.tsv",
        f"{PHENOTYPE} Top Male-Specific Genes",
        FIGURES_DIR / f"{PHENOTYPE}_top_male_specific_genes.png"
    )

    plot_top_genes(
        TABLES_DIR / f"{PHENOTYPE}_female_specific_gene_table.tsv",
        f"{PHENOTYPE} Top Female-Specific Genes",
        FIGURES_DIR / f"{PHENOTYPE}_top_female_specific_genes.png"
    )

    plot_top_genes(
        TABLES_DIR / f"{PHENOTYPE}_shared_gene_table.tsv",
        f"{PHENOTYPE} Top Shared Genes",
        FIGURES_DIR / f"{PHENOTYPE}_top_shared_genes.png"
    )

    print("-" * 60)
    print("Visualization completed successfully.")


if __name__ == "__main__":
    main()
