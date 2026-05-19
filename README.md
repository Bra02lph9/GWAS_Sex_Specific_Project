# Sex-Specific GWAS Analysis Pipeline

A bioinformatics pipeline for analyzing sex-specific genetic mechanisms in complex diseases using GWAS, pathway enrichment, network analysis, and hub-gene identification.

The pipeline can be reused for multiple diseases such as:
- CAD
- Diabetes
- Alzheimer’s disease
- Hypertension
- Cancer-related GWAS
- Other complex traits

---

# Workflow

```text
Raw GWAS Data
    ↓
Preprocessing & SNP Filtering
    ↓
Male/Female SNP Comparison
    ↓
FUMA / MAGMA Analysis
    ↓
Pathway Enrichment
    ↓
STRING / Cytoscape Networks
    ↓
Hub Gene Analysis
    ↓
Visualization & Biological Interpretation
```

---

# Main Tools

- Python
- Pandas
- Matplotlib
- NetworkX
- FUMA
- MAGMA
- STRING
- Cytoscape
- g:Profiler

---

# Run Pipeline

## Pre-FUMA

```bash
python run_pre_fuma.py
```

Includes:
- preprocessing
- SNP filtering
- SNP comparison
- FUMA input preparation

---

## External Tools

Manual analysis using:
- FUMA
- g:Profiler
- STRING
- Cytoscape

Save outputs inside:

```text
3_tools_results/
```

---

## Post-Tools

```bash
python run_post_tools.py
```

Includes:
- gene extraction
- pathway analysis
- hub-gene analysis
- visualization generation

---

# Generated Results

- Manhattan plots
- QQ plots
- Volcano plots
- Venn diagrams
- Pathway enrichment dotplots
- Hub-gene plots
- Cytoscape-ready networks

Results are saved in:

```text
4_results/
```

---

# Authors

Ilyass Ghouaghou & IBrahim Rmili
