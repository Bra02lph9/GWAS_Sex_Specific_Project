# Biological Interpretation — Sex-Specific Network-Based GWAS (CAD)

---

## 1. Male-Specific Network

### Top Hub Genes (by Degree Centrality)

| Gene | Degree Centrality | Betweenness |
|------|-------------------|-------------|
| STAT3 | 0.0402 | 0.1442 |
| IL1B | 0.0333 | 0.0596 |
| BRCA1 | 0.0325 | 0.0651 |
| RPS3 | 0.0308 | 0.0308 |
| H3-3B | 0.0308 | 0.0440 |

### Top Enriched Pathways (by Adjusted P-value)

| Pathway | Source | Adjusted P-value | Intersection Size |
|---------|--------|-----------------|-------------------|
| Response to stimulus | GO:BP | 4.69e-20 | 934 |
| Intracellular signal transduction | GO:BP | 2.14e-14 | 366 |
| Membrane | GO:CC | 5.62e-14 | 979 |
| Intracellular signaling cassette | GO:BP | 3.76e-12 | 248 |
| Regulation of signal transduction | GO:BP | 5.04e-11 | 358 |
| Regulation of cell communication | GO:BP | 9.87e-11 | 399 |
| Regulation of signaling | GO:BP | 1.16e-10 | 398 |
| Response to stress | GO:BP | 3.39e-10 | 443 |
| Regulation of response to stimulus | GO:BP | 7.62e-10 | 442 |
| Vesicle | GO:CC | 9.68e-10 | 441 |

### Interpretation

The male-specific network is centered around **inflammatory and stress-response signaling**. STAT3 and IL1B, the two most connected hub genes, are well-established mediators of inflammatory cascades. STAT3 is a transcription factor activated by cytokines and growth factors, while IL1B is a pro-inflammatory cytokine directly implicated in atherosclerosis progression.

The top enriched pathways — response to stimulus, intracellular signal transduction, and regulation of signaling — are consistent with these hub genes, forming a coherent biological signal. This convergence between network topology and pathway enrichment suggests that **inflammatory signaling is a central male-specific mechanism in CAD genetic risk**.

The pathway enrichment is highly significant (best adjusted p-value: 4.69e-20), indicating a strong and reliable biological signal in the male cohort.

---

## 2. Female-Specific Network

### Top Hub Genes (by Degree Centrality)

| Gene | Degree Centrality | Betweenness |
|------|-------------------|-------------|
| JAK2 | 0.0290 | 0.1254 |
| UTP4 | 0.0228 | 0.0260 |
| MTREX | 0.0228 | 0.0554 |
| CS | 0.0228 | 0.1109 |
| STAT2 | 0.0207 | 0.0154 |

### Top Enriched Pathways (by Adjusted P-value)

| Pathway | Source | Adjusted P-value | Intersection Size |
|---------|--------|-----------------|-------------------|
| Regulation of DNA-templated transcription | GO:BP | 0.00192 | 213 |
| Regulation of RNA biosynthetic process | GO:BP | 0.00248 | 213 |
| Regulation of RNA metabolic process | GO:BP | 0.00290 | 226 |
| DNA-templated transcription | GO:BP | 0.00313 | 219 |
| Regulation of nucleobase-containing compound metabolic process | GO:BP | 0.00353 | 241 |
| Regulation of multicellular organismal development | GO:BP | 0.00953 | 101 |
| Skeletal muscle; myocytes [≥Low] | HPA | 0.01198 | 388 |
| Skeletal muscle | HPA | 0.01198 | 388 |
| Lung; macrophages [≥Low] | HPA | 0.01708 | 436 |
| Transport vesicle membrane | GO:CC | 0.01765 | 27 |

### Interpretation

The female-specific network is dominated by **transcriptional regulation and RNA biosynthesis**. JAK2, the top hub gene by both degree centrality and betweenness, is a kinase central to the JAK-STAT signaling pathway and plays a key role in regulating gene expression in response to cytokines. STAT2, also a top hub, is a direct transcriptional effector downstream of JAK2.

The top enriched pathways — regulation of DNA-templated transcription, RNA biosynthetic process, and RNA metabolic process — are directly consistent with the function of JAK2 and STAT2, confirming that **transcriptional regulation is the primary female-specific mechanism in CAD genetic risk**.

```JAK2 (Janus kinase 2) and STAT2 (Signal transducer and activator of transcription 2) have a strong relationship with hormones, particularly in mediating signals from hormone receptors to the cell nucleus. While JAK2 is a central, well-documented actor in hormone signaling, STAT2 is often involved alongside other STAT proteins (like STAT1, 3, and 5) in mediating responses to cytokine-like hormones.```

---

## 3. Shared Network

### Top Hub Genes (by Degree Centrality)

| Gene | Degree Centrality | Betweenness |
|------|-------------------|-------------|
| AKT1 | 0.0746 | 0.1391 |
| APOE | 0.0634 | 0.1255 |
| FN1 | 0.0597 | 0.1404 |
| APOB | 0.0560 | 0.0463 |
| UBC | 0.0560 | 0.1358 |

### Top Enriched Pathways (by Adjusted P-value)

| Pathway | Source | Adjusted P-value | Intersection Size |
|---------|--------|-----------------|-------------------|
| Protein binding | GO:MF | 4.50e-15 | 867 |
| Biological regulation | GO:BP | 1.93e-10 | 701 |
| Regulation of biological process | GO:BP | 4.80e-09 | 677 |
| Regulation of cellular process | GO:BP | 2.99e-08 | 654 |
| Regulation of primary metabolic process | GO:BP | 5.27e-05 | 316 |
| Localization | GO:BP | 5.36e-04 | 324 |
| Multicellular organismal process | GO:BP | 6.45e-04 | 414 |
| Developmental process | GO:BP | 6.47e-04 | 375 |
| System development | GO:BP | 9.19e-04 | 248 |
| Positive regulation of biological process | GO:BP | 1.36e-03 | 354 |

### Interpretation

The shared network captures the **classical core biology of CAD**, driven by well-established disease genes. AKT1 is a central survival kinase regulating cell proliferation and apoptosis. APOE is the primary lipoprotein involved in cholesterol transport and is one of the most replicated CAD risk genes. FN1 encodes fibronectin, a key extracellular matrix protein involved in vascular remodeling. APOB is the primary apolipoprotein of LDL particles. UBC encodes ubiquitin C, a central regulator of protein degradation and cellular homeostasis.

The top shared pathways — protein binding and broad biological regulation — reflect the general regulatory roles of these hub genes across multiple cellular processes, consistent with their known pleiotropic effects in cardiovascular disease.

---

## 4. Summary — Key Biological Finding

| Level | Males | Females | Shared |
|-------|-------|---------|--------|
| Top hub gene | STAT3 | JAK2 | AKT1 |
| Core mechanism | Inflammatory signaling | Transcriptional regulation | Classical CAD biology |
| Top pathway | Response to stimulus | Regulation of transcription | Protein binding |
| Signal strength | Very strong (p ~ 1e-20) | Moderate (p ~ 0.002) | Strong (p ~ 1e-15) |

The central finding of this analysis is that **CAD genetic risk operates through distinct molecular mechanisms in males and females**. Males show a convergence of hub genes and enriched pathways around inflammatory and stress-response signaling, while females show convergence around transcriptional regulation via the JAK-STAT axis. The shared network recovers established CAD biology including lipid metabolism (APOE, APOB) and cell survival signaling (AKT1), providing biological validation of the analysis.

This sex-specific divergence in molecular mechanisms has potential implications for sex-stratified therapeutic strategies in CAD.