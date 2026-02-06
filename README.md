# Bioinformatics Gene Expression Analysis Pipeline

## Project Overview

This bioinformatics pipeline implements the methodology from the research paper:

**"Characterization of gene expression profiles in Alzheimer's disease and osteoarthritis: A bioinformatics study"**  
*Liu et al. 2025 (PMC11805404)*

https://pmc.ncbi.nlm.nih.gov/articles/PMC11805404/

Published in PLoS One, this study explores the common molecular mechanisms between Alzheimer's Disease (AD) and Osteoarthritis (OA) through comprehensive gene expression analysis.

---

## Differences

The divergence between our results and the paper stems from two critical methodological differences: database evolution and data scope.

First, Protein-Protein Interaction (PPI) databases are not static; they grow daily. The paper relied on a snapshot of biological knowledge from the past. Your analysis utilizes a modern "map" where immune genes like CD74 and CSF1R likely have newly documented interactions. These recent additions can drastically increase a gene's connectivity score, causing it to leapfrog over the structural genes found in the older study.

Second, using a smaller data subset with different significance parameters fundamentally shifts the biological focus. Large datasets tend to amplify the strongest, most common signals—in this case, the ubiquitous structural matrix genes (COL1A2, VCAN) that hold tissue together. By analyzing a subset, we effectively "zoomed in." We reduced the global noise of the structural tissue, allowing the specific, highly interconnected immune signal to rise to the top.

This does not make our findings "wrong"; rather, it makes them context-specific. While the paper described the general "building architecture" (stroma), your parameters successfully isolated the "occupants" (immune cells). You have captured a valid, high-resolution view of the inflammatory microenvironment that was likely diluted in the paper’s broader dataset.

## Authors

**Group Members:**
- Tako Gelashvili
- Nika Gvalia
- Giorgi Bachaliashvili
---

## HOW TO RUN
pip install -r requirements.txt 
python pipeline/main.py

## Scientific Background

### Research Question
What are the shared molecular pathways and potential biomarkers between Alzheimer's disease and osteoarthritis?
---

## Pipeline Methodology

This implementation follows the exact workflow from Liu et al. 2025:

### 1. **Data Loading**
   - **Real data (default):** GEO datasets from the paper (PMC11805404):
     - **AD training:** GSE5281, GSE28146, GSE29378 (brain)
     - **AD validation:** GSE122063 (brain)
     - **OA training:** GSE55235, GSE206848, GSE82107 (synovial)
     - **OA validation:** GSE55457 (synovial)
   - Requires `GEOparse` and network access for first run (data cached in `geo_cache/`).

### 2. **Quality Control (PCA)**
   - Principal Component Analysis to assess sample clustering
   - Batch effect detection
   - Visualization of sample distribution

### 3. **Differential Expression Analysis**
   - Statistical testing using t-tests
   - Multiple testing correction (Benjamini-Hochberg FDR)
   - Thresholds: adj. P < 0.05 and |log2(FC)| > 0.5
   - Volcano plot visualization

### 4. **Common DEG Identification**
   - Venn diagram analysis
   - Identification of upregulated and downregulated Co-DEGs
   - Statistical overlap assessment

### 5. **Protein-Protein Interaction Network**
   - Network construction (simulated; real analysis would use STRING database)
   - Hub gene identification via degree centrality
   - Core gene selection

### 6. **GO/KEGG Enrichment Analysis**
   - Biological Process, Molecular Function, Cellular Component
   - KEGG pathway enrichment
   - Visualization of top enriched terms

### 7. **Biomarker Validation**
   - Expression level comparison (box plots)
   - ROC curve analysis
   - AUC calculation for diagnostic potential

### 8. **Report Generation**
   - Comprehensive text report
   - Summary of findings
   - Biological interpretation

---

## Installation & Setup

### Prerequisites
- Python 3.8 or higher
- pip package manager

### Step 1: Install Dependencies

```bash
pip install -r requirements.txt
```

### Step 2: Create Results Directory

The script will automatically create this, but you can also do it manually:

```bash
mkdir results
```

---

## Usage

### Running the Complete Pipeline (Real GEO Data)

By default the pipeline loads the **real GEO datasets** from the paper:

```bash
pip install -r requirements.txt   # includes GEOparse
python pipeline/main.py
```

First run will download SOFT files from NCBI GEO into `geo_cache/` (requires network). Subsequent runs use the cache.

```bash
python pipeline/main.py 
```

Roughly **~10–30 MB total** and **~1–3 minutes** on first run. The pipeline still runs DEG, Co-DEG, PPI, and enrichment; validation will use these training samples if you skip validation (see below).

To also **skip validation sets** (no GSE122063, GSE55457) and save more download:

```bash
python bioinformatics_pipeline.py
```

Then only the two small training GSEs are downloaded; validation plots use the same training data.

### Optional: GEO Cache Directory

```bash
export GEO_CACHE_DIR=/path/to/cache
python bioinformatics_pipeline.py
```

Later runs use the cache and are much faster.

---

## Output Files

The pipeline generates comprehensive results in the `results/` directory:

### Visualizations
| File | Description |
|------|-------------|
| `pca_alzheimers.png` | PCA plot for AD dataset |
| `pca_osteoarthritis.png` | PCA plot for OA dataset |
| `volcano_alzheimers.png` | Volcano plot showing AD DEGs |
| `volcano_osteoarthritis.png` | Volcano plot showing OA DEGs |
| `venn_diagram.png` | Venn diagram of common DEGs |
| `ppi_network.png` | Protein-protein interaction network |
| `enrichment_analysis.png` | GO/KEGG enrichment results |
| `validation_ad_expression.png` / `validation_oa_expression.png` | Core gene expression (AD/OA validation sets) |
| `validation_ad_roc.png` / `validation_oa_roc.png` | ROC curves for AD/OA validation |

### Data Files (CSV)
| File | Description |
|------|-------------|
| `degs_alzheimers.csv` | Complete DEG statistics for AD |
| `degs_osteoarthritis.csv` | Complete DEG statistics for OA |
| `common_degs.csv` | List of common DEGs between AD and OA |
| `core_genes.csv` | Ranked list of core/hub genes |
| `enrichment_GO_BP.csv` | Gene Ontology Biological Process results |
| `enrichment_GO_MF.csv` | Gene Ontology Molecular Function results |
| `enrichment_GO_CC.csv` | Gene Ontology Cellular Component results |
| `enrichment_KEGG.csv` | KEGG pathway enrichment results |
| `auc_scores.csv` | AUC scores for core gene validation |

### Report
| File | Description |
|------|-------------|
| `analysis_report.txt` | Comprehensive analysis summary |


## Understanding the Results

### Key Metrics

1. **DEG Count**: Number of significantly differentially expressed genes
2. **Co-DEGs**: Genes showing differential expression in BOTH conditions
3. **Core Genes**: Hub genes with high connectivity in PPI network
4. **AUC Scores**: Diagnostic accuracy (0.5 = random, 1.0 = perfect)

### Interpreting Visualizations

#### PCA Plot
- Clustering = Good separation between control and disease
- Overlap = Potential confounding factors or disease heterogeneity

#### Volcano Plot
- Red points (right side) = Upregulated genes
- Blue points (left side) = Downregulated genes
- Higher points = More statistically significant

#### Venn Diagram
- Overlap area = Common DEGs between conditions
- Larger overlap = More shared molecular mechanisms

#### ROC Curves
- Higher AUC = Better diagnostic potential
- AUC > 0.7 = Acceptable biomarker
- AUC > 0.8 = Good biomarker
- AUC > 0.9 = Excellent biomarker

### Accessing Real Datasets

The original paper used GEO datasets:
- **AD**: GSE5281, GSE28146, GSE29378, GSE122063
- **OA**: GSE55235, GSE55457, GSE206848, GSE82107

Download from: https://www.ncbi.nlm.nih.gov/geo/

---

### Simplified Analyses
- **PPI Network**: Real analysis would query STRING database for validated interactions
- **Enrichment**: Real analysis would use tools like DAVID, Enrichr, or clusterProfiler
- **No p-value correction** in some visualizations for clarity

### Statistical Considerations
- Multiple testing is a major concern in genomics
- Effect sizes (fold changes) are as important as p-values
- Biological replication is essential for validation
- 
---


