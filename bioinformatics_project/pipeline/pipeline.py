import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

from .differential_expression import differential_expression_analysis, find_common_degs
from .ppi_network import construct_ppi_network
from .enrichment import run_enrichment_analysis
from .visualization import perform_pca, plot_volcano, plot_ppi_network, plot_enrichment
from .validation import validate_core_genes

plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")


class BioinformaticsPipeline:
    
    def __init__(self, random_state=42):
        self.random_state = random_state
        np.random.seed(random_state)
        self.degs_disease1 = None
        self.degs_disease2 = None
        self.co_degs = None
        self.core_genes = None
    
    def perform_pca(self, expression_data, labels, title="PCA Analysis", filename=None):
        return perform_pca(expression_data, labels, title, filename, self.random_state)
    
    def differential_expression_analysis(self, expression_data, labels, 
                                        pvalue_threshold=0.05, 
                                        logfc_threshold=0.5):
        return differential_expression_analysis(expression_data, labels, 
                                               pvalue_threshold, logfc_threshold)
    
    def plot_volcano(self, deg_results, title="Volcano Plot", filename=None):
        return plot_volcano(deg_results, title, filename)
    
    def find_common_degs(self, degs1, degs2, disease1_name="AD", disease2_name="OA"):
        fig, co_degs = find_common_degs(degs1, degs2, disease1_name, disease2_name)
        self.co_degs = co_degs
        return fig, co_degs
    
    def construct_ppi_network(self, genes, k=10):
        G, core_genes = construct_ppi_network(genes, k)
        self.core_genes = core_genes
        return G, core_genes
    
    def plot_ppi_network(self, G, core_genes, filename=None):
        return plot_ppi_network(G, core_genes, filename, self.random_state)
    
    def run_enrichment_analysis(self, genes, cutoff=0.05):
        return run_enrichment_analysis(genes, cutoff)
    
    def plot_enrichment(self, enrichment_results, top_n=10, filename=None):
        return plot_enrichment(enrichment_results, top_n, filename)
    
    def validate_core_genes(self, expression_data, labels, core_genes, filename_prefix=None):
        return validate_core_genes(expression_data, labels, core_genes, filename_prefix)
    
    def generate_report(self, output_file="analysis_report.txt"):
        print(f"\n{'='*70}")
        print("Generating Analysis Report")
        print(f"{'='*70}")
        
        report = []
        report.append("="*70)
        report.append("BIOINFORMATICS ANALYSIS REPORT")
        report.append("Gene Expression Analysis Pipeline")
        report.append("Based on Liu et al. 2025 (PMC11805404)")
        report.append("="*70)
        report.append("")
        
        report.append("STUDY OVERVIEW:")
        report.append("-" * 70)
        report.append("This analysis implements the bioinformatics pipeline from:")
        report.append("'Characterization of gene expression profiles in Alzheimer's")
        report.append("disease and osteoarthritis: A bioinformatics study'")
        report.append("")
        
        report.append("METHODOLOGY:")
        report.append("-" * 70)
        report.append("1. Data: GEO training/validation sets (AD: GSE5281, GSE28146, GSE29378, GSE122063; OA: GSE55235, GSE206848, GSE82107, GSE55457)")
        report.append("2. Quality Control: PCA analysis for batch effects")
        report.append("3. Differential Expression: DEG identification (adj. P < 0.05, |log2FC| > 0.5)")
        report.append("4. Co-DEG Analysis: Identification of common DEGs between conditions")
        report.append("5. PPI Network: Protein-protein interaction network construction")
        report.append("6. Enrichment: GO/KEGG pathway enrichment (Enrichr)")
        report.append("7. Validation: ROC curve analysis for biomarker potential")
        report.append("")
        
        if self.co_degs:
            report.append("KEY FINDINGS:")
            report.append("-" * 70)
            report.append(f"Common upregulated DEGs: {len(self.co_degs['upregulated'])}")
            report.append(f"Common downregulated DEGs: {len(self.co_degs['downregulated'])}")
            report.append(f"Total Co-DEGs: {len(self.co_degs['all'])}")
            report.append("")
        
        if self.core_genes:
            report.append("CORE GENES IDENTIFIED:")
            report.append("-" * 70)
            for i, gene in enumerate(self.core_genes, 1):
                report.append(f"{i}. {gene}")
            report.append("")
        
        report.append("BIOLOGICAL INTERPRETATION:")
        report.append("-" * 70)
        report.append("The identified core genes represent potential:")
        report.append("• Biomarkers for disease diagnosis")
        report.append("• Therapeutic targets for drug development")
        report.append("• Shared pathological mechanisms between conditions")
        report.append("")
        
        report.append("SUGGESTED NEXT STEPS:")
        report.append("-" * 70)
        report.append("1. Validate findings in independent datasets")
        report.append("2. Perform experimental validation (qRT-PCR, Western blot)")
        report.append("3. Investigate gene functions through literature review")
        report.append("4. Design functional studies to elucidate mechanisms")
        report.append("5. Explore therapeutic potential through drug repurposing")
        report.append("")
        
        report.append("="*70)
        report.append("Analysis completed successfully!")
        report.append("="*70)
        
        with open(output_file, 'w') as f:
            f.write('\n'.join(report))
        
        print(f"✓ Report saved: {output_file}")
        
        print("\n" + '\n'.join(report))
