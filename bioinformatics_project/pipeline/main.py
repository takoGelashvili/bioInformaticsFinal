import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Allow running this file directly: ensure project root is on sys.path BEFORE package imports
if __name__ == "__main__":
    try:
        # This file lives at: <project_root>/bioinformatics_project/pipeline/main.py
        _THIS_FILE = os.path.abspath(__file__)
        _PIPELINE_DIR = os.path.dirname(_THIS_FILE)                 # .../bioinformatics_project/pipeline
        _PKG_DIR = os.path.dirname(_PIPELINE_DIR)                   # .../bioinformatics_project
        _PROJECT_ROOT = os.path.dirname(_PKG_DIR)                   # .../ (repo root)
        if _PROJECT_ROOT not in sys.path:
            sys.path.insert(0, _PROJECT_ROOT)
    except Exception:
        # If anything goes wrong, proceed — imports below may still work in packaged context
        pass

# Use absolute package imports to avoid relative-import issues when run via script
from bioinformatics_project.pipeline.pipeline import BioinformaticsPipeline
from bioinformatics_project.pipeline.data_loader import load_paper_geo_data


def main(geo_cache_dir=None):
    print("\n" + "="*70)
    print("BIOINFORMATICS PIPELINE - GENE EXPRESSION ANALYSIS")
    print("Based on: Liu et al. 2025 (PMC11805404)")
    print("="*70)
    print("\nThis pipeline implements the methodology from the paper:")
    print("'Characterization of gene expression profiles in")
    print("Alzheimer's disease and osteoarthritis'")
    print("\n" + "="*70 + "\n")
    
    pipeline = BioinformaticsPipeline(random_state=42)
    
    print("\n### STEP 1: DATA LOADING ###")
    try:
        from bioinformatics_project.geo_loader.utils import get_ad_gse_ids, get_oa_gse_ids
        ad_ids, oa_ids = get_ad_gse_ids(), get_oa_gse_ids()
        print(f"AD GSE IDs: {', '.join(ad_ids)}")
        print(f"OA GSE IDs: {', '.join(oa_ids)}")
        print("(Edit geo_ids.txt or set AD_GSE_IDS / OA_GSE_IDS to change. Cache used if present; else download.)")
    except Exception:
        pass
    print("Loading GEO data...")
    (expr_ad, labels_ad, expr_oa, labels_oa,
     expr_ad_val, labels_ad_val, expr_oa_val, labels_oa_val) = load_paper_geo_data(geo_cache_dir)
    print(f"  AD training: {expr_ad.shape[0]} samples, {expr_ad.shape[1]} genes")
    print(labels_ad.value_counts().to_string())
    print(f"  OA training: {expr_oa.shape[0]} samples, {expr_oa.shape[1]} genes")
    print(labels_oa.value_counts().to_string())
    
    print("\n### STEP 2: QUALITY CONTROL - PCA ###")
    
    fig_pca_ad, pca_ad = pipeline.perform_pca(
        expr_ad, labels_ad, 
        title="PCA: Alzheimer's Disease Dataset",
        filename="results/pca_alzheimers.png"
    )
    
    fig_pca_oa, pca_oa = pipeline.perform_pca(
        expr_oa, labels_oa,
        title="PCA: Osteoarthritis Dataset",
        filename="results/pca_osteoarthritis.png"
    )
    
    print("\n### STEP 3: DIFFERENTIAL EXPRESSION ANALYSIS ###")
    
    degs_ad = pipeline.differential_expression_analysis(
        expr_ad, labels_ad,
        pvalue_threshold=0.05,
        logfc_threshold=0.5
    )
    
    degs_oa = pipeline.differential_expression_analysis(
        expr_oa, labels_oa,
        pvalue_threshold=0.05,
        logfc_threshold=0.5
    )
    
    degs_ad.to_csv("results/degs_alzheimers.csv", index=False)
    degs_oa.to_csv("results/degs_osteoarthritis.csv", index=False)
    print("✓ DEG results saved to CSV files")
    
    print("\n### STEP 4: VOLCANO PLOTS ###")
    
    fig_volcano_ad = pipeline.plot_volcano(
        degs_ad,
        title="Volcano Plot: Alzheimer's Disease",
        filename="results/volcano_alzheimers.png"
    )
    
    fig_volcano_oa = pipeline.plot_volcano(
        degs_oa,
        title="Volcano Plot: Osteoarthritis",
        filename="results/volcano_osteoarthritis.png"
    )
    
    print("\n### STEP 5: COMMON DEG ANALYSIS ###")
    
    pipeline.degs_disease1 = degs_ad
    pipeline.degs_disease2 = degs_oa
    
    fig_venn, co_degs = pipeline.find_common_degs(
        degs_ad, degs_oa,
        disease1_name="AD",
        disease2_name="OA"
    )
    plt.savefig("results/venn_diagram.png", dpi=300, bbox_inches='tight')
    print("✓ Venn diagram saved: results/venn_diagram.png")
    
    co_degs_df = pd.DataFrame({
        'Gene': co_degs['all'],
        'Regulation': ['Upregulated'] * len(co_degs['upregulated']) + 
                     ['Downregulated'] * len(co_degs['downregulated'])
    })
    co_degs_df.to_csv("results/common_degs.csv", index=False)
    print("✓ Common DEGs saved: results/common_degs.csv")
    
    print("\n### STEP 6: PPI NETWORK CONSTRUCTION ###")
    
    core_genes = []
    G = None
    
    if len(co_degs['all']) > 0:
        G, core_genes = pipeline.construct_ppi_network(
            co_degs['all'][:min(len(co_degs['all']), 30)],
            k=10
        )
        
        fig_ppi = pipeline.plot_ppi_network(
            G, core_genes,
            filename="results/ppi_network.png"
        )
        
        core_genes_df = pd.DataFrame({
            'Gene': core_genes,
            'Rank': range(1, len(core_genes) + 1)
        })
        core_genes_df.to_csv("results/core_genes.csv", index=False)
        print("✓ Core genes saved: results/core_genes.csv")
    else:
        print("⚠ Warning: No common DEGs found. Skipping PPI network construction.")
    
    print("\n### STEP 7: ENRICHMENT ANALYSIS ###")
    
    enrichment_results = pipeline.run_enrichment_analysis(co_degs['all'])
    
    fig_enrichment = pipeline.plot_enrichment(
        enrichment_results,
        top_n=10,
        filename="results/enrichment_analysis.png"
    )
    
    for category, df in enrichment_results.items():
        df.to_csv(f"results/enrichment_{category}.csv", index=False)
    print("✓ Enrichment results saved to CSV files")
    
    print("\n### STEP 8: CORE GENE VALIDATION ###")
    
    if len(core_genes) > 0:
        has_ad_val = expr_ad_val is not None and getattr(expr_ad_val, "shape", (0,))[0] > 0
        has_oa_val = expr_oa_val is not None and getattr(expr_oa_val, "shape", (0,))[0] > 0
        if has_ad_val or has_oa_val:
            auc_ad, auc_oa = {}, {}
            if has_ad_val:
                core_in_ad = [g for g in core_genes if g in expr_ad_val.columns]
                if core_in_ad:
                    _, _, auc_ad = pipeline.validate_core_genes(
                        expr_ad_val, labels_ad_val, core_in_ad,
                        filename_prefix="results/validation_ad"
                    )
            if has_oa_val:
                core_in_oa = [g for g in core_genes if g in expr_oa_val.columns]
                if core_in_oa:
                    _, _, auc_oa = pipeline.validate_core_genes(
                        expr_oa_val, labels_oa_val, core_in_oa,
                        filename_prefix="results/validation_oa"
                    )
            rows = []
            for g in core_genes:
                rows.append({
                    "Gene": g,
                    "AUC_AD": auc_ad.get(g, np.nan),
                    "AUC_OA": auc_oa.get(g, np.nan),
                })
            auc_df = pd.DataFrame(rows)
            auc_df.to_csv("results/auc_scores.csv", index=False)
        else:
            expr_combined = pd.concat([expr_ad, expr_oa], axis=0)
            labels_combined = pd.concat([labels_ad, labels_oa], axis=0)
            fig_expr, fig_roc, auc_scores = pipeline.validate_core_genes(
                expr_combined, labels_combined, core_genes,
                filename_prefix="results/validation"
            )
            auc_df = pd.DataFrame(list(auc_scores.items()), columns=["Gene", "AUC"])
            auc_df = auc_df.sort_values("AUC", ascending=False)
            auc_df.to_csv("results/auc_scores.csv", index=False)
        print("✓ AUC scores saved: results/auc_scores.csv")
    else:
        print("⚠ Warning: No core genes available for validation.")
    
    print("\n### STEP 9: GENERATING FINAL REPORT ###")
    
    pipeline.generate_report("results/analysis_report.txt")
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE!")
    print("="*70)
    print("\nAll results have been saved to the 'results/' directory:")
    print("  • PCA plots")
    print("  • Volcano plots")
    print("  • Venn diagrams")
    print("  • PPI network visualization")
    print("  • Enrichment analysis plots")
    print("  • Gene expression validation")
    print("  • ROC curves")
    print("  • CSV files with detailed results")
    print("  • Comprehensive analysis report")
    print("\n" + "="*70)
    
    plt.show()


if __name__ == "__main__":
    os.makedirs("results", exist_ok=True)
    main(geo_cache_dir=os.environ.get("GEO_CACHE_DIR"))
