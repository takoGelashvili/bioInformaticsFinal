import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.metrics import roc_curve, auc


def validate_core_genes(expression_data, labels, core_genes, filename_prefix=None):
    print(f"\n{'='*70}")
    print("Validating Core Genes")
    print(f"{'='*70}")
    
    n_genes = len(core_genes)
    n_cols = 5
    n_rows = (n_genes + n_cols - 1) // n_cols
    
    fig1, axes = plt.subplots(n_rows, n_cols, figsize=(18, n_rows * 3))
    axes = axes.flatten() if n_genes > 1 else [axes]
    
    for i, gene in enumerate(core_genes):
        if gene not in expression_data.columns:
            continue
        
        data_to_plot = []
        labels_to_plot = []
        
        for condition in labels.unique():
            mask = labels == condition
            data_to_plot.append(expression_data.loc[mask, gene].values)
            labels_to_plot.append(condition)
        
        bp = axes[i].boxplot(data_to_plot, labels=labels_to_plot, patch_artist=True)
        
        colors = ['#3498db', '#e74c3c', '#2ecc71']
        for patch, color in zip(bp['boxes'], colors[:len(labels_to_plot)]):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        axes[i].set_title(gene, fontweight='bold', fontsize=10)
        axes[i].set_ylabel('Expression Level', fontsize=9)
        axes[i].grid(axis='y', alpha=0.3)
        
        control_data = expression_data.loc[labels == 'Control', gene]
        disease_data = expression_data.loc[labels != 'Control', gene]
        t_stat, p_val = stats.ttest_ind(disease_data, control_data)
        
        if p_val < 0.001:
            sig = '***'
        elif p_val < 0.01:
            sig = '**'
        elif p_val < 0.05:
            sig = '*'
        else:
            sig = 'ns'
        
        axes[i].text(0.5, 0.95, f'p-value: {sig}', 
                    transform=axes[i].transAxes, 
                    ha='center', va='top', fontsize=8)
    
    for i in range(n_genes, len(axes)):
        axes[i].axis('off')
    
    plt.suptitle('Core Gene Expression Validation', 
                fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout()
    
    if filename_prefix:
        plt.savefig(f"{filename_prefix}_expression.png", dpi=300, bbox_inches='tight')
        print(f"✓ Expression validation plot saved: {filename_prefix}_expression.png")
    
    print("\nCalculating ROC curves for core genes...")
    
    fig2, ax = plt.subplots(figsize=(10, 8))
    
    y_true = (labels != 'Control').astype(int)
    
    auc_scores = {}
    
    for gene in core_genes[:5]:
        if gene not in expression_data.columns:
            continue
        
        y_score = expression_data[gene].values
        
        fpr, tpr, thresholds = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)
        auc_scores[gene] = roc_auc
        
        ax.plot(fpr, tpr, linewidth=2, 
               label=f'{gene} (AUC = {roc_auc:.3f}')
    
    ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.5, label='Random (AUC = 0.500)')
    
    ax.set_xlabel('False Positive Rate', fontsize=12, fontweight='bold')
    ax.set_ylabel('True Positive Rate', fontsize=12, fontweight='bold')
    ax.set_title('ROC Curves for Core Gene Validation', 
                fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=10, loc='lower right', frameon=True)
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    
    if filename_prefix:
        plt.savefig(f"{filename_prefix}_roc.png", dpi=300, bbox_inches='tight')
        print(f"✓ ROC curve plot saved: {filename_prefix}_roc.png")
    
    print("\nAUC Scores for Core Genes:")
    for gene, score in sorted(auc_scores.items(), key=lambda x: x[1], reverse=True):
        print(f"  {gene}: {score:.3f}")
    
    return fig1, fig2, auc_scores
