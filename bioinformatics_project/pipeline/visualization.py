import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib.patches import Patch


def perform_pca(expression_data, labels, title="PCA Analysis", filename=None, random_state=42):
    print(f"\nPerforming PCA Analysis...")
    
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(expression_data)
    
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_data)
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    colors = {'Control': '#3498db', 'AD': '#e74c3c', 'OA': '#2ecc71'}
    
    for condition in labels.unique():
        mask = labels == condition
        color = colors.get(condition, '#95a5a6')
        ax.scatter(pca_result[mask, 0], pca_result[mask, 1], 
                  label=condition, alpha=0.7, s=120, 
                  edgecolors='black', linewidth=0.5, color=color)
    
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)', 
                 fontsize=12, fontweight='bold')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)', 
                 fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"✓ PCA plot saved: {filename}")
    
    return fig, pca


def plot_volcano(deg_results, title="Volcano Plot", filename=None):
    print(f"\nCreating volcano plot...")
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    x = deg_results['Log2FC']
    y = -np.log10(deg_results['Adj_P_Value'] + 1e-300)
    colors = deg_results['Regulation'].map({
        'Upregulated': '#e74c3c',
        'Downregulated': '#3498db',
        'Not Significant': '#95a5a6'
    })
    
    scatter = ax.scatter(x, y, c=colors, alpha=0.6, s=20, edgecolors='none')
    
    ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', 
              linewidth=1, alpha=0.5, label='P = 0.05')
    ax.axvline(x=0.5, color='black', linestyle='--', 
              linewidth=1, alpha=0.5)
    ax.axvline(x=-0.5, color='black', linestyle='--', 
              linewidth=1, alpha=0.5, label='|Log2FC| = 0.5')
    
    ax.set_xlabel('Log2 Fold Change', fontsize=12, fontweight='bold')
    ax.set_ylabel('-Log10 Adjusted P-value', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    
    legend_elements = [
        Patch(facecolor='#e74c3c', label='Upregulated'),
        Patch(facecolor='#3498db', label='Downregulated'),
        Patch(facecolor='#95a5a6', label='Not Significant')
    ]
    ax.legend(handles=legend_elements, fontsize=10, frameon=True)
    
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"✓ Volcano plot saved: {filename}")
    
    return fig


def plot_ppi_network(G, core_genes, filename=None, random_state=42):
    print(f"\nVisualizing PPI network...")
    
    fig, ax = plt.subplots(figsize=(14, 14))
    
    pos = nx.spring_layout(G, k=0.5, iterations=50, seed=random_state)
    
    node_colors = ['#e74c3c' if node in core_genes else '#3498db' 
                  for node in G.nodes()]
    
    node_sizes = [300 + 100 * G.degree(node) for node in G.nodes()]
    
    nx.draw_networkx_edges(G, pos, alpha=0.2, width=0.5, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, 
                          node_size=node_sizes, alpha=0.8, ax=ax,
                          edgecolors='black', linewidths=1)
    
    core_labels = {node: node for node in core_genes}
    nx.draw_networkx_labels(G, pos, labels=core_labels, 
                           font_size=9, font_weight='bold', ax=ax)
    
    ax.set_title('Protein-Protein Interaction Network\n(Red nodes = Core genes)', 
                fontsize=14, fontweight='bold', pad=20)
    ax.axis('off')
    
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"✓ PPI network plot saved: {filename}")
    
    return fig


def plot_enrichment(enrichment_results, top_n=10, filename=None):
    print(f"\nPlotting enrichment results...")
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    axes = [ax1, ax2, ax3, ax4]
    categories = ['GO_BP', 'GO_MF', 'GO_CC', 'KEGG']
    titles = ['GO: Biological Process', 'GO: Molecular Function', 
             'GO: Cellular Component', 'KEGG Pathways']
    
    for ax, cat, title in zip(axes, categories, titles):
        df = enrichment_results[cat].sort_values('P_value').head(top_n)
        
        if df.empty or len(df) == 0:
            ax.text(0.5, 0.5, 'No significant terms', ha='center', va='center',
                    fontsize=11, transform=ax.transAxes)
            ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
            ax.set_xlabel('Gene Count', fontsize=11, fontweight='bold')
            continue
        
        y_pos = np.arange(len(df))
        p_vals = df['P_value'].replace(0, np.nextafter(0, 1))
        colors = plt.cm.RdYlBu_r(p_vals / p_vals.max())
        
        ax.barh(y_pos, df['Count'], color=colors, edgecolor='black', linewidth=0.5)
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df['Term'], fontsize=9)
        ax.set_xlabel('Gene Count', fontsize=11, fontweight='bold')
        ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=0.3)
        
        sm = plt.cm.ScalarMappable(cmap=plt.cm.RdYlBu_r,
                                   norm=plt.Normalize(vmin=df['P_value'].min(),
                                                    vmax=df['P_value'].max()))
        sm.set_array([])
    
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"✓ Enrichment plot saved: {filename}")
    
    return fig
