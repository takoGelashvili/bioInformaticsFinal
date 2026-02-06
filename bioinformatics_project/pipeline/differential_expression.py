import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests


def differential_expression_analysis(expression_data, labels, 
                                    pvalue_threshold=0.05, 
                                    logfc_threshold=0.5):
    print(f"\n{'='*70}")
    print("Differential Expression Analysis")
    print(f"{'='*70}")
    
    control_mask = labels == 'Control'
    disease_conditions = labels[labels != 'Control'].unique()
    
    if len(disease_conditions) == 0:
        raise ValueError("No disease samples found!")
    
    disease_name = disease_conditions[0]
    disease_mask = labels == disease_name
    
    results = []
    
    for gene in expression_data.columns:
        control_expr = expression_data.loc[control_mask, gene]
        disease_expr = expression_data.loc[disease_mask, gene]
        
        mean_control = control_expr.mean()
        mean_disease = disease_expr.mean()
        std_control = control_expr.std()
        std_disease = disease_expr.std()
        
        log2fc = np.log2((mean_disease + 0.001) / (mean_control + 0.001))
        
        t_stat, p_value = stats.ttest_ind(disease_expr, control_expr)
        
        results.append({
            'Gene': gene,
            'Mean_Control': mean_control,
            'Mean_Disease': mean_disease,
            'Std_Control': std_control,
            'Std_Disease': std_disease,
            'Log2FC': log2fc,
            'T_Statistic': t_stat,
            'P_Value': p_value
        })
    
    df_results = pd.DataFrame(results)
    
    _, adj_pvals, _, _ = multipletests(df_results['P_Value'], 
                                       method='fdr_bh')
    df_results['Adj_P_Value'] = adj_pvals
    
    df_results['Regulation'] = 'Not Significant'
    
    upregulated = (df_results['Adj_P_Value'] < pvalue_threshold) & \
                 (df_results['Log2FC'] > logfc_threshold)
    downregulated = (df_results['Adj_P_Value'] < pvalue_threshold) & \
                   (df_results['Log2FC'] < -logfc_threshold)
    
    df_results.loc[upregulated, 'Regulation'] = 'Upregulated'
    df_results.loc[downregulated, 'Regulation'] = 'Downregulated'
    
    n_up = (df_results['Regulation'] == 'Upregulated').sum()
    n_down = (df_results['Regulation'] == 'Downregulated').sum()
    
    print(f"✓ Analysis complete for {disease_name}")
    print(f"  - Upregulated genes: {n_up}")
    print(f"  - Downregulated genes: {n_down}")
    print(f"  - Total DEGs: {n_up + n_down}")
    
    return df_results


def find_common_degs(degs1, degs2, disease1_name="AD", disease2_name="OA"):
    print(f"\n{'='*70}")
    print(f"Finding Common DEGs between {disease1_name} and {disease2_name}")
    print(f"{'='*70}")
    
    degs1_up = set(degs1[degs1['Regulation'] == 'Upregulated']['Gene'])
    degs1_down = set(degs1[degs1['Regulation'] == 'Downregulated']['Gene'])
    
    degs2_up = set(degs2[degs2['Regulation'] == 'Upregulated']['Gene'])
    degs2_down = set(degs2[degs2['Regulation'] == 'Downregulated']['Gene'])
    
    co_degs_up = degs1_up.intersection(degs2_up)
    co_degs_down = degs1_down.intersection(degs2_down)
    
    print(f"✓ Common upregulated DEGs: {len(co_degs_up)}")
    print(f"✓ Common downregulated DEGs: {len(co_degs_down)}")
    print(f"✓ Total Co-DEGs: {len(co_degs_up) + len(co_degs_down)}")
    
    from matplotlib_venn import venn2
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    venn2([degs1_up, degs2_up], 
          set_labels=(disease1_name, disease2_name), ax=ax1)
    ax1.set_title(f'Upregulated DEGs\n(Common: {len(co_degs_up)})', 
                 fontsize=12, fontweight='bold')
    
    venn2([degs1_down, degs2_down], 
          set_labels=(disease1_name, disease2_name), ax=ax2)
    ax2.set_title(f'Downregulated DEGs\n(Common: {len(co_degs_down)})', 
                 fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    
    co_degs = {
        'upregulated': list(co_degs_up),
        'downregulated': list(co_degs_down),
        'all': list(co_degs_up) + list(co_degs_down)
    }
    
    return fig, co_degs
