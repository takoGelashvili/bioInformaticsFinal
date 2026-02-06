import pandas as pd
import gseapy as gp


def run_enrichment_analysis(genes, cutoff=0.05):
    print(f"\n{'='*70}")
    print("GO/KEGG Enrichment Analysis (Enrichr)")
    print(f"{'='*70}")
    
    empty_df = lambda cat: pd.DataFrame(
        columns=['Term', 'Count', 'P_value', 'Category']
    )
    enrichment_results = {
        'GO_BP': empty_df('Biological Process'),
        'GO_MF': empty_df('Molecular Function'),
        'GO_CC': empty_df('Cellular Component'),
        'KEGG': empty_df('KEGG Pathway'),
    }
    
    if not genes or len(genes) == 0:
        print("⚠ No genes provided. Skipping enrichment.")
        return enrichment_results
    
    gene_list = [str(g).strip() for g in genes if pd.notna(g) and str(g).strip()]
    if not gene_list:
        print("⚠ No valid gene symbols. Skipping enrichment.")
        return enrichment_results
    
    gene_sets = [
        'GO_Biological_Process_2018',
        'GO_Molecular_Function_2018',
        'GO_Cellular_Component_2018',
        'KEGG_2021_Human',
    ]
    gene_set_to_key = {
        'GO_Biological_Process_2018': ('GO_BP', 'Biological Process'),
        'GO_Molecular_Function_2018': ('GO_MF', 'Molecular Function'),
        'GO_Cellular_Component_2018': ('GO_CC', 'Cellular Component'),
        'KEGG_2021_Human': ('KEGG', 'KEGG Pathway'),
        'KEGG_2016': ('KEGG', 'KEGG Pathway'),
    }
    
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism='human',
            outdir=None,
            cutoff=cutoff,
        )
    except Exception as e:
        print(f"⚠ Enrichr API error: {e}")
        return enrichment_results
    
    if not hasattr(enr, 'results') or enr.results is None or enr.results.empty:
        print("⚠ No enrichment results returned from Enrichr.")
        return enrichment_results
    
    df = enr.results
    pcol = 'P-value' if 'P-value' in df.columns else 'P_value'
    if pcol not in df.columns:
        print("⚠ Enrichment results missing P-value column.")
        return enrichment_results
    df = df.rename(columns={pcol: 'P_value'})
    
    if 'Overlap' in df.columns:
        df['Count'] = df['Overlap'].astype(str).str.split('/').str[0].fillna(0).astype(int)
    else:
        df['Count'] = 0
    
    if 'Gene_set' not in df.columns:
        print("⚠ Enrichment results missing Gene_set column.")
        return enrichment_results
    
    key_to_cat = {
        'GO_BP': 'Biological Process',
        'GO_MF': 'Molecular Function',
        'GO_CC': 'Cellular Component',
        'KEGG': 'KEGG Pathway',
    }
    for key, cat_name in key_to_cat.items():
        gs_names = [gs for gs, (k, _) in gene_set_to_key.items() if k == key]
        sub = df.loc[df['Gene_set'].isin(gs_names), ['Term', 'Count', 'P_value']].copy()
        sub['P_value'] = sub['P_value'].astype(float)
        sub['Category'] = cat_name
        enrichment_results[key] = sub.drop_duplicates(subset=['Term']).reset_index(drop=True)
    
    print("✓ Enrichment analysis completed (Enrichr)")
    for key, title in [('GO_BP', 'GO Biological Process'), ('GO_MF', 'GO Molecular Function'),
                      ('GO_CC', 'GO Cellular Component'), ('KEGG', 'KEGG Pathways')]:
        n = len(enrichment_results[key])
        print(f"  - {title}: {n} terms")
    
    return enrichment_results
