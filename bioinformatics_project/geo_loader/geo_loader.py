import os
import numpy as np
import pandas as pd
import warnings
from typing import Tuple, Optional

try:
    import GEOparse
    HAS_GEOPARSE = True
except ImportError:
    HAS_GEOPARSE = False

from .utils import (
    infer_value_column, infer_id_column, infer_gene_symbol_column,
    extract_gene_from_assignment, collapse_probes_to_genes,
    get_sample_labels, soft_file_path
)


def load_single_gse(
    gse_id: str,
    destdir: str,
    disease_label: str,
    value_col: Optional[str] = None,
    id_col: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.Series, Optional[pd.Series]]:
    if not HAS_GEOPARSE:
        raise ImportError("GEOparse is required. Install with: pip install GEOparse")

    os.makedirs(destdir, exist_ok=True)
    path = soft_file_path(gse_id, destdir)
    if path:
        gse = GEOparse.get_GEO(filepath=path, silent=True)
    else:
        gse = GEOparse.get_GEO(geo=gse_id, destdir=destdir, silent=True)

    probe_to_gene = None
    if gse.gpls:
        gpl = list(gse.gpls.values())[0]
        gpl_table = gpl.table
        id_col_pl = infer_id_column(gpl_table)
        gene_col = infer_gene_symbol_column(gpl_table)

        if gene_col:
            probe_to_gene = gpl_table.set_index(id_col_pl)[gene_col].astype(str)
        else:
            if "gene_assignment" in gpl_table.columns:
                gene_syms = gpl_table[id_col_pl].to_frame()
                gene_syms["Gene"] = gpl_table["gene_assignment"].map(extract_gene_from_assignment)
                probe_to_gene = gene_syms.set_index(id_col_pl)["Gene"]

    id_col_s = None
    val_col_s = None
    expr_list = []
    labels_list = []

    for gsm_name, gsm in gse.gsms.items():
        if gsm.table is None or gsm.table.empty:
            continue
        tbl = gsm.table
        if id_col_s is None:
            id_col_s = id_col or infer_id_column(tbl)
            val_col_s = value_col or infer_value_column(tbl)
        if id_col_s not in tbl.columns or val_col_s not in tbl.columns:
            continue
        expr_series = tbl.set_index(id_col_s)[val_col_s]
        expr_series = pd.to_numeric(expr_series, errors="coerce").dropna()
        expr_list.append(expr_series)
        label = get_sample_labels(gsm.metadata, gsm_name, disease_label)
        labels_list.append(label)

    if not expr_list:
        return pd.DataFrame(), pd.Series(dtype=object), probe_to_gene

    expr = pd.DataFrame(expr_list, index=[f"{gse_id}_{i}" for i in range(len(expr_list))])
    labels = pd.Series(labels_list, index=expr.index, name="Condition")

    if probe_to_gene is not None and not probe_to_gene.empty:
        expr = collapse_probes_to_genes(expr, probe_to_gene)

    return expr, labels, probe_to_gene


def merge_and_batch_correct(
    expr_list: list,
    label_list: list,
    batch_ids: list,
) -> Tuple[pd.DataFrame, pd.Series]:
    if not expr_list:
        return pd.DataFrame(), pd.Series(dtype=object)

    common_genes = expr_list[0].columns.tolist()
    for e in expr_list[1:]:
        common_genes = list(set(common_genes) & set(e.columns.tolist()))

    merged = []
    merged_labels = []
    for i, (expr, labels) in enumerate(zip(expr_list, label_list)):
        keep = [g for g in common_genes if g in expr.columns]
        if not keep:
            continue
        merged.append(expr[keep])
        merged_labels.append(labels)

    expr_merged = pd.concat(merged, axis=0)
    labels_merged = pd.concat(merged_labels, axis=0)

    if len(expr_list) > 1 and len(common_genes) > 0:
        try:
            from combat.py_combat import py_combat
            batch_vec = []
            for j, (expr, _) in enumerate(zip(expr_list, label_list)):
                n = expr[ [g for g in common_genes if g in expr.columns] ].shape[0]
                batch_vec.extend([j] * n)
            if len(batch_vec) == expr_merged.shape[0]:
                dat = expr_merged.values.T
                expr_corrected = py_combat(dat, np.array(batch_vec))
                expr_merged = pd.DataFrame(
                    expr_corrected.T,
                    index=expr_merged.index,
                    columns=expr_merged.columns
                )
        except Exception as e:
            warnings.warn(f"Batch correction skipped: {e}")

    return expr_merged, labels_merged
