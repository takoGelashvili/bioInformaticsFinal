import os
import re
import pandas as pd
from typing import List, Optional


def parse_gse_list(s: str) -> List[str]:
    if not s or not str(s).strip():
        return []
    return [x.strip() for x in str(s).split(",") if x.strip()]


def get_ad_gse_ids() -> List[str]:
    ids = os.environ.get("AD_GSE_IDS")
    if ids:
        return parse_gse_list(ids)
    from .constants import GEO_IDS_FILE, AD_TRAINING_IDS
    path = os.path.join(os.path.dirname(os.path.dirname(__file__)), GEO_IDS_FILE)
    if os.path.isfile(path):
        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if line.upper().startswith("AD:"):
                    return parse_gse_list(line.split(":", 1)[1])
    return list(AD_TRAINING_IDS)


def get_oa_gse_ids() -> List[str]:
    ids = os.environ.get("OA_GSE_IDS")
    if ids:
        return parse_gse_list(ids)
    from .constants import GEO_IDS_FILE, OA_TRAINING_IDS
    path = os.path.join(os.path.dirname(os.path.dirname(__file__)), GEO_IDS_FILE)
    if os.path.isfile(path):
        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if line.upper().startswith("OA:"):
                    return parse_gse_list(line.split(":", 1)[1])
    return list(OA_TRAINING_IDS)


def infer_value_column(table: pd.DataFrame) -> str:
    for name in ["VALUE", "Signal", "CH1_INTENSITY", "RPKM", "TPM", "raw_value", "norm_value"]:
        if name in table.columns:
            return name
    for c in table.columns:
        if c.upper() in ("ID_REF", "ID", "PROBEID", "REF"):
            continue
        if pd.api.types.is_numeric_dtype(table[c]):
            return c
    return table.columns[1]


def infer_id_column(table: pd.DataFrame) -> str:
    for name in ["ID_REF", "ID", "PROBEID", "REF", "ProbeID"]:
        if name in table.columns:
            return name
    return table.columns[0]


def infer_gene_symbol_column(gpl_table: pd.DataFrame) -> Optional[str]:
    candidates = [
        "Gene symbol", "Gene Symbol", "GENE_SYMBOL", "Gene_Symbol",
        "GeneSymbol", "gene_assignment", "GENE_NAME", "Gene name",
        "Symbol", "SYMBOL", "Gene", "GENE", "ILMN_Gene", "GeneName",
    ]
    for c in candidates:
        if c in gpl_table.columns:
            return c
    for col in gpl_table.columns:
        if "symbol" in col.lower() or "gene" in col.lower():
            return col
    return None


def extract_gene_from_assignment(assignment: str) -> str:
    if pd.isna(assignment) or not isinstance(assignment, str):
        return ""
    parts = re.split(r"\s*///\s*|\s*//\s*|\t", str(assignment).strip())
    for p in parts:
        p = p.strip()
        if p and not p.startswith("NM_") and not p.startswith("NR_") and len(p) < 20:
            return p
    return ""


def collapse_probes_to_genes(expr: pd.DataFrame, probe_to_gene: pd.Series) -> pd.DataFrame:
    def _map_probe(p):
        if p in probe_to_gene.index:
            g = probe_to_gene.loc[p]
            return str(g).strip() if pd.notna(g) and str(g).strip() else ""
        return ""
    gene_map = expr.columns.map(_map_probe)
    expr = expr.copy()
    expr.columns = gene_map
    expr = expr.loc[:, (expr.columns != "") & (~expr.columns.isna())]
    expr = expr.T.groupby(expr.columns).mean().T
    return expr


def get_sample_labels(gsm_metadata: dict, gsm_name: str, disease_label: str) -> str:
    for key in ["characteristics_ch1", "title", "source_name_ch1"]:
        if key not in gsm_metadata:
            continue
        vals = gsm_metadata[key]
        if not isinstance(vals, list):
            vals = [vals]
        text = " ".join(str(v).lower() for v in vals)
        if "control" in text or "normal" in text or "non-demented" in text or "healthy" in text:
            return "Control"
        if "alzheimer" in text or "ad " in text or "demented" in text or "disease" in text:
            return disease_label
        if "osteoarthritis" in text or " oa " in text or "oa patient" in text:
            return disease_label
        if "normal knee" in text or "non-oa" in text:
            return "Control"
    title = " ".join(gsm_metadata.get("title", []))
    if "control" in title.lower():
        return "Control"
    return disease_label


def soft_file_path(gse_id: str, destdir: str) -> Optional[str]:
    if not os.path.isdir(destdir):
        return None
    for name in (
        f"{gse_id}.soft.gz",
        f"{gse_id}.soft",
        f"{gse_id}_family.soft.gz",
        f"{gse_id}_family.soft",
    ):
        path = os.path.join(destdir, name)
        if os.path.isfile(path):
            return path
    for name in os.listdir(destdir):
        if name.startswith(gse_id) and ".soft" in name and os.path.isfile(os.path.join(destdir, name)):
            return os.path.join(destdir, name)
    return None
