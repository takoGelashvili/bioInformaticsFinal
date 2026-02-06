import warnings
from typing import Tuple, List, Optional, Dict
import pandas as pd

from .constants import (
    DEFAULT_CACHE_DIR, AD_TRAINING_IDS, AD_VALIDATION_IDS,
    OA_TRAINING_IDS, OA_VALIDATION_IDS,
    AD_TRAINING_MINIMAL, OA_TRAINING_MINIMAL
)
from .utils import get_ad_gse_ids, get_oa_gse_ids
from .geo_loader import load_single_gse, merge_and_batch_correct


def load_ad_training_data(
    cache_dir: str = DEFAULT_CACHE_DIR,
    ad_gse_ids: Optional[List[str]] = None,
    minimal: bool = False,
) -> Tuple[pd.DataFrame, pd.Series]:
    ids = ad_gse_ids if ad_gse_ids is not None else (AD_TRAINING_MINIMAL if minimal else get_ad_gse_ids())
    expr_list = []
    label_list = []
    for gse_id in ids:
        try:
            expr, labels, _ = load_single_gse(gse_id, cache_dir, "AD")
            if expr.shape[0] > 0 and expr.shape[1] > 0:
                expr_list.append(expr)
                label_list.append(labels)
        except Exception as e:
            warnings.warn(f"Failed to load {gse_id}: {e}")
    return merge_and_batch_correct(
        expr_list, label_list,
        batch_ids=ids[: len(expr_list)],
    )


def load_ad_validation_data(cache_dir: str = DEFAULT_CACHE_DIR) -> Tuple[pd.DataFrame, pd.Series]:
    expr, labels, _ = load_single_gse(AD_VALIDATION_IDS[0], cache_dir, "AD")
    return expr, labels


def load_oa_training_data(
    cache_dir: str = DEFAULT_CACHE_DIR,
    oa_gse_ids: Optional[List[str]] = None,
    minimal: bool = False,
) -> Tuple[pd.DataFrame, pd.Series]:
    ids = oa_gse_ids if oa_gse_ids is not None else (OA_TRAINING_MINIMAL if minimal else get_oa_gse_ids())
    expr_list = []
    label_list = []
    for gse_id in ids:
        try:
            expr, labels, _ = load_single_gse(gse_id, cache_dir, "OA")
            if expr.shape[0] > 0 and expr.shape[1] > 0:
                expr_list.append(expr)
                label_list.append(labels)
        except Exception as e:
            warnings.warn(f"Failed to load {gse_id}: {e}")
    return merge_and_batch_correct(
        expr_list, label_list,
        batch_ids=ids[: len(expr_list)],
    )


def load_oa_validation_data(cache_dir: str = DEFAULT_CACHE_DIR) -> Tuple[pd.DataFrame, pd.Series]:
    expr, labels, _ = load_single_gse(OA_VALIDATION_IDS[0], cache_dir, "OA")
    return expr, labels


def load_paper_data(
    cache_dir: Optional[str] = None,
    ad_gse_ids: Optional[List[str]] = None,
    oa_gse_ids: Optional[List[str]] = None,
    minimal: bool = False,
    skip_validation: bool = False,
) -> Dict[str, Tuple[pd.DataFrame, pd.Series]]:
    cache_dir = cache_dir or DEFAULT_CACHE_DIR
    out = {
        "ad_training": load_ad_training_data(cache_dir, ad_gse_ids=ad_gse_ids, minimal=minimal),
        "oa_training": load_oa_training_data(cache_dir, oa_gse_ids=oa_gse_ids, minimal=minimal),
    }
    if skip_validation:
        out["ad_validation"] = (pd.DataFrame(), pd.Series(dtype=object))
        out["oa_validation"] = (pd.DataFrame(), pd.Series(dtype=object))
    else:
        out["ad_validation"] = load_ad_validation_data(cache_dir)
        out["oa_validation"] = load_oa_validation_data(cache_dir)
    return out
