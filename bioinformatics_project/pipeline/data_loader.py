def load_paper_geo_data(geo_cache_dir=None):
    try:
        from bioinformatics_project.geo_loader import (
            load_ad_training_data,
            load_ad_validation_data,
            load_oa_training_data,
            load_oa_validation_data,
            DEFAULT_CACHE_DIR,
        )
    except ImportError as e:
        raise ImportError(
            "Real GEO data requires the in-repo package 'bioinformatics_project.geo_loader'. "
            "Ensure the project root is on PYTHONPATH (we add it in run_pipeline.py) and that dependencies like GEOparse are installed."
        ) from e
    cache = geo_cache_dir or DEFAULT_CACHE_DIR
    expr_ad, labels_ad = load_ad_training_data(cache)
    expr_oa, labels_oa = load_oa_training_data(cache)
    expr_ad_val, labels_ad_val = load_ad_validation_data(cache)
    expr_oa_val, labels_oa_val = load_oa_validation_data(cache)
    return (
        expr_ad, labels_ad,
        expr_oa, labels_oa,
        expr_ad_val, labels_ad_val,
        expr_oa_val, labels_oa_val,
    )
