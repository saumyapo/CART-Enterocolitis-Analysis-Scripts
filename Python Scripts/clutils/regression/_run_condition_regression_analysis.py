from pathlib import Path

import pandas as pd

from ._filter_valid_interactions import filter_valid_interactions
from ._fit_one_model import fit_one_model
from ._postprocess_interaction_results import postprocess_interaction_results
from ._prepare_interaction_regression_df import prepare_interaction_regression_df
from ._regression_frameworks import MODEL_FRAMEWORKS
from ._run_interaction_ols_loop import run_interaction_ols_loop


def run_condition_regression_analysis(
    df=None,
    input_path=None,
    output_path=None,
    framework="nv_vs_cartctrl",
    condition_min_samples=None,
    min_total=6,
    key_cols=None,
    max_valid_permutations=None,
    required_effects=None,
    include_ts=None,
    center_by_sample=False,
    sample_center_method="mean",
    verbose=True,
    **fit_kwargs,
):
    """Run the full interaction regression pipeline for one framework."""

    key_cols = key_cols or ["source", "target", "interaction_name"]

    if framework not in MODEL_FRAMEWORKS:
        raise ValueError(
            f"Unknown framework '{framework}'. Use one of: {list(MODEL_FRAMEWORKS)}"
        )

    if df is None:
        if input_path is None:
            raise ValueError("Provide either df or input_path.")
        df = pd.read_parquet(input_path)
    else:
        df = df.copy()

    framework_conditions = MODEL_FRAMEWORKS[framework]["conditions"]

    if condition_min_samples is None:
        condition_min_samples = {condition: 3 for condition in framework_conditions}

    if include_ts is None:
        include_ts = "CAR-TEC" in framework_conditions

    df_valid, n_valid = filter_valid_interactions(
        df,
        condition_min_samples=condition_min_samples,
        min_total=min_total,
        key_cols=key_cols,
        verbose=verbose,
    )

    if max_valid_permutations is not None:
        if max_valid_permutations < 0:
            raise ValueError("max_valid_permutations must be non-negative or None.")

        valid_keys = df_valid[key_cols].drop_duplicates().head(max_valid_permutations)
        df_valid = df_valid.merge(valid_keys, on=key_cols, how="inner")
        n_valid = len(valid_keys)

        if verbose:
            print(f"Running first {n_valid} valid permutations.")

    df_valid_formatted = prepare_interaction_regression_df(
        df_valid,
        key_cols=key_cols,
        include_ts=include_ts,
        center_by_sample=center_by_sample,
        sample_center_method=sample_center_method,
        extra_covariates=fit_kwargs.get("extra_covariates"),
    )

    results_raw = run_interaction_ols_loop(
        df_valid_formatted,
        fit_func=fit_one_model,
        framework=framework,
        expected_n_groups=n_valid,
        **fit_kwargs,
    )

    if required_effects is None:
        required_effects = list(MODEL_FRAMEWORKS[framework]["terms"])

    results_pp, summary = postprocess_interaction_results(
        results_raw,
        required_effects=required_effects,
        verbose=verbose,
    )

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        results_pp.to_parquet(output_path)

    return {
        "df_valid": df_valid,
        "df_valid_formatted": df_valid_formatted,
        "results_raw": results_raw,
        "results_pp": results_pp,
        "summary": summary,
        "n_valid": n_valid,
    }


def run_default_condition_regression_suite(
    df=None,
    input_path=None,
    output_dir=None,
    framework_options=None,
    center_by_sample=False,
    sample_center_method="mean",
    max_valid_permutations=None,
    required_effects=None,
    include_ts=None,
    verbose=True,
    **fit_kwargs,
):
    """Run the default set of condition regression analyses and optionally save each output."""

    framework_options = framework_options or {
        "nv_vs_cartctrl": {
            "condition_min_samples": {"NV": 3, "CAR-TCTRL": 3},
            "min_total": 6,
            "output_name": "results_regr_NV_vs_CARTCTRL_pp.parquet",
        },
        "nv_vs_cartec": {
            "condition_min_samples": {"NV": 3, "CAR-TEC": 3},
            "min_total": 6,
            "output_name": "results_regr_NV_vs_CARTEC_pp.parquet",
        },
        "three_group": {
            "condition_min_samples": {"NV": 3, "CAR-TCTRL": 3, "CAR-TEC": 3},
            "min_total": 9,
            "output_name": "results_regr_NV_vs_CARTCTRL_vs_CARTEC_pp.parquet",
        },
        "cartctrl_vs_cartec": {
            "condition_min_samples": {"CAR-TCTRL": 3, "CAR-TEC": 3},
            "min_total": 6,
            "output_name": "results_regr_CARTCTRL_vs_CARTEC_pp.parquet",
        },
    }

    suite_results = {}

    for framework, options in framework_options.items():
        framework_output_path = None
        if output_dir is not None:
            framework_output_path = Path(output_dir) / options["output_name"]

        framework_fit_kwargs = {
            **fit_kwargs,
            **options.get("fit_kwargs", {}),
        }

        if "include_chemistry" in options:
            framework_fit_kwargs["include_chemistry"] = options["include_chemistry"]

        suite_results[framework] = run_condition_regression_analysis(
            df=df,
            input_path=input_path,
            output_path=framework_output_path,
            framework=framework,
            condition_min_samples=options.get("condition_min_samples"),
            min_total=options.get("min_total", 6),
            max_valid_permutations=options.get(
                "max_valid_permutations", max_valid_permutations
            ),
            required_effects=options.get("required_effects", required_effects),
            include_ts=options.get("include_ts", include_ts),
            center_by_sample=center_by_sample,
            sample_center_method=sample_center_method,
            verbose=verbose,
            **framework_fit_kwargs,
        )

    return suite_results
