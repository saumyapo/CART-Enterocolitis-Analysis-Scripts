
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from ._regression_frameworks import MODEL_FRAMEWORKS


def fit_one_model(
    subdf,
    framework="three_group",
    prob_col="logprob",
    condition_col="condition",
    sample_col="sample",
    chemistry_col="chemistry",
    include_chemistry=True,
    hc_label="NV",
    tw_label="CAR-TCTRL",
    ts_label="CAR-TEC",
    cov_type="HC3",
    extra_covariates=None,
    return_model=False,
):
    if framework not in MODEL_FRAMEWORKS:
        raise ValueError(
            f"Unknown framework '{framework}'. Use one of: {list(MODEL_FRAMEWORKS)}"
        )

    config = MODEL_FRAMEWORKS[framework]
    data = subdf.copy()
    extra_covariates = list(extra_covariates or [])

    label_map = {
        "NV": hc_label,
        "CAR-TCTRL": tw_label,
        "CAR-TEC": ts_label,
    }

    conditions = [label_map[label] for label in config["conditions"]]
    reference_label = label_map[config["reference_label"]]
    positive_label = label_map[config["positive_label"]]

    required_cols = [prob_col, condition_col, sample_col]
    if include_chemistry:
        required_cols.append(chemistry_col)
    required_cols.extend(extra_covariates)

    missing_cols = [col for col in required_cols if col not in data.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    data[condition_col] = data[condition_col].astype(str).str.strip()
    data[sample_col] = data[sample_col].astype(str).str.strip()
    data[prob_col] = np.asarray(data[prob_col], dtype=float)
    data = data.loc[data[condition_col].isin(conditions)].copy()

    if include_chemistry:
        data[chemistry_col] = data[chemistry_col].astype(str).str.strip()

    for covariate in extra_covariates:
        data[covariate] = pd.to_numeric(data[covariate], errors="coerce")

    data = data.dropna(subset=required_cols)

    if "is_treated" in config["terms"]:
        data["is_treated"] = (data[condition_col] != hc_label).astype(int)

    if "has_colitis" in config["terms"]:
        data["has_colitis"] = (data[condition_col] == ts_label).astype(int)

    terms = list(config["terms"])
    if include_chemistry:
        terms.append(f"C({chemistry_col})")
    terms.extend(extra_covariates)

    formula = f"{prob_col} ~ {' + '.join(terms)}" if terms else f"{prob_col} ~ 1"

    with warnings.catch_warnings(record=True) as caught_warnings:
        warnings.simplefilter("always")
        model = smf.ols(formula, data=data).fit(cov_type=cov_type)

    warning_messages = [str(warning.message) for warning in caught_warnings]
    divzero_messages = [
        message for message in warning_messages if "divide by zero" in message.lower()
    ]

    result = {
        "framework": framework,
        "formula": formula,
        "conditions": ",".join(conditions),
        "reference_label": reference_label,
        "positive_label": positive_label,
        "n_rows": len(data),
        "n_samples": data[sample_col].nunique(),
        "coef_intercept": model.params.get("Intercept", np.nan),
        "r2": model.rsquared,
        "adj_r2": model.rsquared_adj,
        "aic": model.aic,
        "bic": model.bic,
        "warn_divzero": bool(divzero_messages),
        "warn_msg": "; ".join(divzero_messages) if divzero_messages else None,
    }

    for condition in conditions:
        safe_condition = condition.replace("-", "_").replace(" ", "_").replace("+", "_")
        result[f"n_condition_{safe_condition}"] = int((data[condition_col] == condition).sum())

    for term in [*config["terms"], *extra_covariates]:
        result[f"coef_{term}"] = model.params.get(term, np.nan)
        result[f"pval_{term}"] = model.pvalues.get(term, np.nan)
        result[f"stderr_{term}"] = model.bse.get(term, np.nan)
        result[f"t_{term}"] = model.tvalues.get(term, np.nan)

    if return_model:
        result["model"] = model

    return result
