import pandas as pd


def prepare_interaction_regression_df(
    df,
    hc_label="NV",
    tw_label="CAR-TCTRL",
    ts_label="CAR-TEC",
    key_cols=None,
    response_col="logprob",
    sample_col="sample",
    condition_col="condition",
    chemistry_col="chemistry",
    include_ts=True,
    center_by_sample=False,
    sample_center_method="mean",
    extra_covariates=None,
):
    """
    Prepare post-filtered interaction data for regression modeling.

    Adds:
        is_treated : 1 if condition is not HC, else 0
        has_colitis : 1 if condition is TS, else 0

    Optionally recenters the response within each sample using the requested
    summary statistic.
    """

    key_cols = key_cols or ["source", "target", "interaction_name"]
    extra_covariates = list(extra_covariates or [])

    conditions = [hc_label, tw_label]
    if include_ts:
        conditions.append(ts_label)

    required_cols = key_cols + [
        response_col,
        sample_col,
        condition_col,
        chemistry_col,
    ] + extra_covariates

    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    data = df.loc[df[condition_col].astype(str).str.strip().isin(conditions)].copy()
    data[response_col] = pd.to_numeric(data[response_col], errors="coerce")
    data[sample_col] = data[sample_col].astype(str).str.strip()
    data[condition_col] = data[condition_col].astype(str).str.strip()
    data[chemistry_col] = data[chemistry_col].astype(str).str.strip()

    for covariate in extra_covariates:
        data[covariate] = pd.to_numeric(data[covariate], errors="coerce")

    data = data.dropna(subset=required_cols)

    if center_by_sample:
        if sample_center_method not in {"mean", "median"}:
            raise ValueError("sample_center_method must be either 'mean' or 'median'.")

        sample_center = data.groupby(sample_col, observed=True)[response_col].transform(
            sample_center_method
        )
        data[response_col] = data[response_col] - sample_center

    data["is_treated"] = (data[condition_col] != hc_label).astype(int)
    data["has_colitis"] = 0

    if include_ts:
        data["has_colitis"] = (data[condition_col] == ts_label).astype(int)

    return data
