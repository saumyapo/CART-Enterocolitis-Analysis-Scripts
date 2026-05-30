import pandas as pd


def filter_valid_interactions(
    df,
    condition_min_samples=None,
    min_total=6,
    key_cols=None,
    condition_col="condition",
    sample_col="sample",
    verbose=True,
):
    """
    Filter to source/target/interaction combinations with enough sample support
    in selected conditions.
    """

    key_cols = key_cols or ["source", "target", "interaction_name"]

    if condition_min_samples is None:
        condition_min_samples = {
            "NV": 2,
            "CAR-TCTRL": 2,
            "CAR-TEC": 2,
        }

    conditions = list(condition_min_samples)

    data = df.loc[df[condition_col].astype(str).str.strip().isin(conditions)].copy()
    data[condition_col] = data[condition_col].astype(str).str.strip()
    data[sample_col] = data[sample_col].astype(str).str.strip()

    counts = (
        data[key_cols + [sample_col, condition_col]]
        .drop_duplicates()
        .groupby(key_cols + [condition_col], observed=True)
        .size()
        .unstack(condition_col, fill_value=0)
        .reindex(columns=conditions, fill_value=0)
    )

    valid_mask = pd.Series(True, index=counts.index)

    for condition, min_samples in condition_min_samples.items():
        if min_samples is not None:
            valid_mask &= counts[condition] >= min_samples

    if min_total is not None:
        valid_mask &= counts.sum(axis=1) >= min_total

    valid_keys = counts.loc[valid_mask].reset_index()[key_cols]

    if verbose:
        print(f"Valid permutations: {len(valid_keys)}")

    filtered = data.merge(valid_keys, on=key_cols, how="inner")

    return filtered, len(valid_keys)
