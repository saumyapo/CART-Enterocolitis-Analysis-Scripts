import pandas as pd
from tqdm.auto import tqdm


def run_interaction_ols_loop(
    data,
    fit_func,
    group_cols=None,
    metadata_cols=None,
    expected_n_groups=None,
    show_progress=True,
    **fit_kwargs,
):
    """Run interaction-level OLS fitting across grouped rows."""

    group_cols = group_cols or ["source", "target", "interaction_name"]
    metadata_cols = metadata_cols or [
        "ligand",
        "receptor",
        "pathway_name",
        "annotation",
    ]

    missing_group_cols = [col for col in group_cols if col not in data.columns]
    if missing_group_cols:
        raise ValueError(f"Missing group columns: {missing_group_cols}")

    available_metadata_cols = [col for col in metadata_cols if col in data.columns]
    grouped = data.groupby(group_cols, sort=False, observed=True)

    total = grouped.ngroups if expected_n_groups is None else expected_n_groups
    iterator = tqdm(grouped, total=total) if show_progress else grouped
    results = []

    for group_key, subdf in iterator:
        if not isinstance(group_key, tuple):
            group_key = (group_key,)

        fit_result = fit_func(subdf, **fit_kwargs)

        for col, value in zip(group_cols, group_key):
            fit_result[col] = value

        for col in available_metadata_cols:
            fit_result[col] = subdf[col].iloc[0]

        results.append(fit_result)

    return pd.DataFrame(results)
