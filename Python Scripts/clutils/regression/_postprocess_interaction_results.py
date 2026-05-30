import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def postprocess_interaction_results(
    results_raw,
    fdr_method="fdr_bh",
    pval_cols=None,
    required_effects=None,
    warn_col="warn_divzero",
    warn_msg_col="warn_msg",
    verbose=True,
):
    """Clean raw interaction-level regression results and add q-values."""

    results = results_raw.copy()

    if warn_col not in results.columns:
        results[warn_col] = False

    if warn_msg_col not in results.columns:
        results[warn_msg_col] = "N/A"

    results[warn_col] = results[warn_col].fillna(False).astype(bool)
    results[warn_msg_col] = results[warn_msg_col].fillna("N/A")

    n_raw = len(results)
    n_divzero = int(results[warn_col].sum())
    results = results.loc[~results[warn_col]].copy()

    if pval_cols is None:
        pval_cols = sorted(col for col in results.columns if col.startswith("pval_"))

    if required_effects is None:
        required_effects = [col.replace("pval_", "", 1) for col in pval_cols]

    unstable_mask = pd.Series(False, index=results.index)

    for effect in required_effects:
        pval_col = f"pval_{effect}"
        stderr_col = f"stderr_{effect}"

        if pval_col in results.columns:
            unstable_mask |= results[pval_col].isna()
        else:
            unstable_mask |= True

        if stderr_col in results.columns:
            unstable_mask |= results[stderr_col].isna()
            unstable_mask |= np.isinf(results[stderr_col])
        else:
            unstable_mask |= True

    results["is_unstable"] = unstable_mask
    n_unstable = int(results["is_unstable"].sum())
    results = results.loc[~results["is_unstable"]].copy()

    for pval_col in pval_cols:
        if pval_col not in results.columns:
            continue

        effect = pval_col.replace("pval_", "", 1)
        qval_col = f"qval_{effect}"
        valid_pvals = results[pval_col].notna()

        results[qval_col] = np.nan
        if valid_pvals.any():
            results.loc[valid_pvals, qval_col] = multipletests(
                results.loc[valid_pvals, pval_col],
                method=fdr_method,
            )[1]

    summary = {
        "n_raw": n_raw,
        "n_divzero": n_divzero,
        "n_after_divzero_filter": n_raw - n_divzero,
        "n_unstable": n_unstable,
        "n_final": len(results),
        "pval_cols_corrected": pval_cols,
        "required_effects": required_effects,
    }

    if verbose:
        print(f"Number of raw results: {n_raw}")
        print(f"Number of divide-by-zero conditions: {n_divzero}")
        print(f"Number of unstable conditions: {n_unstable}")
        print(f"Number of final results: {len(results)}")

    return results, summary
