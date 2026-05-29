from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional, Union

import numpy as np
import pandas as pd
from scipy.stats import beta
from statsmodels.stats.multitest import multipletests
from tqdm.auto import tqdm


PathLike = Union[str, Path]

def _rra_score(ranks: Iterable[float], n_lists_total: int) -> float:
    """Compute an RRA score from normalized ranks.

    Args:
        ranks: Iterable of normalized ranks where smaller values are better.
            Missing observations are encoded as ``np.nan``.
        n_lists_total: Total number of ranked lists or samples.

    Returns:
        RRA score. Smaller values indicate more consistently high ranks across samples.
    """

    rank_array = np.asarray(list(ranks), dtype=float)
    rank_array = rank_array[~np.isnan(rank_array)]

    if len(rank_array) == 0:
        return np.nan

    eps = 1e-12
    rank_array = np.clip(rank_array, eps, 1 - eps)
    rank_array = np.sort(rank_array)

    pvals = []
    for k, rank_value in enumerate(rank_array, start=1):
        pvals.append(beta.cdf(rank_value, k, n_lists_total + 1 - k))

    rra_pval = np.min(pvals)
    return min(1.0, rra_pval * len(rank_array))


def build_robust_rank_aggregation(
    source_file: PathLike,
    output_file: Optional[PathLike] = None,
    condition: Optional[str] = "CAR-TEC",
    show_progress: bool = True,
) -> pd.DataFrame:
    """Aggregate interaction ranks across samples using RRA.

    Args:
        source_file: Path to merged CellChat results file.
        output_file: Optional destination parquet path for the result table.
        condition: Optional condition label used to subset the merged table
            before ranking. When ``None``, all rows are included.
        show_progress: Whether to show a tqdm progress bar while computing RRA
            scores.

    Returns:
        Dataframe of interaction-level RRA scores and descriptive summaries.

    Raises:
        FileNotFoundError: If the source parquet does not exist.
        ValueError: If required columns are missing or filtering removes all
            rows.
    """
    source_path = Path(source_file)
    output_path = Path(output_file) if output_file is not None else None

    if not source_path.exists():
        raise FileNotFoundError(f"Source parquet does not exist: {source_path}")

    df = pd.read_parquet(source_path)

    required_columns = {
        "sample",
        "source",
        "target",
        "interaction_name",
        "logprob",
    }
    if condition is not None:
        required_columns.add("condition")

    missing_columns = sorted(required_columns.difference(df.columns))
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    if condition is not None:
        df = df[df["condition"].eq(condition)].copy()
    else:
        df = df.copy()

    if df.empty:
        raise ValueError("No rows remain after filtering the input dataframe.")

    df["interaction_tuple"] = (
        df["source"].astype(str)
        + "<->"
        + df["target"].astype(str)
        + ": "
        + df["interaction_name"].astype(str)
    )

    df["rank_pct"] = df.groupby("sample")["logprob"].rank(pct=True)
    df["rra_rank"] = 1 - df["rank_pct"]

    rank_mat = df.pivot_table(
        index="sample",
        columns="interaction_tuple",
        values="rra_rank",
        aggfunc="mean",
    )

    if rank_mat.empty:
        raise ValueError("No interaction ranks were generated from the input data.")

    interaction_tuples = rank_mat.columns.to_numpy()
    rank_array = rank_mat.to_numpy(dtype=float)
    n_samples_total = rank_array.shape[0]

    iterator = range(rank_array.shape[1])
    if show_progress:
        iterator = tqdm(iterator, desc="RRA scores")

    rra_pvals = np.array(
        [_rra_score(rank_array[:, column_index], n_samples_total) for column_index in iterator]
    )

    valid = ~np.isnan(rra_pvals)
    rra_fdr = np.full_like(rra_pvals, np.nan, dtype=float)
    if valid.any():
        rra_fdr[valid] = multipletests(rra_pvals[valid], method="fdr_bh")[1]

    rank_pct_mat = 1 - rank_array

    results = pd.DataFrame(
        {
            "interaction_tuple": interaction_tuples,
            "rra_pval": rra_pvals,
            "rra_fdr": rra_fdr,
            "median_rank_pct": np.nanmedian(rank_pct_mat, axis=0),
            "mean_rank_pct": np.nanmean(rank_pct_mat, axis=0),
            "sd_rank_pct": np.nanstd(rank_pct_mat, axis=0),
            "n_samples_observed": np.sum(~np.isnan(rank_pct_mat), axis=0),
            "top1_freq": np.nanmean(rank_pct_mat >= 0.99, axis=0),
            "top5_freq": np.nanmean(rank_pct_mat >= 0.95, axis=0),
            "top10_freq": np.nanmean(rank_pct_mat >= 0.90, axis=0),
        }
    )

    logprob_mat = df.pivot_table(
        index="sample",
        columns="interaction_tuple",
        values="logprob",
        aggfunc="mean",
    ).reindex(columns=interaction_tuples)

    results["median_logprob"] = logprob_mat.median(axis=0).to_numpy()
    results["mean_logprob"] = logprob_mat.mean(axis=0).to_numpy()
    results["nonzero_freq"] = (logprob_mat > 0).mean(axis=0).to_numpy()

    interaction_metadata = (
        df[["interaction_tuple", "source", "target", "interaction_name"]]
        .drop_duplicates(subset=["interaction_tuple"])
        .set_index("interaction_tuple")
    )

    results.insert(0, "interaction_name", interaction_metadata["interaction_name"].reindex(results["interaction_tuple"]).to_numpy())
    results.insert(0, "target", interaction_metadata["target"].reindex(results["interaction_tuple"]).to_numpy())
    results.insert(0, "source", interaction_metadata["source"].reindex(results["interaction_tuple"]).to_numpy())

    results = results.sort_values(
        ["rra_fdr", "rra_pval", "median_rank_pct", "top5_freq"],
        ascending=[True, True, False, False],
    ).reset_index(drop=True)

    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        results.to_parquet(output_path)

    return results