from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Union

import pandas as pd
from tqdm.auto import tqdm


PathLike = Union[str, Path]


def run_median_rank_analysis(
    source_file: PathLike,
    output_file: Optional[PathLike] = None,
    condition: Optional[Union[str, Sequence[str]]] = "NV",
    min_samples: int = 3,
    show_progress: bool = True,
) -> pd.DataFrame:
    """Summarize interaction ranks across samples using median-ranks.

    Args:
        source_file: Path to the merged CellChat parquet file.
        output_file: Optional destination parquet path for the summary table.
        condition: Optional condition label or labels used to subset the merged
            table before ranking. When ``None``, all rows are included.
        min_samples: Minimum number of observations required for an interaction
            to be retained in the summary table.
        show_progress: Whether to show a tqdm progress bar while aggregating
            interaction groups.

    Returns:
        Dataframe of interaction-level median-rank summaries.

    Raises:
        FileNotFoundError: If the source parquet does not exist.
        ValueError: If required columns are missing, filtering removes all rows,
            or the minimum sample threshold is invalid.
    """
    source_path = Path(source_file)
    output_path = Path(output_file) if output_file is not None else None

    if not source_path.exists():
        raise FileNotFoundError(f"Source parquet does not exist: {source_path}")

    if min_samples < 1:
        raise ValueError("min_samples must be at least 1.")

    df = pd.read_parquet(source_path)

    required_columns = {
        "sample",
        "source",
        "target",
        "interaction_name",
        "ligand",
        "receptor",
        "pathway_name",
        "prob",
        "logprob",
    }
    if condition is not None:
        required_columns.add("condition")

    missing_columns = sorted(required_columns.difference(df.columns))
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    if condition is None:
        df = df.copy()
    elif isinstance(condition, str):
        df = df[df["condition"].eq(condition)].copy()
    else:
        conditions = tuple(condition)
        df = df[df["condition"].isin(conditions)].copy()

    if df.empty:
        raise ValueError("No rows remain after filtering the input dataframe.")

    df["rank"] = df.groupby("sample")["logprob"].rank(pct=True)

    group_cols = ["source", "target", "interaction_name"]
    grouped = df.groupby(group_cols, sort=False)
    iterator = grouped
    if show_progress:
        iterator = tqdm(grouped, total=grouped.ngroups, desc="Median-rank summaries")

    results_raw = []
    for (source, target, interaction_name), subdf in iterator:
        if subdf.shape[0] < min_samples:
            continue

        results_raw.append(
            {
                "source": source,
                "target": target,
                "interaction_name": interaction_name,
                "ligand": subdf["ligand"].iloc[0],
                "receptor": subdf["receptor"].iloc[0],
                "pathway_name": subdf["pathway_name"].iloc[0],
                "rank_median": subdf["rank"].median(),
                "rank_mean": subdf["rank"].mean(),
                "rank_std": subdf["rank"].std(),
                "prob": subdf["prob"].mean(),
                "logprob": subdf["logprob"].mean(),
                "n_samples": int(subdf["sample"].nunique()),
            }
        )

    results = pd.DataFrame(results_raw)

    if not results.empty:
        results = results.sort_values(
            ["rank_mean", "rank_median", "logprob"],
            ascending=[False, False, False],
        ).reset_index(drop=True)

    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        results.to_parquet(output_path)

    return results