from __future__ import annotations

from collections.abc import Mapping
import importlib.util
from pathlib import Path
import re
from typing import Union

import numpy as np
import pandas as pd


PathLike = Union[str, Path]


def _count_positive(values: pd.Series) -> int:
    """Count positive coefficients within a grouped series."""
    return int((values > 0).sum())


def _count_negative(values: pd.Series) -> int:
    """Count negative coefficients within a grouped series."""
    return int((values < 0).sum())


def _percent(numerator, denominator):
    """Compute percentages while avoiding divide-by-zero artifacts."""
    if np.isscalar(denominator):
        return 100 * numerator / denominator if denominator else np.nan

    safe_denominator = denominator.replace(0, np.nan)
    return 100 * numerator / safe_denominator


def _require_columns(df: pd.DataFrame, columns: list[str]) -> None:
    """Raise a clear error when required regression-summary columns are missing."""
    missing_columns = [column for column in columns if column not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")


def _build_direction_column(coef: pd.Series, positive_label: str) -> pd.Series:
    """Label coefficients by direction of change in the positive condition.

    Zero-valued coefficients are kept in the non-positive bucket to preserve the
    historical behavior of this summarizer.
    """
    return pd.Series(
        np.where(
            coef > 0,
            f"Increased in {positive_label}",
            f"Decreased in {positive_label}",
        ),
        index=coef.index,
        name="direction",
    )


def _build_overall_summary(
    df: pd.DataFrame,
    coef_col: str,
    positive_label: str,
) -> pd.DataFrame:
    """Summarize overall counts and effect-size extremes across all rows."""
    total = len(df)
    n_increased = _count_positive(df[coef_col])
    n_decreased = _count_negative(df[coef_col])

    return pd.DataFrame(
        {
            "metric": [
                "Total significant interactions",
                f"Increased in {positive_label}",
                f"Decreased in {positive_label}",
                "Median effect size",
                "Mean effect size",
                "Max positive effect",
                "Max negative effect",
            ],
            "value": [
                total,
                n_increased,
                n_decreased,
                df[coef_col].median(),
                df[coef_col].mean(),
                df[coef_col].max(),
                df[coef_col].min(),
            ],
            "percent": [
                100,
                _percent(n_increased, total),
                _percent(n_decreased, total),
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ],
        }
    )


def _build_group_summary(
    df: pd.DataFrame,
    group_cols,
    coef_col: str,
    q_col: str,
    total_significant: int,
    include_mean_effect: bool = True,
) -> pd.DataFrame:
    """Aggregate counts, effect sizes, and minimum q-values by group."""
    aggregation = {
        "n_interactions": (coef_col, "size"),
        "n_increased": (coef_col, _count_positive),
        "n_decreased": (coef_col, _count_negative),
        "median_effect": (coef_col, "median"),
    }
    if include_mean_effect:
        aggregation["mean_effect"] = (coef_col, "mean")
    aggregation["min_q"] = (q_col, "min")

    summary = df.groupby(group_cols).agg(**aggregation).reset_index()
    summary["percent_of_all_significant"] = _percent(
        summary["n_interactions"],
        total_significant,
    )
    summary["percent_increased"] = _percent(
        summary["n_increased"],
        summary["n_interactions"],
    )
    summary["percent_decreased"] = _percent(
        summary["n_decreased"],
        summary["n_interactions"],
    )
    return summary.sort_values("n_interactions", ascending=False)


def _build_directional_summary(
    df: pd.DataFrame,
    group_col: str,
    top_n: int | None = None,
) -> pd.DataFrame:
    """Count increased and decreased interactions within each category."""
    summary = (
        df.groupby([group_col, "direction"])
        .size()
        .unstack(fill_value=0)
        .assign(total=lambda table: table.sum(axis=1))
        .sort_values("total", ascending=False)
    )
    return summary.head(top_n) if top_n is not None else summary


def _select_top_interactions(
    df: pd.DataFrame,
    *,
    sort_col: str,
    ascending: bool,
    columns: list[str],
    top_n: int,
    mask: pd.Series | None = None,
) -> pd.DataFrame:
    """Return the most relevant interaction rows after sorting and filtering."""
    subset = df.loc[mask] if mask is not None else df
    return subset.sort_values(sort_col, ascending=ascending).loc[:, columns].head(top_n)


def _filter_ta_targets(df: pd.DataFrame, target_col: str) -> pd.DataFrame:
    """Keep only interactions whose target label contains TA cells."""
    ta_mask = df[target_col].astype(str).str.contains(r"\bTA cells\b", case=False, na=False)
    return df.loc[ta_mask].copy()


def _make_excel_sheet_name(name: str, used_names: set[str]) -> str:
    """Normalize sheet names so each summary key maps to a valid Excel tab."""
    cleaned_name = re.sub(r"[\\/*?:\[\]]", "_", str(name)).strip().strip("'")
    cleaned_name = cleaned_name[:31] or "Sheet"

    candidate = cleaned_name
    suffix = 1
    while candidate in used_names:
        suffix_text = f"_{suffix}"
        candidate = f"{cleaned_name[: 31 - len(suffix_text)]}{suffix_text}"
        suffix += 1

    used_names.add(candidate)
    return candidate


def _resolve_excel_writer_engine(engine: str | None) -> str:
    """Choose an installed pandas Excel writer engine for xlsx output."""
    if engine is not None:
        return engine

    for candidate in ("xlsxwriter", "openpyxl"):
        if importlib.util.find_spec(candidate) is not None:
            return candidate

    raise ImportError(
        "Writing .xlsx files requires either 'xlsxwriter' or 'openpyxl'. "
        "Install one of those packages or pass a supported engine explicitly."
    )


def export_regression_analysis_summaries_to_excel(
    summaries: Mapping[str, pd.DataFrame],
    output_file: PathLike,
    include_index: bool = True,
    engine: str | None = None,
) -> Path:
    """Write regression summary tables to a multi-sheet Excel workbook.

    Each key in ``summaries`` becomes its own worksheet. Sheet names are based
    on the dictionary keys and are normalized only when needed to satisfy Excel
    naming rules or avoid duplicate tabs.

    Args:
        summaries: Mapping of summary names to dataframes, such as the output of
            ``summarize_regression_analysis``.
        output_file: Destination path for the ``.xlsx`` workbook.
        include_index: Whether to write dataframe indexes into each worksheet.
        engine: Optional pandas Excel writer engine, for example
            ``"xlsxwriter"`` or ``"openpyxl"``.

    Returns:
        Path to the workbook that was written.

    Raises:
        ValueError: If ``summaries`` is empty.
        TypeError: If any summary value is not a pandas dataframe.
    """
    if not summaries:
        raise ValueError("summaries must contain at least one dataframe.")

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    writer_kwargs = {"engine": _resolve_excel_writer_engine(engine)}
    used_sheet_names: set[str] = set()

    with pd.ExcelWriter(output_path, **writer_kwargs) as writer:
        for summary_name, summary_df in summaries.items():
            if not isinstance(summary_df, pd.DataFrame):
                raise TypeError(
                    "Each summary must be a pandas DataFrame. "
                    f"Key '{summary_name}' has type {type(summary_df).__name__}."
                )

            sheet_name = _make_excel_sheet_name(summary_name, used_sheet_names)
            summary_df.to_excel(writer, sheet_name=sheet_name, index=include_index)

    return output_path


def summarize_regression_analysis(
    df,
    coef_col="coef_is_treated",
    p_col="pval_is_treated",
    q_col="qval_is_treated",
    source_col="source",
    target_col="target",
    pathway_col="pathway_name",
    annotation_col="annotation",
    interaction_col="interaction_name",
    positive_label="CAR-TCTRL",
    reference_label="NV",
    top_n=10,
):
    """Build a compact set of summaries for regression results.

    This helper expects a dataframe that has already been filtered to the
    regression hits that you consider significant. Positive coefficients are
    interpreted as stronger interactions in ``positive_label`` relative to
    ``reference_label``. Negative coefficients are interpreted as weaker
    interactions in ``positive_label`` relative to ``reference_label``.

    Args:
        df: Significant interaction-level regression results.
        coef_col: Column containing the treatment coefficient or effect size.
        p_col: Column containing nominal p-values.
        q_col: Column containing multiple-testing adjusted q-values.
        source_col: Column identifying the sender cell type.
        target_col: Column identifying the receiver cell type.
        pathway_col: Column identifying the signaling pathway.
        annotation_col: Column describing the interaction annotation.
        interaction_col: Column identifying the ligand-receptor interaction.
        positive_label: Name of the condition associated with positive effects.
        reference_label: Name of the baseline condition used to interpret the
            sign of the coefficient. It is kept in the public API for context
            and compatibility, even though it does not change output column
            names.
        top_n: Maximum number of rows to retain for rank-ordered summary tables.

    Returns:
        Dictionary of summary tables with stable keys:
            ``overall``: Global counts and effect-size statistics.
            ``annotation_summary``: Top annotations by interaction count.
            ``pathway_summary``: Top pathways by interaction count.
            ``source_summary``: Top sender cell types by interaction count.
            ``target_summary``: Top receiver cell types by interaction count.
            ``source_target_summary``: Top sender-receiver pairs.
            ``directional_by_pathway``: Increased vs. decreased counts by pathway.
            ``directional_by_annotation``: Increased vs. decreased counts by annotation.
            ``top_by_q``: Most statistically significant interactions.
            ``top_positive_effects``: Largest positive coefficients.
            ``top_negative_effects``: Most negative coefficients.
            ``ta_summary``: Overall summary for TA-cell-directed interactions.
            ``ta_top_sources``: Top senders targeting TA cells.
            ``ta_top_interactions``: Most significant TA-cell-directed interactions.

    Raises:
        ValueError: If any required summary columns are missing.
    """
    required_columns = [
        coef_col,
        p_col,
        q_col,
        source_col,
        target_col,
        pathway_col,
        annotation_col,
        interaction_col,
    ]
    _require_columns(df, required_columns)

    summary_df = df.copy()
    summary_df["direction"] = _build_direction_column(summary_df[coef_col], positive_label)
    total_significant = len(summary_df)

    interaction_columns = [
        source_col,
        target_col,
        interaction_col,
        pathway_col,
        annotation_col,
        coef_col,
        p_col,
        q_col,
    ]
    interaction_columns_with_direction = [*interaction_columns, "direction"]

    source_summary = _build_group_summary(
        summary_df,
        source_col,
        coef_col,
        q_col,
        total_significant,
    ).head(top_n)
    target_summary = _build_group_summary(
        summary_df,
        target_col,
        coef_col,
        q_col,
        total_significant,
    ).head(top_n)
    pathway_summary = _build_group_summary(
        summary_df,
        pathway_col,
        coef_col,
        q_col,
        total_significant,
    ).head(top_n)
    annotation_summary = _build_group_summary(
        summary_df,
        annotation_col,
        coef_col,
        q_col,
        total_significant,
    ).head(top_n)

    source_target_summary = _build_group_summary(
        summary_df,
        [source_col, target_col],
        coef_col,
        q_col,
        total_significant,
    ).head(top_n)

    directional_by_pathway = _build_directional_summary(
        summary_df,
        pathway_col,
        top_n=top_n,
    )
    directional_by_annotation = _build_directional_summary(summary_df, annotation_col)

    top_by_q = _select_top_interactions(
        summary_df,
        sort_col=q_col,
        ascending=True,
        columns=interaction_columns_with_direction,
        top_n=top_n,
    )
    top_positive_effects = _select_top_interactions(
        summary_df,
        sort_col=coef_col,
        ascending=False,
        columns=interaction_columns,
        top_n=top_n,
        mask=summary_df[coef_col] > 0,
    )
    top_negative_effects = _select_top_interactions(
        summary_df,
        sort_col=coef_col,
        ascending=True,
        columns=interaction_columns,
        top_n=top_n,
        mask=summary_df[coef_col] < 0,
    )


    return {
        "overall": _build_overall_summary(summary_df, coef_col, positive_label),
        "annotation_summary": annotation_summary,
        "pathway_summary": pathway_summary,
        "source_summary": source_summary,
        "target_summary": target_summary,
        "source_target_summary": source_target_summary,
        "directional_by_pathway": directional_by_pathway,
        "directional_by_annotation": directional_by_annotation,
        "top_by_q": top_by_q,
        "top_positive_effects": top_positive_effects,
        "top_negative_effects": top_negative_effects,
    }