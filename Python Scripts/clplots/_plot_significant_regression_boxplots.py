"""Boxplot helper for significant regression effects by celltype."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ResultsLike = Union[pd.DataFrame, str, Path]


def _load_results(results: ResultsLike) -> pd.DataFrame:
    if isinstance(results, pd.DataFrame):
        return results.copy()

    result_path = Path(results)
    if not result_path.exists():
        raise FileNotFoundError(f"Results file does not exist: {result_path}")

    if result_path.suffix == ".parquet":
        return pd.read_parquet(result_path)
    if result_path.suffix == ".csv":
        return pd.read_csv(result_path)

    raise ValueError("results must be a DataFrame or a .parquet/.csv file path")


def _validate_required_columns(data: pd.DataFrame, required: Iterable[str]) -> None:
    missing = [col for col in required if col not in data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def _normalize_order_by(order_by: str, treated_label: str, colitis_label: str) -> str:
    aliases = {
        "count": "count",
        "median_abs": "median_abs",
        "both": "both",
        "combined": "both",
        "median": "both",
        "median_signed": "both",
        "treated": "treated",
        "treated_median": "treated",
        treated_label.casefold(): "treated",
        "colitis": "colitis",
        "colitis_median": "colitis",
        colitis_label.casefold(): "colitis",
    }

    normalized = aliases.get(order_by.strip().casefold())
    if normalized is None:
        raise ValueError(
            "order_by must be one of: 'median_abs', 'count', 'treated', 'colitis', 'both'"
        )

    return normalized


def _ordered_celltypes(primary: pd.Series, fallback: Optional[pd.Series] = None) -> list[str]:
    ordered = primary.dropna().sort_values(ascending=False).index.tolist()

    if fallback is not None:
        for celltype in fallback.dropna().sort_values(ascending=False).index.tolist():
            if celltype not in ordered:
                ordered.append(celltype)

    return ordered


def plot_significant_regression_boxplots(
    results: ResultsLike,
    celltype_col: str = "target",
    q_thresh: float = 0.05,
    treated_coef_col: str = "coef_is_treated",
    treated_qval_col: str = "qval_is_treated",
    colitis_coef_col: str = "coef_has_colitis",
    colitis_qval_col: str = "qval_has_colitis",
    treated_label: str = "Treated",
    colitis_label: str = "Colitis",
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[float] = None,
    title: Optional[str] = None,
    xlabel: str = "Regression coefficient",
    ylabel: Optional[str] = None,
    palette: Optional[dict] = None,
    order: Optional[Sequence[str]] = None,
    order_by: str = "median_abs",
    show_points: bool = False,
    point_size: float = 2.5,
    point_alpha: float = 0.35,
    box_linewidth: float = 0.85,
    fliersize: float = 0,
    legend: bool = True,
    despine_trim: bool = True,
    despine_offset: float = 5,
    ax=None,
):
    """Plot significant treated and colitis regression coefficients by celltype.

    The regression results are reshaped from wide to long format so both effects
    share a single coefficient column and an effect label suitable for hue.

    order_by controls the default y-axis order when order is not supplied.
    Supported values are:
    - 'median_abs': median absolute coefficient across both effects.
    - 'count': number of significant rows per celltype.
    - 'treated': median treated coefficient, highest median at the top.
    - 'colitis': median colitis coefficient, highest median at the top.
    - 'both': median signed coefficient across both effects, highest median at the top.

    Returns
    -------
    fig, ax, plot_df
        Figure, axis, and the reshaped significant rows used for plotting.
    """

    data = _load_results(results)
    _validate_required_columns(
        data,
        [
            celltype_col,
            treated_coef_col,
            treated_qval_col,
            colitis_coef_col,
            colitis_qval_col,
        ],
    )

    order_by = _normalize_order_by(order_by, treated_label=treated_label, colitis_label=colitis_label)

    if q_thresh <= 0:
        raise ValueError("q_thresh must be > 0")

    if palette is None:
        palette = {
            treated_label: "#2B6CB0",
            colitis_label: "#C93C37",
        }

    effect_specs = [
        (treated_label, treated_coef_col, treated_qval_col),
        (colitis_label, colitis_coef_col, colitis_qval_col),
    ]

    optional_cols = [
        col for col in ["source", "target", "interaction_name", "pathway_name"] if col in data.columns
    ]

    plot_frames = []
    for effect_label, coef_col, qval_col in effect_specs:
        sig = data.loc[data[qval_col].lt(q_thresh)].copy()
        if sig.empty:
            continue

        sig["effect"] = effect_label
        sig["coefficient"] = pd.to_numeric(sig[coef_col], errors="coerce")
        sig["qvalue"] = pd.to_numeric(sig[qval_col], errors="coerce")
        sig["celltype"] = sig[celltype_col]
        plot_frames.append(sig[["effect", "coefficient", "qvalue", "celltype", *optional_cols]])

    if not plot_frames:
        raise ValueError(f"No significant rows found at q < {q_thresh}")

    plot_df = pd.concat(plot_frames, ignore_index=True)
    plot_df = plot_df.dropna(subset=["coefficient", "qvalue", "celltype"]).copy()

    if plot_df.empty:
        raise ValueError("No valid rows remain after reshaping and dropping missing values.")

    if order is None:
        if order_by == "count":
            order = (
                plot_df["celltype"]
                .value_counts()
                .sort_values(ascending=False)
                .index
                .tolist()
            )
        elif order_by == "median_abs":
            order = (
                plot_df.groupby("celltype", observed=True)["coefficient"]
                .apply(lambda values: values.abs().median())
                .sort_values(ascending=False)
                .index
                .tolist()
            )
        elif order_by == "both":
            order = (
                plot_df.groupby("celltype", observed=True)["coefficient"]
                .median()
                .sort_values(ascending=False)
                .index
                .tolist()
            )
        else:
            overall_median = (
                plot_df.groupby("celltype", observed=True)["coefficient"]
                .median()
            )

            effect_label = treated_label if order_by == "treated" else colitis_label
            effect_median = (
                plot_df.loc[plot_df["effect"].eq(effect_label)]
                .groupby("celltype", observed=True)["coefficient"]
                .median()
            )
            order = _ordered_celltypes(effect_median, fallback=overall_median)

    if figsize is None:
        figsize = (9, max(4, 0.35 * len(order)))

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, constrained_layout=True)
    else:
        fig = ax.figure

    sns.boxplot(
        data=plot_df,
        x="coefficient",
        y="celltype",
        hue="effect",
        order=order,
        hue_order=[treated_label, colitis_label],
        palette=palette,
        linewidth=box_linewidth,
        fliersize=fliersize,
        ax=ax,
    )

    if show_points:
        sns.stripplot(
            data=plot_df,
            x="coefficient",
            y="celltype",
            hue="effect",
            order=order,
            hue_order=[treated_label, colitis_label],
            dodge=True,
            palette=palette,
            size=point_size,
            alpha=point_alpha,
            linewidth=0,
            ax=ax,
        )

    ax.axvline(0, color="0.45", linestyle="--", linewidth=1.0, zorder=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel if ylabel is not None else celltype_col.capitalize())
    ax.set_title(
        title if title is not None else f"Significant interaction effects by {celltype_col}"
    )
    ax.grid(axis="x", alpha=0.14, linewidth=0.8)
    ax.grid(axis="y", visible=False)
    sns.despine(ax=ax, trim=despine_trim, offset=despine_offset)

    if legend:
        handles, labels = ax.get_legend_handles_labels()
        deduped = []
        seen = set()
        for handle, label in zip(handles, labels):
            if label in seen or label not in {treated_label, colitis_label}:
                continue
            deduped.append((handle, label))
            seen.add(label)

        if deduped:
            legend_handles, legend_labels = zip(*deduped)
            ax.legend(legend_handles, legend_labels, frameon=False, title=None)
    else:
        legend_obj = ax.get_legend()
        if legend_obj is not None:
            legend_obj.remove()

    return fig, ax, plot_df