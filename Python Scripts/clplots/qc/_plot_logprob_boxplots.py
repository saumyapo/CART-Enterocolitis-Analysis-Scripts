"""QC boxplot helpers for raw interaction log-probabilities."""

from __future__ import annotations

from typing import Iterable, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def _validate_required_columns(data: pd.DataFrame, required: Iterable[str]) -> None:
    missing = [column for column in required if column not in data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def _prepare_plot_df(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
) -> pd.DataFrame:
    _validate_required_columns(df, [x_col, y_col])

    plot_df = df[[x_col, y_col]].dropna().copy()
    plot_df[x_col] = pd.to_numeric(plot_df[x_col], errors="coerce")
    plot_df = plot_df.dropna(subset=[x_col]).copy()

    if plot_df.empty:
        raise ValueError("No valid rows remain after dropping missing values.")

    plot_df[y_col] = plot_df[y_col].astype(str)
    return plot_df


def _resolve_order(
    plot_df: pd.DataFrame,
    y_col: str,
    order: Optional[Sequence[str]],
    order_by: str,
) -> Sequence[str]:
    if order is not None:
        return list(order)

    if order_by not in {"median", "mean", "count", "alphabetical"}:
        raise ValueError("order_by must be one of: 'median', 'mean', 'count', 'alphabetical'")

    if order_by == "alphabetical":
        return sorted(plot_df[y_col].unique().tolist())

    if order_by == "count":
        return plot_df[y_col].value_counts().sort_values(ascending=False).index.tolist()

    metric = "median" if order_by == "median" else "mean"
    return (
        plot_df.groupby(y_col, observed=True)["logprob"]
        .agg(metric)
        .sort_values(ascending=False)
        .index.tolist()
    )


def _default_figsize(n_groups: int, width: float, min_height: float, step_height: float) -> Tuple[float, float]:
    return (width, max(min_height, step_height * n_groups))


def _plot_logprob_boxplot(
    df: pd.DataFrame,
    y_col: str,
    x_col: str = "logprob",
    order: Optional[Sequence[str]] = None,
    order_by: str = "median",
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[float] = None,
    title: Optional[str] = None,
    xlabel: str = "log(CellChat Interaction Probability)",
    ylabel: Optional[str] = None,
    color: str = "#D9D9D9",
    box_linewidth: float = 0.9,
    median_linewidth: float = 1.1,
    whisker_linewidth: float = 0.9,
    cap_linewidth: float = 0.9,
    fliersize: float = 0.0,
    show_points: bool = False,
    point_size: float = 2.0,
    point_alpha: float = 0.35,
    point_color: str = "#4A4A4A",
    grid_alpha: float = 0.14,
    despine_offset: float = 5,
    ax=None,
):
    plot_df = _prepare_plot_df(df=df, x_col=x_col, y_col=y_col)
    plot_df = plot_df.rename(columns={x_col: "logprob", y_col: "group"})

    resolved_order = _resolve_order(
        plot_df=plot_df,
        y_col="group",
        order=order,
        order_by=order_by,
    )

    plot_df["group"] = pd.Categorical(
        plot_df["group"],
        categories=resolved_order,
        ordered=True,
    )
    plot_df = plot_df.sort_values("group").reset_index(drop=True)

    if figsize is None:
        figsize = _default_figsize(
            n_groups=len(resolved_order),
            width=4.6,
            min_height=1.6,
            step_height=0.34,
        )

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, constrained_layout=True)
    else:
        fig = ax.figure

    sns.boxplot(
        data=plot_df,
        x="logprob",
        y="group",
        order=resolved_order,
        color=color,
        width=0.62,
        linewidth=box_linewidth,
        fliersize=fliersize,
        boxprops={"edgecolor": "#333333", "facecolor": color},
        medianprops={"color": "#111111", "linewidth": median_linewidth},
        whiskerprops={"color": "#333333", "linewidth": whisker_linewidth},
        capprops={"color": "#333333", "linewidth": cap_linewidth},
        ax=ax,
    )

    if show_points:
        sns.stripplot(
            data=plot_df,
            x="logprob",
            y="group",
            order=resolved_order,
            color=point_color,
            size=point_size,
            alpha=point_alpha,
            jitter=0.18,
            linewidth=0,
            ax=ax,
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel if ylabel is not None else y_col.replace("_", " ").title())
    if title is not None:
        ax.set_title(title, fontsize=12, fontweight="bold", pad=4)

    ax.grid(axis="x", alpha=grid_alpha, linewidth=0.8)
    ax.grid(axis="y", visible=False)
    ax.tick_params(axis="both", labelsize=9)
    sns.despine(ax=ax, trim=True, offset=despine_offset)

    return fig, ax, plot_df.rename(columns={"logprob": x_col, "group": y_col})


def plot_logprob_condition_qc(
    df: pd.DataFrame,
    logprob_col: str = "logprob",
    condition_col: str = "condition",
    order: Optional[Sequence[str]] = None,
    order_by: str = "median",
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[float] = None,
    title: str = "Logprob by condition",
    xlabel: str = "log(CellChat Interaction Probability)",
    ylabel: str = "Condition",
    color: str = "#D9D9D9",
    box_linewidth: float = 0.9,
    median_linewidth: float = 1.1,
    whisker_linewidth: float = 0.9,
    cap_linewidth: float = 0.9,
    fliersize: float = 0.0,
    show_points: bool = False,
    point_size: float = 2.0,
    point_alpha: float = 0.35,
    point_color: str = "#4A4A4A",
    grid_alpha: float = 0.14,
    despine_offset: float = 5,
    ax=None,
):
    """Plot a publication-style QC boxplot of logprob by condition."""

    return _plot_logprob_boxplot(
        df=df,
        y_col=condition_col,
        x_col=logprob_col,
        order=order,
        order_by=order_by,
        figsize=figsize,
        dpi=dpi,
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        color=color,
        box_linewidth=box_linewidth,
        median_linewidth=median_linewidth,
        whisker_linewidth=whisker_linewidth,
        cap_linewidth=cap_linewidth,
        fliersize=fliersize,
        show_points=show_points,
        point_size=point_size,
        point_alpha=point_alpha,
        point_color=point_color,
        grid_alpha=grid_alpha,
        despine_offset=despine_offset,
        ax=ax,
    )


def plot_logprob_sample_qc(
    df: pd.DataFrame,
    logprob_col: str = "logprob",
    sample_col: str = "sample",
    order: Optional[Sequence[str]] = None,
    order_by: str = "median",
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[float] = None,
    title: str = "Logprob by sample",
    xlabel: str = "log(CellChat Interaction Probability)",
    ylabel: str = "Sample",
    color: str = "#D9D9D9",
    box_linewidth: float = 0.9,
    median_linewidth: float = 1.1,
    whisker_linewidth: float = 0.9,
    cap_linewidth: float = 0.9,
    fliersize: float = 0.0,
    show_points: bool = False,
    point_size: float = 2.0,
    point_alpha: float = 0.3,
    point_color: str = "#4A4A4A",
    grid_alpha: float = 0.14,
    despine_offset: float = 5,
    ax=None,
):
    """Plot a publication-style QC boxplot of logprob by sample."""

    return _plot_logprob_boxplot(
        df=df,
        y_col=sample_col,
        x_col=logprob_col,
        order=order,
        order_by=order_by,
        figsize=figsize,
        dpi=dpi,
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        color=color,
        box_linewidth=box_linewidth,
        median_linewidth=median_linewidth,
        whisker_linewidth=whisker_linewidth,
        cap_linewidth=cap_linewidth,
        fliersize=fliersize,
        show_points=show_points,
        point_size=point_size,
        point_alpha=point_alpha,
        point_color=point_color,
        grid_alpha=grid_alpha,
        despine_offset=despine_offset,
        ax=ax,
    )
