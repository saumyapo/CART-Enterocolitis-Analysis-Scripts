
from __future__ import annotations

import math
import inspect
from typing import Iterable, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
import pandas as pd
import seaborn as sns


def _validate_required_columns(data: pd.DataFrame, required: Iterable[str]) -> None:
    missing = [column for column in required if column not in data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def _prepare_barplot_df(df: pd.DataFrame, x_col: str, y_col: str) -> pd.DataFrame:
    _validate_required_columns(df, [x_col, y_col])

    plot_df = df[[x_col, y_col]].dropna().copy()
    plot_df[x_col] = pd.to_numeric(plot_df[x_col], errors="coerce")
    plot_df = plot_df.dropna(subset=[x_col]).copy()

    if plot_df.empty:
        raise ValueError("No valid rows remain after dropping missing values.")

    plot_df[y_col] = plot_df[y_col].astype(str)
    return plot_df


def _resolve_symmetric_limit(values: Sequence[float]) -> float:
    max_abs = max(abs(float(value)) for value in values)
    if max_abs == 0:
        return 1.0
    if max_abs >= 1.0:
        return float(math.ceil(max_abs))
    return math.ceil(max_abs * 2.0) / 2.0


def _union_order(left_order: Sequence[str], right_order: Sequence[str]) -> list[str]:
    union = []
    seen = set()
    for label in list(left_order) + list(right_order):
        if label in seen:
            continue
        seen.add(label)
        union.append(label)
    return union


def plot_publication_barplot(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    figsize: Tuple[float, float] = (3.0, 6.0),
    dpi: Optional[float] = None,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    axis_labelsize: float = 10,
    tick_labelsize: float = 9,
    x_tick_labelsize: Optional[float] = None,
    y_tick_labelsize: Optional[float] = None,
    title_fontsize: float = 12,
    color: str = "#8E9AAF",
    edgecolor: str = "#2F2F2F",
    linewidth: float = 0.8,
    saturation: float = 1.0,
    errorbar=None,
    order: Optional[Sequence[str]] = None,
    grid_axis: str = "x",
    grid_alpha: float = 0.14,
    despine_offset: float = 5,
    ax=None,
):
    """Plot a clean publication-style barplot from a summarized dataframe."""

    if x_tick_labelsize is None:
        x_tick_labelsize = tick_labelsize
    if y_tick_labelsize is None:
        y_tick_labelsize = tick_labelsize

    plot_df = _prepare_barplot_df(df=df, x_col=x_col, y_col=y_col)

    if order is None:
        resolved_order = plot_df[y_col].drop_duplicates().tolist()
    else:
        resolved_order = list(order)

    plot_df[y_col] = pd.Categorical(plot_df[y_col], categories=resolved_order, ordered=True)
    plot_df = plot_df.sort_values(y_col).reset_index(drop=True)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, constrained_layout=True)
    else:
        fig = ax.figure

    barplot_kwargs = {
        "data": plot_df,
        "x": x_col,
        "y": y_col,
        "order": resolved_order,
        "color": color,
        "saturation": saturation,
        "edgecolor": edgecolor,
        "linewidth": linewidth,
        "ax": ax,
    }

    if "errorbar" in inspect.signature(sns.barplot).parameters:
        barplot_kwargs["errorbar"] = errorbar
    else:
        barplot_kwargs["ci"] = errorbar

    sns.barplot(**barplot_kwargs)

    ax.set_xlabel(xlabel if xlabel is not None else x_col.replace("_", " ").title(), fontsize=axis_labelsize)
    ax.set_ylabel(ylabel if ylabel is not None else y_col.replace("_", " ").title(), fontsize=axis_labelsize)
    ax.tick_params(axis="x", labelsize=x_tick_labelsize)
    ax.tick_params(axis="y", labelsize=y_tick_labelsize)

    if title is not None:
        ax.set_title(title, fontsize=title_fontsize, fontweight="bold", pad=4)

    ax.grid(axis=grid_axis, alpha=grid_alpha, linewidth=0.8)
    if grid_axis == "x":
        ax.grid(axis="y", visible=False)
    elif grid_axis == "y":
        ax.grid(axis="x", visible=False)

    sns.despine(ax=ax, trim=True, offset=despine_offset)

    return fig, ax, plot_df


def plot_publication_barplot_pair(
    left_df: pd.DataFrame,
    right_df: pd.DataFrame,
    x_col: str,
    left_y_col: str,
    right_y_col: str,
    figsize: Tuple[float, float] = (3.8, 6.0),
    dpi: Optional[float] = None,
    title: Optional[str] = None,
    left_xlabel: str = "Median Effect",
    right_xlabel: str = "Median Effect",
    left_ylabel: Optional[str] = None,
    axis_labelsize: float = 10,
    tick_labelsize: float = 9,
    x_tick_labelsize: Optional[float] = None,
    y_tick_labelsize: Optional[float] = None,
    title_fontsize: float = 12,
    left_color: str = "#8E9AAF",
    right_color: str = "#C3CBD6",
    edgecolor: str = "#2F2F2F",
    linewidth: float = 0.8,
    left_order: Optional[Sequence[str]] = None,
    right_order: Optional[Sequence[str]] = None,
    left_panel_title: Optional[str] = None,
    right_panel_title: Optional[str] = None,
    panel_title_fontsize: float = 10,
    grid_alpha: float = 0.14,
    despine_offset: float = 5,
    wspace: float = 0.08,
    x_tick_pad: float = 8,
    match_x_limits: bool = True,
):
    """Plot paired publication-style horizontal barplots with shared y positions."""

    if x_tick_labelsize is None:
        x_tick_labelsize = tick_labelsize
    if y_tick_labelsize is None:
        y_tick_labelsize = tick_labelsize

    left_plot_df = _prepare_barplot_df(df=left_df, x_col=x_col, y_col=left_y_col)
    right_plot_df = _prepare_barplot_df(df=right_df, x_col=x_col, y_col=right_y_col)

    if left_order is None:
        left_order = left_plot_df[left_y_col].drop_duplicates().tolist()
    else:
        left_order = list(left_order)

    if right_order is None:
        right_order = right_plot_df[right_y_col].drop_duplicates().tolist()
    else:
        right_order = list(right_order)

    union_order = _union_order(left_order=left_order, right_order=right_order)

    left_values = left_plot_df.drop_duplicates(subset=[left_y_col]).set_index(left_y_col)[x_col].reindex(union_order)
    right_values = right_plot_df.drop_duplicates(subset=[right_y_col]).set_index(right_y_col)[x_col].reindex(union_order)

    left_plot_df = pd.DataFrame({left_y_col: union_order, x_col: left_values.values})
    right_plot_df = pd.DataFrame({right_y_col: union_order, x_col: right_values.values})

    n_rows = len(union_order)
    positions = list(range(n_rows))
    left_mask = left_plot_df[x_col].notna()
    right_mask = right_plot_df[x_col].notna()

    fig, axes = plt.subplots(
        1,
        2,
        figsize=figsize,
        dpi=dpi,
        gridspec_kw={"wspace": wspace},
        constrained_layout=False,
    )
    left_ax, right_ax = axes

    left_ax.barh(
        left_plot_df.index[left_mask],
        left_plot_df.loc[left_mask, x_col],
        color=left_color,
        edgecolor=edgecolor,
        linewidth=linewidth,
        height=0.72,
    )
    right_ax.barh(
        right_plot_df.index[right_mask],
        right_plot_df.loc[right_mask, x_col],
        color=right_color,
        edgecolor=edgecolor,
        linewidth=linewidth,
        height=0.72,
    )

    for ax in axes:
        ax.set_ylim(-0.5, n_rows - 0.5)
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=grid_alpha, linewidth=0.8)
        ax.grid(axis="y", visible=False)
        ax.tick_params(axis="x", labelsize=x_tick_labelsize, pad=x_tick_pad)

    left_ax.set_yticks(positions)
    left_ax.set_yticklabels(union_order, fontsize=y_tick_labelsize)
    left_ax.tick_params(
        axis="y",
        length=2.5,
        width=0.8,
        pad=2,
        labelleft=True,
        labelright=False,
        labelsize=y_tick_labelsize,
    )
    for tick_label in left_ax.get_yticklabels():
        tick_label.set_fontsize(y_tick_labelsize)

    right_ax.set_yticks(positions)
    right_ax.tick_params(axis="y", left=False, right=False, length=0, pad=1, labelleft=False, labelright=False)
    right_ax.set_yticklabels([])

    left_ax.set_xlabel(left_xlabel, fontsize=axis_labelsize)
    right_ax.set_xlabel(right_xlabel, fontsize=axis_labelsize)
    left_ax.set_ylabel(left_ylabel if left_ylabel is not None else left_y_col.replace("_", " ").title(), fontsize=axis_labelsize)
    right_ax.set_ylabel("")

    if left_panel_title is not None:
        left_ax.set_title(left_panel_title, fontsize=panel_title_fontsize, pad=4)
    if right_panel_title is not None:
        right_ax.set_title(right_panel_title, fontsize=panel_title_fontsize, pad=4)
    if title is not None:
        fig.suptitle(title, fontsize=title_fontsize, fontweight="bold", y=0.98)

    if match_x_limits:
        symmetric_limit = _resolve_symmetric_limit(
            [left_plot_df[x_col].min(skipna=True), left_plot_df[x_col].max(skipna=True), right_plot_df[x_col].min(skipna=True), right_plot_df[x_col].max(skipna=True)]
        )
        x_ticks = [-symmetric_limit, 0.0, symmetric_limit]
        left_ax.set_xlim(-symmetric_limit, symmetric_limit)
        right_ax.set_xlim(-symmetric_limit, symmetric_limit)
        left_ax.set_xticks(x_ticks)
        right_ax.set_xticks(x_ticks)
        tick_labels = [f"{tick:g}" for tick in x_ticks]
        left_ax.set_xticklabels(tick_labels)
        right_ax.set_xticklabels(tick_labels)
        left_ticklabels = left_ax.get_xticklabels()
        right_ticklabels = right_ax.get_xticklabels()
        if left_ticklabels:
            left_ticklabels[-1].set_ha("right")
            left_ticklabels[-1].set_transform(
                left_ticklabels[-1].get_transform() + ScaledTranslation(-3 / 72, 0, fig.dpi_scale_trans)
            )
        if right_ticklabels:
            right_ticklabels[0].set_ha("left")
            right_ticklabels[0].set_transform(
                right_ticklabels[0].get_transform() + ScaledTranslation(3 / 72, 0, fig.dpi_scale_trans)
            )

    sns.despine(ax=left_ax, trim=True, offset=despine_offset, right=True)
    sns.despine(ax=right_ax, trim=True, offset=despine_offset, left=True)
    fig.subplots_adjust(wspace=wspace)

    return fig, axes, {"left": left_plot_df, "right": right_plot_df, "order": union_order}