from __future__ import annotations

import textwrap
from typing import Iterable, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def _validate_required_columns(data: pd.DataFrame, required: Iterable[str]) -> None:
    missing = [column for column in required if column not in data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def _merge_order(order: Optional[Sequence[str]], observed: Sequence[str]) -> list[str]:
    observed_values = [str(value) for value in observed]
    if order is None:
        return observed_values

    resolved = []
    seen = set()

    for value in order:
        text = str(value)
        if text not in seen:
            resolved.append(text)
            seen.add(text)

    for value in observed_values:
        if value not in seen:
            resolved.append(value)
            seen.add(value)

    return resolved


def _wrap_labels(values: Sequence[str], width: Optional[int]) -> list[str]:
    if width is None:
        return [str(value) for value in values]

    return ["\n".join(textwrap.wrap(str(value), width=width)) for value in values]


def _max_label_line_length(values: Sequence[str]) -> int:
    if not values:
        return 0

    return max(
        max(len(line) for line in str(value).splitlines() or [""])
        for value in values
    )


def _default_figsize(n_rows: int, n_cols: int, orientation: str) -> Tuple[float, float]:
    if orientation == "horizontal":
        return (
            max(6.0, 0.22 * max(n_cols, 1) + 1.8),
            max(2.4, 0.45 * max(n_rows, 1) + 1.2),
        )

    return (
        max(3.8, 0.6 * max(n_cols, 1) + 1.8),
        max(4.5, 0.22 * max(n_rows, 1) + 1.4),
    )


def plot_sample_count_heatmap_qc(
    df: pd.DataFrame,
    value_col: str = "prob",
    condition_col: str = "condition",
    sample_col: str = "sample",
    feature_col: str = "target",
    presence_threshold: float = 0.0,
    condition_order: Optional[Sequence[str]] = None,
    feature_order: Optional[Sequence[str]] = None,
    orientation: str = "horizontal",
    cluster_rows: Optional[bool] = None,
    cluster_cols: Optional[bool] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[float] = None,
    title: str = "Samples with at least 1 interaction",
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    cmap: str = "inferno",
    annot: bool = True,
    fmt: str = ".0f",
    annot_fontsize: float = 6,
    annot_color: Optional[str] = None,
    linewidths: float = 0.1,
    linecolor: str = "white",
    row_label_wrap: Optional[int] = None,
    col_label_wrap: Optional[int] = None,
    yticklabel_fontsize: float = 8,
    xticklabel_fontsize: float = 7,
    dendrogram_size_inches: float = 0.28,
    dendrogram_gap_inches: float = 0.02,
    side_label_space_inches: float = 0.8,
    bottom_space_inches: float = 0.8,
    top_space_inches: float = 0.5,
    cbar_size_inches: float = 0.14,
    cbar_gap_inches: float = 0.12,
    cbar_span_fraction: Optional[float] = None,
    cbar_orientation: str = "auto",
    cbar_label: str = "Number of samples with at least 1 interaction",
):
    """Plot a QC clustermap of per-condition sample counts.

    The count in each cell is the number of unique samples within a condition
    that have at least one row with ``value_col > presence_threshold`` for the
    requested feature.

    Parameters
    ----------
    orientation:
        ``"horizontal"`` plots conditions as rows and features as columns.
        ``"vertical"`` plots features as rows and conditions as columns.

    cbar_span_fraction:
        Optional fraction of the full figure span used by the colorbar along its
        major axis. Horizontal colorbars use figure width; vertical colorbars
        use figure height. When omitted, the colorbar matches the heatmap span.

    Returns
    -------
    g, plot_mat, count_df
        The seaborn ClusterGrid, the plotted matrix, and the long-form count
        table used to build the heatmap.
    """

    required = [value_col, condition_col, sample_col, feature_col]
    _validate_required_columns(df, required)

    if orientation not in {"horizontal", "vertical"}:
        raise ValueError("orientation must be one of: 'horizontal', 'vertical'")

    if cbar_orientation not in {"auto", "horizontal", "vertical"}:
        raise ValueError(
            "cbar_orientation must be one of: 'auto', 'horizontal', 'vertical'"
        )

    if cbar_span_fraction is not None and not (0 < cbar_span_fraction <= 1):
        raise ValueError("cbar_span_fraction must be within (0, 1].")

    plot_df = df[required].dropna().copy()
    plot_df[value_col] = pd.to_numeric(plot_df[value_col], errors="coerce")
    plot_df = plot_df.dropna(subset=[value_col]).copy()

    if plot_df.empty:
        raise ValueError("No valid rows remain after dropping missing values.")

    for column in [condition_col, sample_col, feature_col]:
        plot_df[column] = plot_df[column].astype(str)

    observed_conditions = plot_df[condition_col].drop_duplicates().tolist()
    observed_features = plot_df[feature_col].drop_duplicates().tolist()

    resolved_condition_order = _merge_order(condition_order, observed_conditions)
    resolved_feature_order = _merge_order(feature_order, observed_features)

    plot_df["present"] = plot_df[value_col].gt(presence_threshold).astype(int)

    sample_presence = (
        plot_df.groupby([condition_col, sample_col, feature_col], observed=True)["present"]
        .max()
        .reset_index()
    )

    if sample_presence.empty:
        raise ValueError("No rows remain after collapsing to sample-level presence.")

    count_df = (
        sample_presence.groupby([condition_col, feature_col], observed=True)["present"]
        .sum()
        .reset_index(name="sample_count")
    )

    mat = count_df.pivot_table(
        index=feature_col,
        columns=condition_col,
        values="sample_count",
        fill_value=0,
        aggfunc="sum",
    )

    mat = mat.reindex(index=resolved_feature_order, columns=resolved_condition_order, fill_value=0)
    mat = mat.astype(int)
    mat = mat.loc[mat.sum(axis=1) > 0, mat.sum(axis=0) > 0]

    if mat.empty:
        raise ValueError("The plotted matrix is empty after dropping zero rows and columns.")

    if orientation == "horizontal":
        plot_mat = mat.T.copy()
        resolved_xlabel = xlabel or feature_col.replace("_", " ").title()
        resolved_ylabel = ylabel or condition_col.replace("_", " ").title()
        auto_cbar_orientation = "horizontal"
    else:
        plot_mat = mat.copy()
        resolved_xlabel = xlabel or condition_col.replace("_", " ").title()
        resolved_ylabel = ylabel or feature_col.replace("_", " ").title()
        auto_cbar_orientation = "vertical"

    plot_mat.index = _wrap_labels(plot_mat.index.tolist(), row_label_wrap)
    plot_mat.columns = _wrap_labels(plot_mat.columns.tolist(), col_label_wrap)

    if figsize is None:
        figsize = _default_figsize(
            n_rows=plot_mat.shape[0],
            n_cols=plot_mat.shape[1],
            orientation=orientation,
        )

    row_label_chars = _max_label_line_length(plot_mat.index.tolist())
    col_label_chars = _max_label_line_length(plot_mat.columns.tolist())

    side_label_space_inches = max(side_label_space_inches, 0.08 * row_label_chars + 0.25)
    bottom_space_inches = max(bottom_space_inches, 0.09 * col_label_chars + 0.35)

    fig_w, fig_h = figsize
    row_dendro_ratio = min(dendrogram_size_inches / fig_w, 0.18)
    col_dendro_ratio = min(dendrogram_size_inches / fig_h, 0.18)

    if cluster_rows is None:
        cluster_rows = orientation == "vertical"
    if cluster_cols is None:
        cluster_cols = orientation == "horizontal"

    cluster_rows = bool(cluster_rows) and plot_mat.shape[0] > 1
    cluster_cols = bool(cluster_cols) and plot_mat.shape[1] > 1

    resolved_cbar_orientation = (
        auto_cbar_orientation if cbar_orientation == "auto" else cbar_orientation
    )

    annot_kws = {"fontsize": annot_fontsize}
    if annot_color is not None:
        annot_kws["color"] = annot_color

    g = sns.clustermap(
        plot_mat,
        cmap=cmap,
        figsize=figsize,
        xticklabels=True,
        yticklabels=True,
        row_cluster=cluster_rows,
        col_cluster=cluster_cols,
        method="average",
        metric="euclidean",
        linewidths=linewidths,
        linecolor=linecolor,
        annot=annot,
        fmt=fmt,
        annot_kws=annot_kws,
        dendrogram_ratio=(row_dendro_ratio, col_dendro_ratio),
        cbar_pos=(0.02, 0.02, 0.2, 0.02),
        cbar_kws={
            "orientation": resolved_cbar_orientation,
            "label": cbar_label,
        },
    )

    left_margin = min((side_label_space_inches / fig_w) if orientation == "vertical" else (0.3 / fig_w), 0.6)
    right_margin = min((side_label_space_inches / fig_w) if orientation == "horizontal" else (0.7 / fig_w), 0.35)
    bottom_margin = min(bottom_space_inches / fig_h, 0.55)
    top_margin = min(top_space_inches / fig_h, 0.22)

    g.fig.set_dpi(dpi if dpi is not None else g.fig.get_dpi())
    g.fig.subplots_adjust(
        left=left_margin,
        right=1 - right_margin,
        bottom=bottom_margin,
        top=1 - top_margin,
    )
    g.fig.canvas.draw()

    heatmap_pos = g.ax_heatmap.get_position()

    if cluster_rows:
        dendro_width = dendrogram_size_inches / fig_w
        dendro_gap_x = dendrogram_gap_inches / fig_w
        g.ax_row_dendrogram.set_position([
            heatmap_pos.x0 - dendro_width - dendro_gap_x,
            heatmap_pos.y0,
            dendro_width,
            heatmap_pos.height,
        ])
    else:
        g.ax_row_dendrogram.set_visible(False)

    if cluster_cols:
        dendro_height = dendrogram_size_inches / fig_h
        dendro_gap_y = dendrogram_gap_inches / fig_h
        g.ax_col_dendrogram.set_position([
            heatmap_pos.x0,
            heatmap_pos.y1 + dendro_gap_y,
            heatmap_pos.width,
            dendro_height,
        ])
    else:
        g.ax_col_dendrogram.set_visible(False)

    if resolved_cbar_orientation == "horizontal":
        cbar_height = cbar_size_inches / fig_h
        cbar_gap = cbar_gap_inches / fig_h
        cbar_width = heatmap_pos.width
        cbar_x0 = heatmap_pos.x0

        if cbar_span_fraction is not None:
            cbar_width = cbar_span_fraction
            heatmap_center_x = heatmap_pos.x0 + (heatmap_pos.width / 2)
            cbar_x0 = min(max(heatmap_center_x - (cbar_width / 2), 0.0), 1 - cbar_width)

        g.cax.set_position([
            cbar_x0,
            heatmap_pos.y0 - cbar_gap - cbar_height,
            cbar_width,
            cbar_height,
        ])
        g.cax.xaxis.set_label_position("bottom")
        g.cax.xaxis.tick_bottom()
        g.cax.tick_params(axis="x", labelsize=8, length=2, pad=2)
    else:
        cbar_width = cbar_size_inches / fig_w
        cbar_gap = cbar_gap_inches / fig_w
        cbar_height = heatmap_pos.height
        cbar_y0 = heatmap_pos.y0

        if cbar_span_fraction is not None:
            cbar_height = cbar_span_fraction
            heatmap_center_y = heatmap_pos.y0 + (heatmap_pos.height / 2)
            cbar_y0 = min(max(heatmap_center_y - (cbar_height / 2), 0.0), 1 - cbar_height)

        g.cax.set_position([
            heatmap_pos.x1 + cbar_gap,
            cbar_y0,
            cbar_width,
            cbar_height,
        ])
        g.cax.yaxis.set_label_position("right")
        g.cax.yaxis.tick_right()
        g.cax.tick_params(axis="y", labelsize=8, length=2, pad=2)

    g.ax_heatmap.set_xlabel(resolved_xlabel)
    g.ax_heatmap.set_ylabel(resolved_ylabel)
    g.ax_heatmap.set_aspect("auto")

    g.ax_heatmap.tick_params(axis="x", labelsize=xticklabel_fontsize, length=0, pad=4)
    g.ax_heatmap.tick_params(axis="y", labelsize=yticklabel_fontsize, length=0, pad=4)

    if orientation == "horizontal":
        g.ax_heatmap.yaxis.tick_right()
        g.ax_heatmap.yaxis.set_label_position("right")

        for label in g.ax_heatmap.get_xticklabels():
            label.set_rotation(90)
            label.set_horizontalalignment("center")
            label.set_verticalalignment("top")

        for label in g.ax_heatmap.get_yticklabels():
            label.set_rotation(0)
            label.set_horizontalalignment("left")
            label.set_verticalalignment("center")
    else:
        g.ax_heatmap.yaxis.tick_left()
        g.ax_heatmap.yaxis.set_label_position("left")

        for label in g.ax_heatmap.get_xticklabels():
            label.set_rotation(0)
            label.set_horizontalalignment("center")
            label.set_verticalalignment("top")

        for label in g.ax_heatmap.get_yticklabels():
            label.set_rotation(0)
            label.set_horizontalalignment("right")
            label.set_verticalalignment("center")

    if title:
        g.fig.suptitle(
            title,
            fontsize=12,
            fontweight="bold",
            y=1 - min(0.08 / fig_h, 0.04),
        )

    return g, plot_mat, count_df