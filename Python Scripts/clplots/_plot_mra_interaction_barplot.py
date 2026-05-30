from __future__ import annotations

import inspect
import textwrap
from pathlib import Path
from typing import Optional, Sequence, Tuple, Union

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


def _wrap_label(label: str, width: int) -> str:
    if width <= 0:
        return label
    return textwrap.fill(
        str(label),
        width=width,
        break_long_words=False,
        break_on_hyphens=False,
    )


def _format_target_label(label: str) -> str:
    formatted = str(label).replace("_CARpos", " (CAR+)")
    return " ".join(formatted.split())


def _format_interaction_label(label: str) -> str:
    return str(label).replace("_", "-")


def _resolve_figsize(figsize: Optional[Tuple[float, float]], n_interactions: int) -> Tuple[float, float]:
    if figsize is not None:
        return figsize
    return (8.0, max(4.0, 0.32 * n_interactions + 1.4))


def plot_mra_interaction_barplot(
    results: ResultsLike,
    target_pattern: Optional[str] = "Cytotox",
    target_values: Optional[Sequence[str]] = None,
    target_col: str = "target",
    interaction_col: str = "interaction_name",
    rank_col: str = "rank_median",
    min_rank: float = 0.9,
    rank_agg: str = "median",
    interaction_order_metric: str = "max",
    top_n: Optional[int] = 30,
    target_regex: bool = True,
    target_case: bool = False,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[float] = None,
    title: Optional[str] = "High-ranking interactions across cytotoxic targets",
    xlabel: str = "Median rank percentile",
    ylabel: str = "Ligand-receptor interaction",
    legend_title: str = "Target",
    legend_loc: str = "upper right",
    legend_bbox_to_anchor: Optional[Tuple[float, float]] = None,
    palette: Union[str, Sequence[str]] = "colorblind",
    xlim: Optional[Tuple[float, float]] = None,
    interaction_label_wrap_width: int = 28,
    title_fontsize: float = 13,
    axis_labelsize: float = 11,
    tick_labelsize: float = 9,
    legend_fontsize: float = 9,
    legend_title_fontsize: float = 10,
    linewidth: float = 0.6,
    edgecolor: str = "#2F2F2F",
    saturation: float = 0.95,
    ax=None,
):
    """Plots grouped barplot for median-rank interaction results."""

    results_df = _load_results(results)

    required = [target_col, interaction_col, rank_col]
    missing = [column for column in required if column not in results_df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if target_values is not None and target_pattern is not None:
        raise ValueError("Provide either target_values or target_pattern, not both.")

    if rank_agg not in {"mean", "median", "max", "min"}:
        raise ValueError("rank_agg must be one of: mean, median, max, min.")

    if interaction_order_metric not in {"mean", "median", "max", "min"}:
        raise ValueError("interaction_order_metric must be one of: mean, median, max, min.")

    plot_df = results_df[[target_col, interaction_col, rank_col]].dropna().copy()
    plot_df[rank_col] = pd.to_numeric(plot_df[rank_col], errors="coerce")
    plot_df = plot_df.dropna(subset=[rank_col]).copy()

    if target_values is not None:
        plot_df = plot_df[plot_df[target_col].isin(target_values)].copy()
    elif target_pattern is not None:
        plot_df = plot_df[
            plot_df[target_col].astype(str).str.contains(
                target_pattern,
                regex=target_regex,
                case=target_case,
                na=False,
            )
        ].copy()

    plot_df = plot_df[plot_df[rank_col] >= min_rank].copy()

    if plot_df.empty:
        raise ValueError("No rows remain after applying the requested target and rank filters.")

    aggregated = (
        plot_df.groupby([target_col, interaction_col], observed=True)[rank_col]
        .agg(rank_agg)
        .reset_index()
        .rename(columns={rank_col: "rank_value"})
    )

    interaction_order = (
        aggregated.groupby(interaction_col, observed=True)["rank_value"]
        .agg(interaction_order_metric)
        .sort_values(ascending=False)
        .index.tolist()
    )

    if top_n is not None:
        interaction_order = interaction_order[:top_n]
        aggregated = aggregated[aggregated[interaction_col].isin(interaction_order)].copy()

    target_order = (
        aggregated.groupby(target_col, observed=True)["rank_value"]
        .mean()
        .sort_values(ascending=False)
        .index.tolist()
    )

    aggregated["interaction_label"] = aggregated[interaction_col].map(
        lambda label: _wrap_label(_format_interaction_label(label), interaction_label_wrap_width)
    )
    interaction_label_map = (
        aggregated[[interaction_col, "interaction_label"]]
        .drop_duplicates()
        .set_index(interaction_col)["interaction_label"]
        .to_dict()
    )
    label_order = [interaction_label_map[label] for label in interaction_order if label in interaction_label_map]

    aggregated[target_col] = pd.Categorical(aggregated[target_col], categories=target_order, ordered=True)
    aggregated["interaction_label"] = pd.Categorical(
        aggregated["interaction_label"],
        categories=label_order,
        ordered=True,
    )
    aggregated = aggregated.sort_values(["interaction_label", target_col]).reset_index(drop=True)

    if ax is None:
        fig, ax = plt.subplots(
            figsize=_resolve_figsize(figsize=figsize, n_interactions=len(label_order)),
            dpi=dpi,
            constrained_layout=True,
        )
    else:
        fig = ax.figure

    with sns.axes_style(
        "whitegrid",
        rc={
            "axes.spines.right": False,
            "axes.spines.top": False,
            "grid.color": "#D9D9D9",
            "grid.linewidth": 0.6,
            "grid.alpha": 0.3,
        },
    ):
        barplot_kwargs = {
            "data": aggregated,
            "y": "interaction_label",
            "x": "rank_value",
            "hue": target_col,
            "order": label_order,
            "hue_order": target_order,
            "palette": palette,
            "edgecolor": edgecolor,
            "linewidth": linewidth,
            "saturation": saturation,
            "ax": ax,
        }

        if "errorbar" in inspect.signature(sns.barplot).parameters:
            barplot_kwargs["errorbar"] = None
        else:
            barplot_kwargs["ci"] = None

        sns.barplot(**barplot_kwargs)

    ax.set_xlabel(xlabel, fontsize=axis_labelsize)
    ax.set_ylabel(ylabel, fontsize=axis_labelsize)
    ax.tick_params(axis="x", labelsize=tick_labelsize)
    ax.tick_params(axis="y", labelsize=tick_labelsize)

    if title is not None:
        ax.set_title(title, fontsize=title_fontsize, fontweight="bold", pad=10)

    if xlim is None:
        xlim = (min_rank, 1.0)
    ax.set_xlim(*xlim)

    ax.grid(axis="x", alpha=0.3, linewidth=0.8)
    ax.grid(axis="y", visible=False)
    sns.despine(ax=ax, trim=True, offset=5)

    legend_kwargs = {
        "title": legend_title,
        "frameon": False,
        "loc": legend_loc,
    }
    if legend_bbox_to_anchor is not None:
        legend_kwargs["bbox_to_anchor"] = legend_bbox_to_anchor

    legend = ax.legend(**legend_kwargs)
    if legend is not None:
        legend.set_title(legend_title)
        legend.get_title().set_fontsize(legend_title_fontsize)
        for text, label in zip(legend.get_texts(), target_order):
            text.set_text(_format_target_label(label))
            text.set_fontsize(legend_fontsize)

    return fig, ax, aggregated