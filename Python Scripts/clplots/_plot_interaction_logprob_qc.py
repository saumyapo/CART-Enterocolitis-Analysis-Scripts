from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable, Tuple, Union

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ResultsLike = Union[pd.DataFrame, str, Path]
InteractionSpec = Tuple[str, str, str]


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


def _prepare_plot_df(
    data: pd.DataFrame,
    source: str,
    target: str,
    interaction_name: str,
    source_col: str,
    target_col: str,
    interaction_col: str,
    logprob_col: str,
    condition_col: str,
    chemistry_col: str | None,
) -> pd.DataFrame:
    required = [
        source_col,
        target_col,
        interaction_col,
        logprob_col,
        condition_col,
    ]
    if chemistry_col is not None:
        required.append(chemistry_col)

    _validate_required_columns(data, required)

    plot_df = data.loc[
        data[source_col].eq(source)
        & data[target_col].eq(target)
        & data[interaction_col].eq(interaction_name),
        required,
    ].dropna().copy()

    if plot_df.empty:
        raise ValueError(
            "No rows remain after filtering for the requested source, target, and interaction_name."
        )

    plot_df[logprob_col] = pd.to_numeric(plot_df[logprob_col], errors="coerce")
    plot_df = plot_df.dropna(subset=[logprob_col]).copy()

    if plot_df.empty:
        raise ValueError("No valid logprob values remain after filtering.")

    return plot_df


def _draw_interaction_logprob_qc(
    ax,
    plot_df: pd.DataFrame,
    source: str,
    target: str,
    interaction_name: str,
    logprob_col: str,
    condition_col: str,
    chemistry_col: str | None,
    title: str | None,
    font_size: float,
    xlabel: str,
    ylabel: str,
    jitter: float | bool,
    dodge: bool,
    box_color: str,
    box_linewidth: float,
    strip_size: float,
    strip_alpha: float,
    strip_color: str,
    title_pad: float,
    show_legend: bool,
):
    sns.boxplot(
        data=plot_df,
        x=logprob_col,
        y=condition_col,
        color=box_color,
        linewidth=box_linewidth,
        fliersize=0,
        ax=ax,
    )
    stripplot_kwargs = {
        "data": plot_df,
        "x": logprob_col,
        "y": condition_col,
        "jitter": jitter,
        "size": strip_size,
        "alpha": strip_alpha,
        "ax": ax,
    }
    if chemistry_col is None:
        stripplot_kwargs["color"] = strip_color
        stripplot_kwargs["dodge"] = False
    else:
        stripplot_kwargs["hue"] = chemistry_col
        stripplot_kwargs["dodge"] = dodge

    sns.stripplot(**stripplot_kwargs)

    ax.set_xlabel(xlabel, fontsize=font_size)
    ax.set_ylabel(ylabel, fontsize=font_size)
    ax.set_title(
        title if title is not None else f"{source} -> {target}\n{interaction_name}",
        pad=title_pad,
        fontsize=font_size + 1,
    )
    ax.tick_params(axis="both", labelsize=font_size)
    ax.grid(axis="x", visible=False)
    ax.grid(axis="y", visible=False)
    sns.despine(ax=ax, trim=True, offset=5)

    legend = ax.get_legend()
    if legend is not None and not show_legend:
        legend.remove()

    return ax.get_legend_handles_labels()


def plot_interaction_logprob_qc(
    results: ResultsLike,
    source: str,
    target: str,
    interaction_name: str,
    source_col: str = "source",
    target_col: str = "target",
    interaction_col: str = "interaction_name",
    logprob_col: str = "logprob",
    condition_col: str = "condition",
    chemistry_col: str = "chemistry",
    figsize: Tuple[float, float] = (4, 1.15),
    dpi=None,
    title: str | None = None,
    font_size: float = 6,
    xlabel: str = "log(CellChat Interaction Probability)",
    ylabel: str = "Condition",
    jitter: float | bool = True,
    dodge: bool = True,
    box_color: str = "#E6E6E6",
    box_linewidth: float = 1.0,
    strip_size: float = 3.0,
    strip_alpha: float = 0.9,
    title_pad: float = 4.0,
    legend_anchor_x: float = 1.0,
    legend_loc: str = "center left",
):
    """Plot raw logprob QC for one interaction across condition and chemistry.

    Returns
    -------
    fig, ax, plot_df
        Figure, axis, and the filtered dataframe used for plotting.
    """

    data = _load_results(results)

    plot_df = _prepare_plot_df(
        data=data,
        source=source,
        target=target,
        interaction_name=interaction_name,
        source_col=source_col,
        target_col=target_col,
        interaction_col=interaction_col,
        logprob_col=logprob_col,
        condition_col=condition_col,
        chemistry_col=chemistry_col,
    )

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi, constrained_layout=True)

    _draw_interaction_logprob_qc(
        ax=ax,
        plot_df=plot_df,
        source=source,
        target=target,
        interaction_name=interaction_name,
        logprob_col=logprob_col,
        condition_col=condition_col,
        chemistry_col=chemistry_col,
        title=title,
        font_size=font_size,
        xlabel=xlabel,
        ylabel=ylabel,
        jitter=jitter,
        dodge=dodge,
        box_color=box_color,
        box_linewidth=box_linewidth,
        strip_size=strip_size,
        strip_alpha=strip_alpha,
        strip_color="#111111",
        title_pad=title_pad,
        show_legend=True,
    )

    legend = ax.get_legend()
    if legend is not None:
        sns.move_legend(
            ax,
            legend_loc,
            bbox_to_anchor=(legend_anchor_x, 0.5),
            frameon=False,
            title=chemistry_col,
            borderaxespad=0.0,
        )
        legend = ax.get_legend()
        plt.setp(legend.get_title(), fontsize=font_size)
        plt.setp(legend.get_texts(), fontsize=font_size)

    return fig, ax, plot_df


def plot_interaction_logprob_qc_stack(
    results: ResultsLike,
    interactions: Iterable[InteractionSpec],
    max_stack: int = 3,
    source_col: str = "source",
    target_col: str = "target",
    interaction_col: str = "interaction_name",
    logprob_col: str = "logprob",
    condition_col: str = "condition",
    panel_figsize: Tuple[float, float] = (3.75, 1.15),
    figsize: Tuple[float, float] | None = None,
    dpi=None,
    font_size: float = 6,
    xlabel: str = "log(CellChat Interaction Probability)",
    ylabel: str = "Condition",
    jitter: float | bool = True,
    dodge: bool = True,
    box_color: str = "#E6E6E6",
    box_linewidth: float = 1.0,
    strip_size: float = 3.0,
    strip_alpha: float = 0.9,
    title_pad: float = 2.5,
    xpad_fraction: float = 0.03,
    left_margin: float = 0.12,
    right_margin: float = 0.96,
    top_margin: float = 0.96,
    bottom_margin: float = 0.12,
    hspace: float = 0.06,
    wspace: float = 0.12,
):
    """Plot multiple interaction QC panels with uniform black sample points."""

    interaction_list = list(interactions)
    if not interaction_list:
        raise ValueError("interactions must contain at least one (source, target, interaction_name) tuple")
    if max_stack < 1:
        raise ValueError("max_stack must be >= 1")

    data = _load_results(results)

    prepared = []
    for source, target, interaction_name in interaction_list:
        plot_df = _prepare_plot_df(
            data=data,
            source=source,
            target=target,
            interaction_name=interaction_name,
            source_col=source_col,
            target_col=target_col,
            interaction_col=interaction_col,
            logprob_col=logprob_col,
            condition_col=condition_col,
            chemistry_col=None,
        )
        prepared.append(
            {
                "source": source,
                "target": target,
                "interaction_name": interaction_name,
                "plot_df": plot_df,
            }
        )

    n_panels = len(prepared)
    n_rows = min(max_stack, n_panels)
    n_cols = math.ceil(n_panels / max_stack)

    if figsize is None:
        figsize = (
            panel_figsize[0] * n_cols,
            panel_figsize[1] * n_rows,
        )

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=figsize,
        dpi=dpi,
        sharex=False,
        constrained_layout=False,
        squeeze=False,
    )

    combined_data = []

    for index, item in enumerate(prepared):
        row = index % max_stack
        col = index // max_stack
        ax = axes[row, col]

        _draw_interaction_logprob_qc(
            ax=ax,
            plot_df=item["plot_df"],
            source=item["source"],
            target=item["target"],
            interaction_name=item["interaction_name"],
            logprob_col=logprob_col,
            condition_col=condition_col,
            chemistry_col=None,
            title=None,
            font_size=font_size,
            xlabel=xlabel,
            ylabel=ylabel,
            jitter=jitter,
            dodge=dodge,
            box_color=box_color,
            box_linewidth=box_linewidth,
            strip_size=strip_size,
            strip_alpha=strip_alpha,
            strip_color="#111111",
            title_pad=title_pad,
            show_legend=False,
        )

        panel_x_min = float(item["plot_df"][logprob_col].min())
        panel_x_max = float(item["plot_df"][logprob_col].max())
        panel_x_pad = (panel_x_max - panel_x_min) * xpad_fraction
        if panel_x_pad == 0:
            panel_x_pad = max(abs(panel_x_min) * xpad_fraction, 0.1)
        ax.set_xlim(panel_x_min - panel_x_pad, panel_x_max + panel_x_pad)

        if col != 0:
            ax.set_ylabel("")

        panel_df = item["plot_df"].copy()
        panel_df["panel_source"] = item["source"]
        panel_df["panel_target"] = item["target"]
        panel_df["panel_interaction_name"] = item["interaction_name"]
        combined_data.append(panel_df)

    for index in range(n_panels, n_rows * n_cols):
        row = index % max_stack
        col = index // max_stack
        axes[row, col].set_visible(False)

    for col in range(n_cols):
        visible_rows = [row for row in range(n_rows) if axes[row, col].get_visible()]
        if not visible_rows:
            continue

        bottom_row = max(visible_rows)
        for row in visible_rows:
            ax = axes[row, col]
            if row == bottom_row:
                ax.set_xlabel(xlabel, fontsize=font_size)
            else:
                ax.set_xlabel("")
                ax.tick_params(axis="x", labelbottom=True)

    fig.subplots_adjust(
        left=left_margin,
        right=right_margin,
        top=top_margin,
        bottom=bottom_margin,
        hspace=hspace,
        wspace=wspace,
    )

    combined_df = pd.concat(combined_data, axis=0, ignore_index=True)
    return fig, axes, combined_df