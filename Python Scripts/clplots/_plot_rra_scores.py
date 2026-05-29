from __future__ import annotations

import textwrap
from pathlib import Path
from typing import Mapping, Optional, Sequence, Tuple, Union

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist


ResultsLike = Union[pd.DataFrame, str, Path]

DEFAULT_TARGETS = (
	"Cytotoxic cells",
	"Cytotoxic cells_CARpos",
)

DEFAULT_RRA_SOURCE_GROUPS = {
	"Fibroblasts": (
		"Inflam. fibroblasts",
		"Crypt associated fibroblasts",
		"Tissue fibroblasts",
		"Activated fibroblasts",
		"Myofibroblasts",
	),
	"Endothelial cells": (
		"Endothelial cells",
		"Post-capillary venous endothelium",
		"Capillary endothelium",
		"Arterial endothelium",
		"Arteriolar endothelium",
		"Lymphatic endothelium",
		"Vascular endothelium",
		"Activated endothelium",
		"Angiogenic endothelium",
	),
	"Epithelial cells": (
		"Crypt progenitor",
		"Goblet cells",
		"TA cells",
		"Absorptive epithelium",
		"BEST4+ colonocytes",
		"Mature colonocytes",
		"Enteroendocrine cells",
		"Tuft cells",
	),
	"Monocytes / macrophages": (
		"Resident macrophages",
		"Inflam. monocytes",
		"Inflam. macrophages",
		"Inflam. monocytes/cDC2",
		"Macrophages",
		"Myeloid cells",
	),
	"Dendritic cells": (
		"cDC2",
		"cDC1",
		"LAMP3+ DC",
	),
	"T cells": (
		"Th17",
		"Cytotoxic cells",
		"Naive/effector",
		"CD8+ naive/central memory",
		"CD4+ effector memory",
		"CD4/Tfh",
		"CD8+ Trm",
		"Tregs",
		"Naive/central memory",
	),
	"CAR+ T cells": (
		"Cytotoxic cells_CARpos",
		"Naive/effector_CARpos",
		"CD8+ naive/central memory_CARpos",
		"CD4+ effector memory_CARpos",
	),
	"NK / innate lymphoid": (
		"NK-like",
		"Cytotoxic/NK-like",
		"ILC-like",
		"MAIT-like",
		"γδT",
	),
}


def _build_source_group_lookup(
	source_groups: Optional[Mapping[str, Sequence[str]]],
) -> dict[str, str]:
	"""Build a source-to-group lookup and reject ambiguous assignments."""
	if source_groups is None:
		return {}

	lookup: dict[str, str] = {}
	for group_name, members in source_groups.items():
		for member in members:
			member_label = str(member)
			if member_label in lookup and lookup[member_label] != str(group_name):
				raise ValueError(
					"Each source cell type can only belong to one source group. "
					f"'{member_label}' was assigned to both "
					f"'{lookup[member_label]}' and '{group_name}'."
				)
			lookup[member_label] = str(group_name)

	return lookup


def _target_tick_labelsize(n_targets: int) -> int:
	if n_targets <= 6:
		return 11
	if n_targets <= 12:
		return 9
	if n_targets <= 20:
		return 8
	if n_targets <= 30:
		return 7
	return 6


def _format_percent_legend(ax: plt.Axes) -> None:
	legend = ax.get_legend()
	if legend is None:
		return

	current_section = None
	for text in legend.get_texts():
		label = text.get_text().strip()
		if label == "mean_rank_pct":
			current_section = "mean_rank_pct"
			text.set_text("mean_rank_pct")
			continue
		if label == "top10_freq":
			current_section = "top10_freq"
			text.set_text("\ntop10_freq")
			continue

		if current_section in {"mean_rank_pct", "top10_freq"}:
			try:
				value = float(label)
			except ValueError:
				continue
			text.set_text(f"{value:.0%}")


def _format_percent_values_uniquely(values: Sequence[float], max_precision: int = 3) -> list[str]:
	"""Format percentages with enough precision to keep legend labels distinct."""
	if not values:
		return []

	for precision in range(max_precision + 1):
		formatted = [f"{value:.{precision}%}" for value in values]
		if len(set(formatted)) == len(formatted):
			return formatted

	return [f"{value:.{max_precision}%}" for value in values]


def _format_publication_percent_legend(
	ax: plt.Axes,
	hue_label: str = "Mean Rank Percentile",
	size_label: str = "Frequency Among Top 10 Rankings",
	hue_label_wrap_width: int = 0,
	size_label_wrap_width: int = 0,
	hue_title_fontsize: Optional[Union[str, float]] = None,
	size_title_fontsize: Optional[Union[str, float]] = None,
	hue_value_fontsize: Optional[Union[str, float]] = None,
	size_value_fontsize: Optional[Union[str, float]] = None,
) -> None:
	legend = ax.get_legend()
	if legend is None:
		return

	wrapped_hue_label = _wrap_label(hue_label, hue_label_wrap_width)
	wrapped_size_label = _wrap_label(size_label, size_label_wrap_width)

	hue_texts: list[plt.Text] = []
	hue_values: list[float] = []
	size_texts: list[plt.Text] = []
	size_values: list[float] = []

	current_section = None
	for text in legend.get_texts():
		label = text.get_text().strip()
		if label == "mean_rank_pct":
			current_section = "mean_rank_pct"
			text.set_text(wrapped_hue_label)
			text.set_fontweight("semibold")
			if hue_title_fontsize is not None:
				text.set_fontsize(hue_title_fontsize)
			continue
		if label == "top10_freq":
			current_section = "top10_freq"
			text.set_text(f"\n{wrapped_size_label}")
			text.set_fontweight("semibold")
			if size_title_fontsize is not None:
				text.set_fontsize(size_title_fontsize)
			continue

		if current_section in {"mean_rank_pct", "top10_freq"}:
			try:
				value = float(label)
			except ValueError:
				continue
			if current_section == "mean_rank_pct" and hue_value_fontsize is not None:
				text.set_fontsize(hue_value_fontsize)
			if current_section == "mean_rank_pct":
				hue_texts.append(text)
				hue_values.append(value)
			if current_section == "top10_freq" and size_value_fontsize is not None:
				text.set_fontsize(size_value_fontsize)
			if current_section == "top10_freq":
				size_texts.append(text)
				size_values.append(value)

	for text, label in zip(hue_texts, _format_percent_values_uniquely(hue_values)):
		text.set_text(label)

	for text, label in zip(size_texts, _format_percent_values_uniquely(size_values)):
		text.set_text(label)


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


def _prepare_rra_scores_plot_data(
	results: ResultsLike,
	target_values: Optional[Sequence[str]],
	fdr_thresh: float,
	top5_freq_thresh: float,
	top_n: int,
	x_axis_mode: str,
	log_prefix: str,
	mean_rank_pct_thresh: Optional[float] = None,
	source_groups: Optional[Mapping[str, Sequence[str]]] = None,
) -> Tuple[pd.DataFrame, str, str, Sequence[str], Sequence[str]]:
	results_df = _load_results(results)
	if x_axis_mode not in {"target", "source_target"}:
		raise ValueError("x_axis_mode must be either 'target' or 'source_target'")

	print(f"{log_prefix}: loaded {len(results_df)} rows")

	if target_values is None:
		plot_df = results_df.copy()
		print(f"{log_prefix}: target filter skipped, {len(plot_df)} rows remain")
	else:
		plot_df = results_df[results_df["target"].isin(target_values)].copy()
		print(
			f"{log_prefix}: after target filter "
			f"({len(target_values)} targets), {len(plot_df)} rows remain"
		)

	plot_df = plot_df[plot_df["rra_fdr"] < fdr_thresh].copy()
	print(f"{log_prefix}: after rra_fdr < {fdr_thresh}, {len(plot_df)} rows remain")

	plot_df = plot_df[plot_df["top5_freq"] >= top5_freq_thresh].copy()
	print(
		f"{log_prefix}: after top5_freq >= {top5_freq_thresh}, "
		f"{len(plot_df)} rows remain"
	)

	if source_groups is not None and x_axis_mode == "source_target":
		source_lookup = _build_source_group_lookup(source_groups)
		plot_df["source"] = (
			plot_df["source"].astype(str).map(source_lookup).fillna(plot_df["source"].astype(str))
		)
		print(
			f"{log_prefix}: collapsed sources into "
			f"{plot_df['source'].nunique()} grouped labels"
		)

	if plot_df.empty:
		raise ValueError("No rows remain after applying target and threshold filters.")

	group_cols = ["target", "interaction_name"]
	if x_axis_mode == "source_target":
		group_cols = ["source", "target", "interaction_name"]

	target_df = (
		plot_df.groupby(group_cols, observed=True)
		.agg(
			mean_rank_pct=("mean_rank_pct", "mean"),
			median_rank_pct=("median_rank_pct", "mean"),
			top10_freq=("top10_freq", "mean"),
			top5_freq=("top5_freq", "mean"),
			best_rra_fdr=("rra_fdr", "min"),
			median_logprob=("median_logprob", "mean"),
			n_sources=("source", "nunique"),
		)
		.reset_index()
	)
	print(
		f"{log_prefix}: after aggregation, {len(target_df)} grouped rows across "
		f"{target_df['interaction_name'].nunique()} unique interactions"
	)

	if mean_rank_pct_thresh is not None:
		target_df = target_df[target_df["mean_rank_pct"] >= mean_rank_pct_thresh].copy()
		print(
			f"{log_prefix}: after mean_rank_pct >= {mean_rank_pct_thresh:.0%}, "
			f"{len(target_df)} grouped rows remain"
		)
		if target_df.empty:
			raise ValueError(
				"No rows remain after applying the mean_rank_pct threshold."
			)

	top_interactions = (
		target_df.groupby("interaction_name", observed=True)
		.agg(
			mean_rank_pct=("mean_rank_pct", "mean"),
			best_rra_fdr=("best_rra_fdr", "min"),
		)
		.sort_values(["best_rra_fdr", "mean_rank_pct"], ascending=[True, False])
		.head(top_n)
		.index
	)
	print(
		f"{log_prefix}: keeping top {min(top_n, len(top_interactions))} "
		"unique interactions after ranking"
	)

	target_df = target_df[target_df["interaction_name"].isin(top_interactions)].copy()
	print(
		f"{log_prefix}: after top_n interaction filter, {len(target_df)} grouped rows remain"
	)
	print(
		f"{log_prefix}: final unique interactions remaining = "
		f"{target_df['interaction_name'].nunique()}"
	)
	if target_df.empty:
		raise ValueError("No rows remain after selecting top interactions.")

	if x_axis_mode == "source_target":
		target_df["source_target_pair"] = (
			target_df["source"].astype(str) + " → " + target_df["target"].astype(str)
		)
		x_axis_col = "source_target_pair"
		x_axis_label = "Source → Target"
	else:
		x_axis_col = "target"
		x_axis_label = "Target cell type"

	mat = target_df.pivot_table(
		index="interaction_name",
		columns=x_axis_col,
		values="mean_rank_pct",
		aggfunc="mean",
		fill_value=0,
	)

	if len(mat.index) > 1:
		row_linkage = linkage(pdist(mat.values), method="average")
		row_order = mat.index[leaves_list(row_linkage)].tolist()
	else:
		row_order = mat.index.tolist()

	if len(mat.columns) > 1:
		col_linkage = linkage(pdist(mat.values.T), method="average")
		col_order = mat.columns[leaves_list(col_linkage)].tolist()
	else:
		col_order = mat.columns.tolist()

	target_df["interaction_name"] = pd.Categorical(
		target_df["interaction_name"],
		categories=row_order,
		ordered=True,
	)
	target_df[x_axis_col] = pd.Categorical(
		target_df[x_axis_col],
		categories=col_order,
		ordered=True,
	)

	return target_df, x_axis_col, x_axis_label, col_order, row_order


def _publication_tick_labelsize(n_labels: int, x_axis_mode: str) -> float:
	if x_axis_mode == "source_target":
		if n_labels <= 12:
			return 9
		if n_labels <= 24:
			return 8
		if n_labels <= 40:
			return 7
		return 6
	return _target_tick_labelsize(n_labels)


def _default_publication_figsize(
	n_x_labels: int,
	n_interactions: int,
	x_axis_mode: str,
) -> Tuple[float, float]:
	if x_axis_mode == "source_target":
		width = max(9.5, 0.55 * n_interactions + 4.0)
		height = max(10.0, 0.22 * n_x_labels + 2.8)
	else:
		width = max(8.5, 1.15 * n_x_labels + 4.0)
		height = max(6.0, 0.38 * n_interactions + 1.4)
	return (width, height)


def _format_publication_cell_label(label: str) -> str:
	formatted = str(label).replace("_CARpos", " (CAR+)")
	formatted = formatted.replace("_", " ")
	return " ".join(formatted.split())


def _wrap_label(label: str, width: int) -> str:
	if width <= 0:
		return label
	return textwrap.fill(
		label,
		width=width,
		break_long_words=False,
		break_on_hyphens=False,
	)


def _build_publication_x_tick_labels(
	target_df: pd.DataFrame,
	target_order: Sequence[str],
	x_axis_mode: str,
	x_label_wrap_width: int,
) -> Sequence[str]:
	if x_axis_mode == "source_target":
		label_frame = (
			target_df[["source_target_pair", "source", "target"]]
			.drop_duplicates()
			.set_index("source_target_pair")
		)
		return [
			f"{_wrap_label(_format_publication_cell_label(label_frame.at[pair, 'source']), x_label_wrap_width)}\n→\n"
			f"{_wrap_label(_format_publication_cell_label(label_frame.at[pair, 'target']), x_label_wrap_width)}"
			for pair in target_order
		]
	return [
		_wrap_label(_format_publication_cell_label(label), x_label_wrap_width)
		for label in target_order
	]


def _build_publication_pair_tick_labels(
	target_df: pd.DataFrame,
	pair_order: Sequence[str],
	wrap_width: int,
) -> Sequence[str]:
	label_frame = (
		target_df[["source_target_pair", "source", "target"]]
		.drop_duplicates()
		.set_index("source_target_pair")
	)
	labels = []
	for pair in pair_order:
		source_label = _format_publication_cell_label(label_frame.at[pair, "source"])
		target_label = _format_publication_cell_label(label_frame.at[pair, "target"])
		labels.append(_wrap_label(f"{source_label} → {target_label}", wrap_width))
	return labels


def _publication_y_tick_labelsize(n_labels: int, x_axis_mode: str) -> float:
	if x_axis_mode == "source_target":
		if n_labels <= 18:
			return 8.0
		if n_labels <= 36:
			return 7.0
		if n_labels <= 60:
			return 6.0
		if n_labels <= 90:
			return 5.5
		return 5.5
	return 9


def _scale_marker_area(
	value: float,
	value_min: float,
	value_max: float,
	size_range: Tuple[float, float],
) -> float:
	if abs(value_max - value_min) < 1e-12:
		return float(sum(size_range) / 2.0)
	scaled = (value - value_min) / (value_max - value_min)
	return size_range[0] + scaled * (size_range[1] - size_range[0])


def _legend_values(series: pd.Series, n_values: int) -> Sequence[float]:
	series = series.dropna()
	if series.empty:
		return []
	value_min = float(series.min())
	value_max = float(series.max())
	if abs(value_max - value_min) < 1e-12:
		return [value_min]
	ticks = [value_min + (value_max - value_min) * step / (n_values - 1) for step in range(n_values)]
	return sorted({round(tick, 6) for tick in ticks})


def plot_rra_scores(
	results: ResultsLike,
	target_values: Optional[Sequence[str]] = DEFAULT_TARGETS,
	fdr_thresh: float = 0.05,
	top5_freq_thresh: float = 0.5,
	top_n: int = 60,
	figsize: Tuple[float, float] = (8.5, 12),
	title: str = "Reproducible signaling programs aggregated by target",
	x_axis_mode: str = "target",
	x_tick_labelsize: Optional[Union[str, float]] = "x-small",
	mean_rank_pct_thresh: Optional[float] = None,
	source_groups: Optional[Mapping[str, Sequence[str]]] = None,
):
	"""Plot RRA results with either target or source-target x-axis labels.

	Args:
		results: Result dataframe or path to a saved RRA result table.
		target_values: Optional targets to include on the x-axis. When ``None``,
			no target filtering is applied.
		fdr_thresh: Maximum ``rra_fdr`` retained before collapsing sources.
		top5_freq_thresh: Minimum ``top5_freq`` retained before collapsing.
		top_n: Number of unique ``interaction_name`` values to keep after
			ranking. This does not cap the final number of plotted rows, because
			each retained interaction can still appear across multiple targets or
			multiple ``source -> target`` combinations.
		figsize: Figure size passed to ``plt.subplots``.
		title: Plot title.
		x_axis_mode: Either ``"target"`` to collapse across sources within each
			target, or ``"source_target"`` to keep full ``source -> target`` labels
			on the x-axis.
		x_tick_labelsize: Font size for x tick labels. If ``None``, choose a size
			based on the number of targets.
		source_groups: Optional mapping from group label to source cell types.
			When provided with ``x_axis_mode="source_target"``, listed sources are
			relabeled to the supplied group name before aggregation. Unlisted
			sources keep their original labels.

	Returns:
		``(fig, ax, target_df)`` where ``target_df`` is the dataframe used for
		plotting.
	"""
	target_df, x_axis_col, x_axis_label, target_order, _ = _prepare_rra_scores_plot_data(
		results=results,
		target_values=target_values,
		fdr_thresh=fdr_thresh,
		top5_freq_thresh=top5_freq_thresh,
		top_n=top_n,
		x_axis_mode=x_axis_mode,
		log_prefix="plot_rra_scores",
		mean_rank_pct_thresh=mean_rank_pct_thresh,
		source_groups=source_groups,
	)

	with sns.axes_style(
		"whitegrid",
		rc={
			"axes.spines.right": False,
			"axes.spines.top": False,
			"grid.color": "#d9d9d9",
			"grid.linewidth": 0.6,
			"grid.alpha": 0.35,
		},
	), sns.plotting_context("talk"):
		fig, ax = plt.subplots(figsize=figsize)
		resolved_x_tick_labelsize = x_tick_labelsize
		if resolved_x_tick_labelsize is None:
			resolved_x_tick_labelsize = _target_tick_labelsize(len(target_order))

		sns.scatterplot(
			data=target_df,
			x=x_axis_col,
			y="interaction_name",
			size="top10_freq",
			hue="mean_rank_pct",
			palette="viridis",
			sizes=(40, 240),
			edgecolor="white",
			linewidth=0.6,
			alpha=0.95,
			ax=ax,
		)

		ax.set_xlabel(x_axis_label, labelpad=12)
		ax.set_ylabel("Ligand-receptor / pathway", labelpad=10)
		ax.set_title(title, pad=16, weight="semibold")

		ax.tick_params(axis="x", rotation=90, labelsize=resolved_x_tick_labelsize)
		for label in ax.get_xticklabels():
			label.set_horizontalalignment("center")
			label.set_verticalalignment("top")

		ax.tick_params(axis="y", labelsize=9)
		ax.grid(axis="x", visible=False)
		ax.grid(axis="y", color="#d9d9d9", linewidth=0.5, alpha=0.25)

		sns.despine(ax=ax, offset=5, trim=True)
		sns.move_legend(
			ax,
			"center left",
			bbox_to_anchor=(1.02, 0.5),
			frameon=False,
			title=None,
		)
		_format_percent_legend(ax)

		plt.subplots_adjust(right=0.78)
		fig.tight_layout()

	return fig, ax, target_df


def plot_rra_scores_publication(
	results: ResultsLike,
	target_values: Optional[Sequence[str]] = DEFAULT_TARGETS,
	fdr_thresh: float = 0.05,
	top5_freq_thresh: float = 0.5,
	mean_rank_pct_thresh: Optional[float] = None,
	top_n: int = 60,
	figsize: Optional[Tuple[float, float]] = None,
	title: str = "Highest-ranking reproducible signaling interactions",
	x_axis_mode: str = "target",
	x_axis_label: Optional[str] = None,
	y_axis_label: Optional[str] = None,
	x_tick_labelsize: Optional[Union[str, float]] = None,
	palette: str = "viridis",
	size_range: Tuple[float, float] = (36, 220),
	x_label_wrap_width: int = 16,
	y_label_wrap_width: int = 0,
	hue_legend_label: str = "Mean Rank Percentile",
	size_legend_label: str = "Frequency Among Top 10 Rankings",
	hue_legend_wrap_width: int = 0,
	size_legend_wrap_width: int = 0,
	hue_legend_title_fontsize: Optional[Union[str, float]] = None,
	size_legend_title_fontsize: Optional[Union[str, float]] = None,
	hue_legend_value_fontsize: Optional[Union[str, float]] = None,
	size_legend_value_fontsize: Optional[Union[str, float]] = None,
	source_groups: Optional[Mapping[str, Sequence[str]]] = None,
):
	"""Plot RRA results

	This keeps the raw values in the returned dataframe, but formats the
	axes and legend text for cleaner presentation.

	Legend controls are exposed separately for the mean-rank and top-10
	frequency sections so notebook code can tune wrapping and font sizes
	without editing the helper.

	When ``source_groups`` is provided with ``x_axis_mode="source_target"``, the
	plot relabels matching sources to the supplied group names before averaging.
	This reduces clutter by collapsing related source cell types into a shared
	``group -> target`` visualization while leaving unmatched sources unchanged.

	Returns:
		``(fig, ax, target_df)`` where ``target_df`` is the dataframe used for
		plotting.
	"""
	target_df, x_axis_col, _, target_order, row_order = _prepare_rra_scores_plot_data(
		results=results,
		target_values=target_values,
		fdr_thresh=fdr_thresh,
		top5_freq_thresh=top5_freq_thresh,
		top_n=top_n,
		x_axis_mode=x_axis_mode,
		log_prefix="plot_rra_scores_publication",
		mean_rank_pct_thresh=mean_rank_pct_thresh,
		source_groups=source_groups,
	)

	if figsize is None:
		figsize = _default_publication_figsize(
			n_x_labels=len(target_order),
			n_interactions=len(row_order),
			x_axis_mode=x_axis_mode,
		)

	if x_axis_mode == "source_target":
		plot_x_col = "interaction_name"
		plot_y_col = x_axis_col
		x_order = row_order
		y_order = target_order
		x_display_labels = list(row_order)
		y_display_labels = _build_publication_pair_tick_labels(
			target_df=target_df,
			pair_order=target_order,
			wrap_width=y_label_wrap_width,
		)
		default_x_axis_label = "Ligand-Receptor Interaction or Pathway"
		default_y_axis_label = "Source Cell Type → Target Cell Type"
		x_rotation = 90
		resolved_y_tick_labelsize = _publication_y_tick_labelsize(
			n_labels=len(target_order),
			x_axis_mode=x_axis_mode,
		)
	else:
		plot_x_col = x_axis_col
		plot_y_col = "interaction_name"
		x_order = target_order
		y_order = row_order
		x_display_labels = _build_publication_x_tick_labels(
			target_df=target_df,
			target_order=target_order,
			x_axis_mode=x_axis_mode,
			x_label_wrap_width=x_label_wrap_width,
		)
		y_display_labels = [_wrap_label(str(label), 28) for label in row_order]
		default_x_axis_label = "Target Cell Type"
		default_y_axis_label = "Ligand-Receptor Interaction or Pathway"
		x_rotation = 0
		resolved_y_tick_labelsize = _publication_y_tick_labelsize(
			n_labels=len(row_order),
			x_axis_mode=x_axis_mode,
		)

	resolved_x_axis_label = x_axis_label or default_x_axis_label
	resolved_y_axis_label = y_axis_label or default_y_axis_label

	resolved_x_tick_labelsize = x_tick_labelsize
	if resolved_x_tick_labelsize is None:
		resolved_x_tick_labelsize = _publication_tick_labelsize(
			n_labels=len(x_order),
			x_axis_mode=x_axis_mode,
		)

	color_min = float(target_df["mean_rank_pct"].min())
	color_max = float(target_df["mean_rank_pct"].max())
	if abs(color_max - color_min) < 1e-12:
		color_max = color_min + 1e-6
	color_norm = mcolors.Normalize(vmin=color_min, vmax=color_max)

	size_min = float(target_df["top10_freq"].min())
	size_max = float(target_df["top10_freq"].max())
	if abs(size_max - size_min) < 1e-12:
		size_max = size_min + 1e-6
	size_norm = mcolors.Normalize(vmin=size_min, vmax=size_max)

	max_x_lines = max((len(label.splitlines()) for label in x_display_labels), default=1)
	max_x_chars = max((max(len(line) for line in label.splitlines()) for label in x_display_labels), default=12)
	max_y_chars = max((max(len(line) for line in label.splitlines()) for label in y_display_labels), default=18)

	if x_axis_mode == "source_target":
		left_margin = min(0.52, 0.16 + 0.0055 * max_y_chars)
	else:
		left_margin = min(0.35, 0.12 + 0.0045 * max_y_chars)
	bottom_margin = 0.18
	if x_axis_mode == "source_target":
		bottom_margin = 0.18 + 0.005 * max(0, max_x_chars - 10)
	bottom_margin += 0.04 * max(0, max_x_lines - 1)
	bottom_margin += 0.002 * max(0, max_x_chars - 16)
	bottom_margin = min(bottom_margin, 0.40)

	with sns.axes_style(
		"whitegrid",
		rc={
			"axes.spines.right": False,
			"axes.spines.top": False,
			"grid.color": "#d9d9d9",
			"grid.linewidth": 0.6,
			"grid.alpha": 0.30,
		},
	), sns.plotting_context("paper", font_scale=1.25):
		fig, ax = plt.subplots(figsize=figsize)
		fig.subplots_adjust(left=left_margin, right=0.79, bottom=bottom_margin, top=0.88)

		sns.scatterplot(
			data=target_df,
			x=plot_x_col,
			y=plot_y_col,
			size="top10_freq",
			hue="mean_rank_pct",
			palette=palette,
			hue_norm=color_norm,
			size_norm=size_norm,
			sizes=size_range,
			legend="brief",
			edgecolor="white",
			linewidth=0.6,
			alpha=0.95,
			ax=ax,
		)

		ax.set_axisbelow(True)
		ax.set_xlabel(resolved_x_axis_label, labelpad=12)
		ax.set_ylabel(resolved_y_axis_label, labelpad=12)
		ax.set_title(title, pad=16, weight="semibold")

		ax.set_xticks(range(len(x_order)))
		ax.set_xticklabels(x_display_labels)
		ax.tick_params(axis="x", labelsize=resolved_x_tick_labelsize, rotation=x_rotation, pad=6)
		for label in ax.get_xticklabels():
			label.set_horizontalalignment("center" if x_rotation == 0 else "right")
			label.set_verticalalignment("top")
			label.set_linespacing(1.1)

		ax.set_yticks(range(len(y_order)))
		ax.set_yticklabels(y_display_labels)
		ax.tick_params(
			axis="x",
			which="major",
			labelsize=resolved_x_tick_labelsize,
			rotation=x_rotation,
			pad=6,
			bottom=True,
			length=3.5,
			width=0.8,
			color="black",
			labelcolor="black",
		)
		ax.tick_params(
			axis="y",
			which="major",
			labelsize=resolved_y_tick_labelsize,
			left=True,
			length=3.5,
			width=0.8,
			color="black",
			labelcolor="black",
		)
		ax.grid(axis="x", visible=False)
		ax.grid(axis="y", color="#d9d9d9", linewidth=0.5, alpha=0.25)

		sns.despine(ax=ax, offset=5, trim=True)
		ax.spines["left"].set_color("black")
		ax.spines["bottom"].set_color("black")
		ax.spines["left"].set_linewidth(0.9)
		ax.spines["bottom"].set_linewidth(0.9)
		sns.move_legend(
			ax,
			"center left",
			bbox_to_anchor=(1.02, 0.5),
			frameon=False,
			title=None,
			borderaxespad=0.0,
		)
		_format_publication_percent_legend(
			ax,
			hue_label=hue_legend_label,
			size_label=size_legend_label,
			hue_label_wrap_width=hue_legend_wrap_width,
			size_label_wrap_width=size_legend_wrap_width,
			hue_title_fontsize=hue_legend_title_fontsize,
			size_title_fontsize=size_legend_title_fontsize,
			hue_value_fontsize=hue_legend_value_fontsize,
			size_value_fontsize=size_legend_value_fontsize,
		)

	return fig, ax, target_df
