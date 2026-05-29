
from __future__ import annotations

import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot_sender_receiver_pathway_clustermap(
	results,
	coef_col="coef_is_treated",
	qval_col="qval_is_treated",
	source_col="source",
	target_col="target",
	pathway_col="pathway_name",
	q_thresh=0.05,
	effect_thresh=0.0,
	agg="mean",
	top_n_pairs=50,
	top_n_pathways=30,
	figsize=(12, 14),
	cmap="coolwarm",
	z_score=None,
	standard_scale=None,
	title="Source → target pathway rewiring (clustered)",
	vmin=None,
	vmax=None,
	center=None,
	dendrogram_size_inches=0.30,
	dendrogram_gap_inches=0.015,
	yticklabel_space_inches=3.5,
	bottom_space_inches=1.4,
	top_space_inches=0.65,
	cbar_height_inches=0.10,
	cbar_gap_inches=0.08,
	row_label_wrap=None,
	col_label_wrap=14,
	yticklabel_fontsize=8,
	xticklabel_fontsize=8,
	annotate_counts=False,
	count_fontsize=5,
	count_fmt=".0f",
	count_color="black",
):
	"""
	Clustered heatmap of source→target pathway effects.

	Cell color:
		aggregated signed coefficient

	Optional cell label:
		number of interactions aggregated into that cell
	"""

	required = [coef_col, qval_col, source_col, target_col, pathway_col]
	missing = [col for col in required if col not in results.columns]
	if missing:
		raise ValueError(f"Missing required columns: {missing}")

	if agg not in {"mean", "median", "sum"}:
		raise ValueError("agg must be one of: 'mean', 'median', 'sum'")

	if z_score not in {None, "row", "col"}:
		raise ValueError("z_score must be one of: None, 'row', 'col'")

	if standard_scale not in {None, "row", "col"}:
		raise ValueError("standard_scale must be one of: None, 'row', 'col'")

	if z_score is not None and standard_scale is not None:
		raise ValueError("Use either z_score or standard_scale, not both.")

	data = results[required].dropna().copy()
	data[coef_col] = pd.to_numeric(data[coef_col], errors="coerce")
	data[qval_col] = pd.to_numeric(data[qval_col], errors="coerce")

	data = data[
		np.isfinite(data[coef_col])
		& np.isfinite(data[qval_col])
		& (data[qval_col] <= q_thresh)
		& (data[coef_col].abs() >= effect_thresh)
	].copy()

	if data.empty:
		raise ValueError("No rows remain after filtering.")

	grouped = (
		data.groupby([source_col, target_col, pathway_col], observed=True)
		.agg(
			coef_mean=(coef_col, "mean"),
			coef_median=(coef_col, "median"),
			coef_sum=(coef_col, "sum"),
			n_interactions=(coef_col, "size"),
		)
		.reset_index()
	)

	grouped["effect"] = grouped[
		{
			"mean": "coef_mean",
			"median": "coef_median",
			"sum": "coef_sum",
		}[agg]
	]

	grouped["abs_effect"] = grouped["effect"].abs()
	grouped["pair"] = grouped[source_col].astype(str) + " → " + grouped[target_col].astype(str)

	top_pairs = (
		grouped.groupby("pair")["abs_effect"]
		.sum()
		.sort_values(ascending=False)
		.head(top_n_pairs)
		.index
	)

	top_pathways = (
		grouped.groupby(pathway_col)["abs_effect"]
		.sum()
		.sort_values(ascending=False)
		.head(top_n_pathways)
		.index
	)

	grouped = grouped[
		grouped["pair"].isin(top_pairs)
		& grouped[pathway_col].isin(top_pathways)
	].copy()

	if grouped.empty:
		raise ValueError("No rows remain after top row/column selection.")

	mat = grouped.pivot_table(
		index="pair",
		columns=pathway_col,
		values="effect",
		aggfunc="mean",
		fill_value=0,
	)

	count_mat = grouped.pivot_table(
		index="pair",
		columns=pathway_col,
		values="n_interactions",
		aggfunc="sum",
		fill_value=0,
	)

	count_mat = count_mat.reindex(index=mat.index, columns=mat.columns).fillna(0)

	if mat.empty:
		raise ValueError("The plotted matrix is empty.")

	raw_row_labels = mat.index.astype(str)
	raw_col_labels = mat.columns.astype(str)

	row_labels = [
		"\n".join(textwrap.wrap(value, width=row_label_wrap)) if row_label_wrap is not None else value
		for value in raw_row_labels
	]
	col_labels = [
		"\n".join(textwrap.wrap(value, width=col_label_wrap)) if col_label_wrap is not None else value
		for value in raw_col_labels
	]

	mat.index = row_labels
	mat.columns = col_labels
	count_mat.index = row_labels
	count_mat.columns = col_labels

	if vmin is None and vmax is None:
		vmax = np.nanmax(np.abs(mat.values))
		vmax = max(float(vmax), 1e-9)
		vmin = -vmax
	elif vmin is None:
		vmin = -abs(vmax)
	elif vmax is None:
		vmax = abs(vmin)

	if center is None:
		center = 0

	fig_w, fig_h = figsize
	row_dendro_ratio = min(dendrogram_size_inches / fig_w, 0.20)
	col_dendro_ratio = min(dendrogram_size_inches / fig_h, 0.20)
	right_margin = min(yticklabel_space_inches / fig_w, 0.55)
	bottom_margin = min(bottom_space_inches / fig_h, 0.35)
	top_margin = min(top_space_inches / fig_h, 0.20)

	z_score_axis = 0 if z_score == "row" else 1 if z_score == "col" else None
	standard_scale_axis = 0 if standard_scale == "row" else 1 if standard_scale == "col" else None
	row_cluster = mat.shape[0] > 1
	col_cluster = mat.shape[1] > 1
	annot_data = count_mat if annotate_counts else False

	g = sns.clustermap(
		mat,
		cmap=cmap,
		center=center,
		vmin=-vmax,
		vmax=vmax,
		figsize=figsize,
		xticklabels=True,
		yticklabels=True,
		row_cluster=row_cluster,
		col_cluster=col_cluster,
		method="average",
		metric="cosine",
		z_score=z_score_axis,
		standard_scale=standard_scale_axis,
		dendrogram_ratio=(row_dendro_ratio, col_dendro_ratio),
		linewidths=0.2,
		linecolor="white",
		annot=annot_data,
		fmt=count_fmt,
		annot_kws={
			"fontsize": count_fontsize,
			"color": count_color,
		},
		cbar_pos=(0.2, 0.02, 0.6, 0.02),
		cbar_kws={
			"orientation": "horizontal",
			"label": "Aggregated coefficient\n(red = increased, blue = decreased)",
		},
	)

	g.fig.subplots_adjust(
		right=1 - right_margin,
		bottom=bottom_margin,
		top=1 - top_margin,
	)

	g.fig.canvas.draw()

	heatmap_pos = g.ax_heatmap.get_position()
	dendro_width = dendrogram_size_inches / fig_w
	dendro_height = dendrogram_size_inches / fig_h
	dendro_gap_x = dendrogram_gap_inches / fig_w
	dendro_gap_y = dendrogram_gap_inches / fig_h

	if row_cluster:
		g.ax_row_dendrogram.set_position([
			heatmap_pos.x0 - dendro_width - dendro_gap_x,
			heatmap_pos.y0,
			dendro_width,
			heatmap_pos.height,
		])

	if col_cluster:
		g.ax_col_dendrogram.set_position([
			heatmap_pos.x0,
			heatmap_pos.y1 + dendro_gap_y,
			heatmap_pos.width,
			dendro_height,
		])

	cbar_height = cbar_height_inches / fig_h
	cbar_gap = cbar_gap_inches / fig_h

	g.cax.set_position([
		heatmap_pos.x0,
		heatmap_pos.y0 - cbar_gap - cbar_height,
		heatmap_pos.width,
		cbar_height,
	])

	g.ax_heatmap.yaxis.tick_right()
	g.ax_heatmap.yaxis.set_label_position("right")
	g.ax_heatmap.tick_params(axis="y", labelsize=yticklabel_fontsize, length=0, pad=6)
	g.ax_heatmap.tick_params(axis="x", labelsize=xticklabel_fontsize, length=0, pad=4)

	for label in g.ax_heatmap.get_yticklabels():
		label.set_rotation(0)
		label.set_horizontalalignment("left")

	for label in g.ax_heatmap.get_xticklabels():
		label.set_rotation(90)
		label.set_horizontalalignment("center")
		label.set_verticalalignment("top")

	g.cax.xaxis.set_label_position("bottom")
	g.cax.xaxis.tick_bottom()
	g.cax.tick_params(axis="x", labelsize=8, length=2, pad=2)

	g.fig.suptitle(title, fontsize=14, fontweight="bold", y=1 - (0.15 / fig_h))

	return g, mat, grouped
