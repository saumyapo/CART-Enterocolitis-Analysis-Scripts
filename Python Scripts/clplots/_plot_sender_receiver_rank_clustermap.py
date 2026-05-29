
from __future__ import annotations

import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot_sender_receiver_rank_clustermap(
	results,
	rank_col="rank_mean",
	source_col="source",
	target_col="target",
	feature_col="pathway_name",   # "pathway_name" or "interaction_name"
	rank_thresh=0.99,
	agg="median",
	top_n_pairs=50,
	top_n_features=30,
	figsize=(12, 14),
	cmap="viridis",
	cluster_metric="euclidean",
	cluster_method="average",
	z_score=None,
	standard_scale=None,
	title=None,
	min_unique_features_per_pair=1,

	# Layout
	dendrogram_size_inches=0.30,
	dendrogram_gap_inches=0.015,
	yticklabel_space_inches=3.5,
	bottom_space_inches=1.4,
	top_space_inches=0.65,
	cbar_height_inches=0.10,
	cbar_gap_inches=0.08,

	# Labels
	row_label_wrap=None,
	col_label_wrap=14,
	yticklabel_fontsize=8,
	xticklabel_fontsize=8,

	# Cell annotations
	annotate_counts=False,
	count_fontsize=5,
	count_fmt=".0f",
	count_color="black",
):
	"""
	Clustered heatmap of source→target communication ranks.

	Rows:
		source → target pairs

	Columns:
		feature_col, usually "pathway_name" or "interaction_name"

	Cell color:
		aggregated rank_col value, usually median rank_mean

	Filtering:
		keeps rows where rank_col > rank_thresh

	Returns:
		g, mat, grouped
	"""

	# ---------------------------
	# Validate
	# ---------------------------
	required = [rank_col, source_col, target_col, feature_col]
	missing = [c for c in required if c not in results.columns]
	if missing:
		raise ValueError(f"Missing required columns: {missing}")

	if agg not in {"mean", "median", "max", "min"}:
		raise ValueError("agg must be one of: 'mean', 'median', 'max', 'min'")

	if z_score not in {None, "row", "col"}:
		raise ValueError("z_score must be one of: None, 'row', 'col'")

	if standard_scale not in {None, "row", "col"}:
		raise ValueError("standard_scale must be one of: None, 'row', 'col'")

	if z_score is not None and standard_scale is not None:
		raise ValueError("Use either z_score or standard_scale, not both.")

	# ---------------------------
	# Filter
	# ---------------------------
	d = results[required].dropna().copy()
	d[rank_col] = pd.to_numeric(d[rank_col], errors="coerce")

	d = d[
		np.isfinite(d[rank_col])
		& (d[rank_col] > rank_thresh)
	].copy()

	if d.empty:
		raise ValueError(
			f"No rows remain after filtering {rank_col} > {rank_thresh}."
		)

	# ---------------------------
	# Aggregate
	# ---------------------------
	grouped = (
		d.groupby([source_col, target_col, feature_col], observed=True)
		.agg(
			rank_mean=(rank_col, "mean"),
			rank_median=(rank_col, "median"),
			rank_max=(rank_col, "max"),
			rank_min=(rank_col, "min"),
			n_entries=(rank_col, "size"),
		)
		.reset_index()
	)

	grouped["rank_value"] = grouped[
		{
			"mean": "rank_mean",
			"median": "rank_median",
			"max": "rank_max",
			"min": "rank_min",
		}[agg]
	]

	grouped["pair"] = (
		grouped[source_col].astype(str)
		+ " → "
		+ grouped[target_col].astype(str)
	)

	# ---------------------------
	# Filter source→target pairs by feature diversity
	# ---------------------------
	if min_unique_features_per_pair is not None:
		if min_unique_features_per_pair < 1:
			raise ValueError("min_unique_features_per_pair must be >= 1 or None.")

		pair_feature_counts = grouped.groupby("pair")[feature_col].nunique()

		keep_pairs = pair_feature_counts[
			pair_feature_counts >= min_unique_features_per_pair
		].index

		grouped = grouped[grouped["pair"].isin(keep_pairs)].copy()

		if grouped.empty:
			raise ValueError(
				"No rows remain after filtering source→target pairs with fewer than "
				f"{min_unique_features_per_pair} unique {feature_col} values."
			)

	# ---------------------------
	# Select top rows/columns
	# ---------------------------
	top_pairs = (
		grouped.groupby("pair")["rank_value"]
		.median()
		.sort_values(ascending=False)
		.head(top_n_pairs)
		.index
	)

	top_features = (
		grouped.groupby(feature_col)["rank_value"]
		.median()
		.sort_values(ascending=False)
		.head(top_n_features)
		.index
	)

	grouped = grouped[
		grouped["pair"].isin(top_pairs)
		& grouped[feature_col].isin(top_features)
	].copy()

	if grouped.empty:
		raise ValueError("No rows remain after top row/column selection.")

	# ---------------------------
	# Matrices
	# ---------------------------
	mat = grouped.pivot_table(
		index="pair",
		columns=feature_col,
		values="rank_value",
		aggfunc="mean",
		fill_value=0,
	)

	count_mat = grouped.pivot_table(
		index="pair",
		columns=feature_col,
		values="n_entries",
		aggfunc="sum",
		fill_value=0,
	)

	count_mat = count_mat.reindex(index=mat.index, columns=mat.columns).fillna(0)

	if mat.empty:
		raise ValueError("The plotted matrix is empty.")

	mat = mat.astype(float)
	mat = mat.replace([np.inf, -np.inf], np.nan).fillna(0)

	# Drop completely empty rows/columns before clustering
	mat = mat.loc[mat.sum(axis=1) > 0, mat.sum(axis=0) > 0]
	count_mat = count_mat.reindex(index=mat.index, columns=mat.columns).fillna(0)

	if mat.empty:
		raise ValueError("The plotted matrix is empty after dropping zero rows/columns.")

	# Mask missing/zero cells visually, but keep finite zeroes for clustering
	mask = mat.eq(0)

	# ---------------------------
	# Label wrapping
	# ---------------------------
	raw_row_labels = mat.index.astype(str)
	raw_col_labels = mat.columns.astype(str)

	row_labels = [
		"\n".join(textwrap.wrap(x, width=row_label_wrap))
		if row_label_wrap is not None else x
		for x in raw_row_labels
	]

	col_labels = [
		"\n".join(textwrap.wrap(x, width=col_label_wrap))
		if col_label_wrap is not None else x
		for x in raw_col_labels
	]

	mat.index = row_labels
	mat.columns = col_labels
	count_mat.index = row_labels
	count_mat.columns = col_labels
	mask.index = row_labels
	mask.columns = col_labels

	# ---------------------------
	# Color scale
	# ---------------------------
	positive_values = mat.values[mat.values > 0]
	if positive_values.size == 0:
		raise ValueError("No positive rank values remain for plotting.")

	vmin = max(rank_thresh, float(np.nanmin(positive_values)))
	vmax = float(np.nanmax(positive_values))

	if vmax <= vmin:
		vmax = vmin + 1e-9

	# ---------------------------
	# Layout ratios
	# ---------------------------
	fig_w, fig_h = figsize

	row_dendro_ratio = min(dendrogram_size_inches / fig_w, 0.20)
	col_dendro_ratio = min(dendrogram_size_inches / fig_h, 0.20)

	right_margin = min(yticklabel_space_inches / fig_w, 0.55)
	bottom_margin = min(bottom_space_inches / fig_h, 0.35)
	top_margin = min(top_space_inches / fig_h, 0.20)

	z_score_axis = 0 if z_score == "row" else 1 if z_score == "col" else None
	standard_scale_axis = (
		0 if standard_scale == "row"
		else 1 if standard_scale == "col"
		else None
	)

	row_cluster = mat.shape[0] > 1
	col_cluster = mat.shape[1] > 1

	annot_data = count_mat if annotate_counts else False

	if title is None:
		title = f"Source → target {feature_col} rank"

	# ---------------------------
	# Plot
	# ---------------------------
	g = sns.clustermap(
		mat,
		mask=mask,
		cmap=cmap,
		vmin=vmin,
		vmax=vmax,
		figsize=figsize,
		xticklabels=True,
		yticklabels=True,
		row_cluster=row_cluster,
		col_cluster=col_cluster,
		method=cluster_method,
		metric=cluster_metric,
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
			"label": f"{agg.capitalize()} {rank_col}",
		},
	)

	# ---------------------------
	# Layout repair
	# ---------------------------
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

	# ---------------------------
	# Styling
	# ---------------------------
	g.ax_heatmap.yaxis.tick_right()
	g.ax_heatmap.yaxis.set_label_position("right")

	g.ax_heatmap.tick_params(
		axis="y",
		labelsize=yticklabel_fontsize,
		length=0,
		pad=6,
	)

	g.ax_heatmap.tick_params(
		axis="x",
		labelsize=xticklabel_fontsize,
		length=0,
		pad=4,
	)

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

	g.fig.suptitle(
		title,
		fontsize=14,
		fontweight="bold",
		y=1 - (0.15 / fig_h),
	)

	return g, mat, grouped