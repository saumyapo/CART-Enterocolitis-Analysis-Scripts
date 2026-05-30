from __future__ import annotations

import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch


def plot_interaction_circle_ranks(
	df,
	cutoff_rank=0.995,
	source_col="source",
	target_col="target",
	interaction_col="interaction_name",
	rank_col="rank_median",
	min_n_interactions=1,
	node_radius=1.0,
	label_radius=1.18,
	label_wrap_width=14,
	figsize=(9, 9),
	font_size=9,
	title=None,
	title_x=0.5,
	title_y=0.975,
	subtitle_x=0.5,
	subtitle_y=0.945,
	legend_y=None,
	plot_top=None,
	cluster_nodes=True,
	cluster_metric="cosine",
):
	"""
	Circular directed interaction plot.

	Filters rows with rank_col >= cutoff_rank.
	Aggregates source-target pairs by number of unique interaction_name values.
	Filters source-target edges to retain only pairs with at least
	min_n_interactions unique interactions.
	Edge width and arrow size scale with interaction count.

	If cluster_nodes=True, nodes are ordered by similarity of their sent and
	received interaction profiles.
	"""

	required = [source_col, target_col, interaction_col, rank_col]
	missing = [c for c in required if c not in df.columns]
	if missing:
		raise ValueError(f"Missing required columns: {missing}")

	if min_n_interactions < 1:
		raise ValueError("min_n_interactions must be >= 1.")

	d = df[required].dropna().copy()
	d = d[d[rank_col] >= cutoff_rank].copy()

	if d.empty:
		raise ValueError("No rows left after applying cutoff_rank.")

	edges = (
		d.groupby([source_col, target_col], as_index=False)
		.agg(
			n_interactions=(interaction_col, "nunique"),
			interaction_names=(interaction_col, lambda x: ", ".join(sorted(set(x)))),
			rank_median=(rank_col, "median"),
		)
	)

	edges = edges[edges["n_interactions"] >= min_n_interactions].copy()

	if edges.empty:
		raise ValueError(
			"No source-target edges left after applying "
			f"cutoff_rank={cutoff_rank} and "
			f"min_n_interactions={min_n_interactions}."
		)

	nodes = sorted(set(edges[source_col]).union(edges[target_col]))

	if cluster_nodes:
		try:
			from scipy.cluster.hierarchy import leaves_list, linkage
			from scipy.spatial.distance import pdist

			send_mat = (
				edges.pivot_table(
					index=source_col,
					columns=target_col,
					values="n_interactions",
					aggfunc="sum",
					fill_value=0,
				)
				.reindex(index=nodes, columns=nodes, fill_value=0)
			)

			recv_mat = (
				edges.pivot_table(
					index=target_col,
					columns=source_col,
					values="n_interactions",
					aggfunc="sum",
					fill_value=0,
				)
				.reindex(index=nodes, columns=nodes, fill_value=0)
			)

			profile = pd.concat(
				[
					send_mat.add_prefix("send_to::"),
					recv_mat.add_prefix("recv_from::"),
				],
				axis=1,
			)

			row_sums = profile.sum(axis=1).replace(0, 1)
			profile = profile.div(row_sums, axis=0)

			if len(nodes) > 2 and profile.shape[1] > 1:
				dist = pdist(profile.values, metric=cluster_metric)

				if np.all(np.isfinite(dist)) and np.nanmax(dist) > 0:
					Z = linkage(dist, method="average")
					nodes = [nodes[i] for i in leaves_list(Z)]

		except Exception:
			nodes = sorted(set(edges[source_col]).union(edges[target_col]))

	angles = np.linspace(0, 2 * np.pi, len(nodes), endpoint=False)

	pos = {
		node: node_radius * np.array([np.cos(a), np.sin(a)])
		for node, a in zip(nodes, angles)
	}

	node_counts = (
		pd.concat(
			[
				edges[[source_col, "n_interactions"]].rename(columns={source_col: "node"}),
				edges[[target_col, "n_interactions"]].rename(columns={target_col: "node"}),
			]
		)
		.groupby("node")["n_interactions"]
		.sum()
		.reindex(nodes)
		.fillna(0)
	)

	outgoing = (
		edges.groupby(source_col)["n_interactions"]
		.sum()
		.reindex(nodes)
		.fillna(0)
	)

	incoming = (
		edges.groupby(target_col)["n_interactions"]
		.sum()
		.reindex(nodes)
		.fillna(0)
	)

	net_signal = outgoing - incoming
	max_abs = max(abs(net_signal).max(), 1e-9)

	def node_color(v):
		if v > 0:
			s = min(abs(v) / max_abs, 1.0)
			return (1.0, 0.85 - 0.5 * s, 0.85 - 0.5 * s)
		if v < 0:
			s = min(abs(v) / max_abs, 1.0)
			return (0.85 - 0.5 * s, 0.90 - 0.4 * s, 1.0)
		return (0.9, 0.9, 0.9)

	node_colors = [node_color(net_signal[n]) for n in nodes]

	if node_counts.max() == node_counts.min():
		node_sizes = np.repeat(900, len(nodes))
	else:
		node_sizes = 500 + 1600 * (
			node_counts - node_counts.min()
		) / (node_counts.max() - node_counts.min())

	counts = edges["n_interactions"].to_numpy(dtype=float)

	if counts.max() == counts.min():
		edge_widths = np.repeat(2.5, len(edges))
	else:
		edge_widths = 1.0 + 5.0 * (
			counts - counts.min()
		) / (counts.max() - counts.min())

	fig, ax = plt.subplots(figsize=figsize)

	ax.scatter(
		[pos[n][0] for n in nodes],
		[pos[n][1] for n in nodes],
		s=node_sizes,
		c=node_colors,
		edgecolors="black",
		linewidths=1.2,
		zorder=2,
	)

	pair_seen = {}

	for i, row in edges.reset_index(drop=True).iterrows():
		u = row[source_col]
		v = row[target_col]

		if u == v:
			continue

		reverse_exists = ((edges[source_col] == v) & (edges[target_col] == u)).any()
		key = tuple(sorted([u, v]))

		if reverse_exists:
			pair_seen[key] = pair_seen.get(key, 0) + 1
			rad = 0.18 if pair_seen[key] == 1 else -0.18
		else:
			rad = 0.10

		width = edge_widths[i]

		arrow = FancyArrowPatch(
			posA=pos[u],
			posB=pos[v],
			connectionstyle=f"arc3,rad={rad}",
			arrowstyle="-|>",
			mutation_scale=10 + 2 * width,
			linewidth=width,
			color="black",
			alpha=0.65,
			shrinkA=18,
			shrinkB=18,
			zorder=1,
		)
		ax.add_patch(arrow)

	for node in nodes:
		x, y = pos[node]
		r = np.sqrt(x**2 + y**2)
		xu, yu = x / r, y / r

		lx, ly = label_radius * xu, label_radius * yu
		ha = "left" if lx >= 0 else "right"
		va = "bottom" if ly > 0.15 else ("top" if ly < -0.15 else "center")

		label = "\n".join(
			textwrap.wrap(str(node), width=label_wrap_width, break_long_words=False)
		)

		ax.text(
			lx,
			ly,
			label,
			ha=ha,
			va=va,
			fontsize=font_size,
			fontweight="bold",
		)

	legend = [
		Line2D([0], [0], color="black", lw=1.5, alpha=0.65, label="Fewer interactions"),
		Line2D([0], [0], color="black", lw=5.5, alpha=0.65, label="More interactions"),
		Line2D(
			[0],
			[0],
			marker="o",
			color="w",
			markerfacecolor=(1, 0.4, 0.4),
			markeredgecolor="black",
			markersize=10,
			label="Net signaler",
		),
		Line2D(
			[0],
			[0],
			marker="o",
			color="w",
			markerfacecolor=(0.4, 0.6, 1.0),
			markeredgecolor="black",
			markersize=10,
			label="Net receiver",
		),
	]

	default_subtitle = (
		f"Cell-cell interactions with {rank_col} ≥ {cutoff_rank} "
		f"and ≥ {min_n_interactions} interaction(s) per edge"
	)
	resolved_legend_y = legend_y if legend_y is not None else (0.98 if title is None else 0.91)
	resolved_plot_top = plot_top if plot_top is not None else (0.93 if title is None else 0.82)

	fig.legend(
		handles=legend,
		frameon=False,
		loc="upper center",
		bbox_to_anchor=(0.5, resolved_legend_y),
		ncol=2,
		columnspacing=1.4,
		handlelength=1.6,
		handletextpad=0.6,
		fontsize=max(font_size - 1, 7),
	)

	if title is None:
		ax.set_title(default_subtitle, fontsize=font_size + 3, pad=12)
		layout_rect = (0.02, 0.02, 0.98, resolved_plot_top)
	else:
		fig.suptitle(title, fontsize=font_size + 4, x=title_x, ha="center", y=title_y)
		fig.text(
			subtitle_x,
			subtitle_y,
			default_subtitle,
			ha="center",
			va="top",
			fontsize=font_size + 1,
		)
		layout_rect = (0.02, 0.02, 0.98, resolved_plot_top)
	ax.set_aspect("equal")
	ax.axis("off")

	lim = label_radius + 0.25
	ax.set_xlim(-lim, lim)
	ax.set_ylim(-lim, lim)

	fig.tight_layout(rect=layout_rect)
	position = ax.get_position()
	if position.y1 < resolved_plot_top:
		ax.set_position(
			[position.x0, position.y0 + (resolved_plot_top - position.y1), position.width, position.height]
		)

	return fig, ax, edges