from __future__ import annotations

import textwrap

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch


def plot_celltype_circle_network_publication_updated(
	results,
	source_col="source",
	target_col="target",
	coef_col="coef_has_side_effect",
	qval_col="qval_has_side_effect",
	q_thresh=0.05,
	agg="mean",
	min_abs_coef=0.0,
	min_interactions_per_edge=1,
	top_n_edges=None,
	balance_directions=True,
	node_size_mode="strength",
	node_color_mode="role",
	role_balance_tol=0.15,
	label_outside=True,
	label_wrap_width=12,
	label_radius=1.08,
	node_radius=0.92,
	font_size=9,
	figsize=(10, 10),
	title="Cell-cell interaction rewiring",
	cluster_nodes=True,
	cluster_metric="cosine",
):
	"""
	Circular directed network for aggregated cell-type rewiring.

	Edges
	-----
	source -> target

	Edge color
	----------
	Red  = positive aggregated effect
	Blue = negative aggregated effect

	Edge width
	----------
	Absolute aggregated effect magnitude

	Edge alpha
	----------
	Confidence from q-value

	Node size
	---------
	Communication burden, either total strength or degree.

	Node color
	----------
	If node_color_mode="role":
		Orange = sender-enriched among plotted edges
		Purple = receiver-enriched among plotted edges
		Gray = balanced/mixed among plotted edges
	"""

	def _wrap_label(text, width):
		return "\n".join(textwrap.wrap(str(text), width=width, break_long_words=False))

	def _safe_scale(values, low, high):
		values = np.asarray(values, dtype=float)
		if len(values) == 0:
			return values
		vmin, vmax = np.nanmin(values), np.nanmax(values)
		if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax == vmin:
			return np.repeat((low + high) / 2, len(values))
		return low + (high - low) * (values - vmin) / (vmax - vmin)

	required = [source_col, target_col, coef_col, qval_col]
	missing = [col for col in required if col not in results.columns]
	if missing:
		raise ValueError(f"Missing required columns: {missing}")

	data = results[required].dropna().copy()
	data[coef_col] = pd.to_numeric(data[coef_col], errors="coerce")
	data[qval_col] = pd.to_numeric(data[qval_col], errors="coerce")

	data = data[np.isfinite(data[coef_col]) & np.isfinite(data[qval_col])].copy()
	data = data[data[qval_col] <= q_thresh].copy()

	if min_abs_coef > 0:
		data = data[data[coef_col].abs() >= min_abs_coef].copy()

	if data.empty:
		raise ValueError("No rows left after filtering.")

	grouped = (
		data.groupby([source_col, target_col], as_index=False)
		.agg(
			coef_mean=(coef_col, "mean"),
			coef_median=(coef_col, "median"),
			coef_sum=(coef_col, "sum"),
			n_interactions=(coef_col, "size"),
			min_q=(qval_col, "min"),
		)
	)

	if agg == "mean":
		grouped["weight"] = grouped["coef_mean"]
	elif agg == "median":
		grouped["weight"] = grouped["coef_median"]
	elif agg == "sum":
		grouped["weight"] = grouped["coef_sum"]
	else:
		raise ValueError("agg must be one of: 'mean', 'median', 'sum'")

	grouped = grouped[grouped["n_interactions"] >= min_interactions_per_edge].copy()

	if grouped.empty:
		raise ValueError("No source-target pairs left after aggregation/filtering.")

	grouped["abs_weight"] = grouped["weight"].abs()

	if top_n_edges is not None and len(grouped) > top_n_edges:
		if balance_directions:
			pos = grouped[grouped["weight"] > 0].sort_values("abs_weight", ascending=False)
			neg = grouped[grouped["weight"] < 0].sort_values("abs_weight", ascending=False)

			n_pos = top_n_edges // 2
			n_neg = top_n_edges - n_pos

			keep = pd.concat([pos.head(n_pos), neg.head(n_neg)])

			if len(keep) < top_n_edges:
				extra = (
					grouped.loc[~grouped.index.isin(keep.index)]
					.sort_values("abs_weight", ascending=False)
					.head(top_n_edges - len(keep))
				)
				keep = pd.concat([keep, extra])

			grouped = keep.copy()
		else:
			grouped = grouped.sort_values("abs_weight", ascending=False).head(top_n_edges).copy()

	grouped["neglog10_q"] = -np.log10(grouped["min_q"].clip(lower=1e-300))

	graph = nx.DiGraph()
	nodes = sorted(set(grouped[source_col]).union(set(grouped[target_col])))
	graph.add_nodes_from(nodes)

	for _, row in grouped.iterrows():
		graph.add_edge(
			row[source_col],
			row[target_col],
			weight=row["weight"],
			abs_weight=row["abs_weight"],
			n_interactions=row["n_interactions"],
			min_q=row["min_q"],
			neglog10_q=row["neglog10_q"],
		)

	if graph.number_of_edges() == 0:
		raise ValueError("No edges to plot after filtering.")

	node_rows = []
	for node in graph.nodes():
		outgoing = [graph[node][neighbor]["weight"] for neighbor in graph.successors(node)]
		incoming = [graph[neighbor][node]["weight"] for neighbor in graph.predecessors(node)]

		outgoing_abs = np.array([abs(value) for value in outgoing], dtype=float)
		incoming_abs = np.array([abs(value) for value in incoming], dtype=float)

		out_strength = float(outgoing_abs.sum()) if len(outgoing_abs) else 0.0
		in_strength = float(incoming_abs.sum()) if len(incoming_abs) else 0.0
		total_strength = out_strength + in_strength

		sender_fraction = out_strength / total_strength if total_strength > 0 else 0.5
		role_score = (out_strength - in_strength) / total_strength if total_strength > 0 else 0.0

		if abs(role_score) <= role_balance_tol:
			communication_role = "Balanced"
		elif role_score > 0:
			communication_role = "Sender-enriched"
		else:
			communication_role = "Receiver-enriched"

		node_rows.append(
			{
				"celltype": node,
				"out_degree": graph.out_degree(node),
				"in_degree": graph.in_degree(node),
				"total_degree": graph.degree(node),
				"out_strength": out_strength,
				"in_strength": in_strength,
				"total_strength": total_strength,
				"sender_fraction": sender_fraction,
				"role_score": role_score,
				"communication_role": communication_role,
				"signed_net_effect": float(np.sum(outgoing) - np.sum(incoming)) if outgoing or incoming else 0.0,
			}
		)

	node_df = (
		pd.DataFrame(node_rows)
		.sort_values(["total_strength", "out_strength", "in_strength"], ascending=False)
		.reset_index(drop=True)
	)

	if cluster_nodes:
		try:
			from scipy.cluster.hierarchy import leaves_list, linkage
			from scipy.spatial.distance import pdist

			all_nodes = node_df["celltype"].tolist()

			send_mat = (
				grouped.pivot_table(
					index=source_col,
					columns=target_col,
					values="abs_weight",
					aggfunc="sum",
					fill_value=0,
				)
				.reindex(index=all_nodes, columns=all_nodes, fill_value=0)
			)

			recv_mat = (
				grouped.pivot_table(
					index=target_col,
					columns=source_col,
					values="abs_weight",
					aggfunc="sum",
					fill_value=0,
				)
				.reindex(index=all_nodes, columns=all_nodes, fill_value=0)
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

			if len(all_nodes) > 2 and profile.shape[1] > 1:
				dist = pdist(profile.values, metric=cluster_metric)

				if np.all(np.isfinite(dist)) and np.nanmax(dist) > 0:
					ordered_nodes = [all_nodes[i] for i in leaves_list(linkage(dist, method="average"))]
				else:
					ordered_nodes = node_df["celltype"].tolist()
			else:
				ordered_nodes = node_df["celltype"].tolist()

		except Exception:
			ordered_nodes = node_df["celltype"].tolist()
	else:
		ordered_nodes = node_df["celltype"].tolist()

	angles = np.linspace(0, 2 * np.pi, len(ordered_nodes), endpoint=False)
	pos = {
		node: node_radius * np.array([np.cos(angle), np.sin(angle)])
		for node, angle in zip(ordered_nodes, angles)
	}

	base_lookup = node_df.set_index("celltype")

	if node_size_mode == "strength":
		base = base_lookup["total_strength"].reindex(ordered_nodes)
	elif node_size_mode == "degree":
		base = base_lookup["total_degree"].reindex(ordered_nodes)
	else:
		raise ValueError("node_size_mode must be 'strength' or 'degree'")

	node_size_values = _safe_scale(base.values, 500, 2100)
	node_sizes = dict(zip(ordered_nodes, node_size_values))

	role_colors = {
		"Sender-enriched": "#E69F00",
		"Receiver-enriched": "#7B6FD6",
		"Balanced": "#D9D9D9",
	}

	if node_color_mode == "role":
		node_colors = [role_colors[base_lookup.loc[node, "communication_role"]] for node in ordered_nodes]
	elif node_color_mode == "burden":
		burden = base_lookup["total_strength"].reindex(ordered_nodes).values
		burden_scaled = _safe_scale(burden, 0.15, 0.85)
		node_colors = [
			(1.0 - 0.35 * value, 1.0 - 0.35 * value, 1.0 - 0.35 * value)
			for value in burden_scaled
		]
	elif node_color_mode == "none":
		node_colors = ["#EDEDED" for _ in ordered_nodes]
	else:
		raise ValueError("node_color_mode must be one of: 'role', 'burden', 'none'")

	edge_list = list(graph.edges())
	edge_weights = np.array([graph[u][v]["weight"] for u, v in edge_list], dtype=float)
	abs_w = np.abs(edge_weights)
	q_strength = np.array([graph[u][v]["neglog10_q"] for u, v in edge_list], dtype=float)

	edge_widths = _safe_scale(abs_w, 0.9, 5.5)
	edge_alphas = _safe_scale(q_strength, 0.25, 0.95)
	edge_colors = ["#C93C37" if weight > 0 else "#2B6CB0" for weight in edge_weights]

	fig, ax = plt.subplots(figsize=figsize)

	nx.draw_networkx_nodes(
		graph,
		pos,
		nodelist=ordered_nodes,
		node_size=[node_sizes[node] for node in ordered_nodes],
		node_color=node_colors,
		edgecolors="black",
		linewidths=1.0,
		ax=ax,
	)

	drawn_pairs = set()

	for index, (u, v) in enumerate(edge_list):
		same_pair_reverse = (v, u) in graph.edges()

		if same_pair_reverse:
			rad = 0.16 if (u, v) not in drawn_pairs else -0.16
		else:
			rad = 0.08

		drawn_pairs.add((u, v))

		patch = FancyArrowPatch(
			posA=pos[u],
			posB=pos[v],
			connectionstyle=f"arc3,rad={rad}",
			arrowstyle="-|>",
			mutation_scale=10 + edge_widths[index],
			linewidth=edge_widths[index],
			color=edge_colors[index],
			alpha=edge_alphas[index],
			shrinkA=14,
			shrinkB=14,
			zorder=1,
		)
		ax.add_patch(patch)

	if label_outside:
		for node in ordered_nodes:
			x, y = pos[node]
			r = np.sqrt(x**2 + y**2)
			x_unit, y_unit = x / r, y / r

			x2 = label_radius * x_unit
			y2 = label_radius * y_unit

			ha = "left" if x2 >= 0 else "right"
			va = "bottom" if y2 > 0.15 else ("top" if y2 < -0.15 else "center")

			ax.text(
				x2,
				y2,
				_wrap_label(node, label_wrap_width),
				fontsize=font_size,
				fontweight="bold",
				ha=ha,
				va=va,
				bbox=dict(
					boxstyle="round,pad=0.15",
					fc="white",
					ec="none",
					alpha=0.85,
				),
				zorder=3,
			)
	else:
		wrapped_labels = {node: _wrap_label(node, label_wrap_width) for node in ordered_nodes}
		nx.draw_networkx_labels(
			graph,
			pos,
			labels=wrapped_labels,
			font_size=font_size,
			font_weight="bold",
			ax=ax,
		)

	legend_elements = [
		Line2D([0], [0], color="#C93C37", lw=3, label="Positive effect"),
		Line2D([0], [0], color="#2B6CB0", lw=3, label="Negative effect"),
		Line2D([0], [0], color="black", lw=1.5, alpha=0.30, label="Lower confidence"),
		Line2D([0], [0], color="black", lw=1.5, alpha=0.90, label="Higher confidence"),
	]

	if node_color_mode == "role":
		legend_elements.extend(
			[
				Line2D([0], [0], marker="o", color="w", markerfacecolor=role_colors["Sender-enriched"], markeredgecolor="black", markersize=10, label="Sender-enriched node"),
				Line2D([0], [0], marker="o", color="w", markerfacecolor=role_colors["Receiver-enriched"], markeredgecolor="black", markersize=10, label="Receiver-enriched node"),
				Line2D([0], [0], marker="o", color="w", markerfacecolor=role_colors["Balanced"], markeredgecolor="black", markersize=10, label="Balanced node"),
			]
		)
	elif node_color_mode == "burden":
		legend_elements.append(
			Line2D([0], [0], marker="o", color="w", markerfacecolor="#999999", markeredgecolor="black", markersize=10, label="Darker = higher burden")
		)

	ax.legend(
		handles=legend_elements,
		frameon=False,
		loc="upper right",
		bbox_to_anchor=(1.12, 1.02),
		fontsize=max(font_size - 1, 7),
	)

	ax.text(
		0,
		-(label_radius + 0.28),
		(
			"Edges are filtered before plotting. "
			"Edge width = effect magnitude | Edge opacity = confidence | "
			"Node size = plotted communication burden | "
			"Node color = sender/receiver enrichment among plotted edges"
		),
		ha="center",
		va="center",
		fontsize=max(font_size - 1, 7),
	)

	lim = label_radius + 0.18
	ax.set_xlim(-lim, lim)
	ax.set_ylim(-lim, lim)
	ax.set_aspect("equal")
	ax.set_title(title, fontsize=font_size + 3, pad=10)
	ax.axis("off")
	plt.tight_layout()

	grouped = grouped.sort_values("abs_weight", ascending=False).reset_index(drop=True)
	node_df = node_df.sort_values("total_strength", ascending=False).reset_index(drop=True)

	return fig, ax, grouped, node_df
