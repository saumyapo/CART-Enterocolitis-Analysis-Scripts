import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text


def plot_side_effect_volcano_publication(
    results,
    coef_col="coef_has_side_effect",
    qval_col="qval_has_side_effect",
    label_col="interaction_name",
    source_col="source",
    target_col="target",
    q_thresh=0.05,
    effect_thresh=0.5,
    top_n=12,
    label_mode="full",
    label_interactions=None,
    override_auto_labels=False,
    figsize=(8, 6),
    title="Side-effect interaction volcano plot",
    xlabel="Colitis effect size",
    dpi=None,
    point_size_ns=16,
    point_size_sig=34,
    alpha_ns=0.28,
    alpha_sig=0.92,
    max_labels_per_side=None,
    min_label_separation_x=None,
    ymax=None,
    y_break=None,
    break_gap=5.0,
    upper_height_ratio=0.8,
    lower_height_ratio=4.2,
    hspace=0.10,
    show_break_label=True,
    legend_position="default",
    legend_loc=None,
    legend_bbox_to_anchor=None,
    legend_borderaxespad=None,
):
    """
    Volcano plot with optional broken y-axis.

    Set y_break=None to disable the y-axis break.
    Set y_break=<number> to use a broken y-axis only when max -log10(q) exceeds it.

    Legend placement:
        legend_position="default" keeps the current upper-right placement.
        legend_position="top-left" snaps the legend to the upper-left corner.
        legend_position="manual" uses legend_loc and legend_bbox_to_anchor.

    Label selection:
        label_interactions can be an iterable of (source, target, interaction_name)
        tuples to force specific labels onto the plot.
        override_auto_labels=True disables the automatic top_n label selection.

    Returns
    -------
    If y_break is used:
        fig, (ax, ax_top), plotted_dataframe
    Else:
        fig, ax, plotted_dataframe
    """

    required = [coef_col, qval_col, label_col, source_col, target_col]
    missing = [col for col in required if col not in results.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    data = results[required].dropna().copy()
    data[coef_col] = pd.to_numeric(data[coef_col], errors="coerce")
    data[qval_col] = pd.to_numeric(data[qval_col], errors="coerce")

    data = data[np.isfinite(data[coef_col]) & np.isfinite(data[qval_col])]
    data = data[data[qval_col] > 0].copy()

    if data.empty:
        raise ValueError("No valid rows remain after filtering.")

    data["neglog10_q"] = -np.log10(data[qval_col].clip(lower=1e-300))
    data["is_sig"] = (
        (data[qval_col] < q_thresh)
        & (np.abs(data[coef_col]) >= effect_thresh)
    )

    data["direction"] = "Not significant"
    data.loc[data["is_sig"] & (data[coef_col] > 0), "direction"] = "Increased"
    data.loc[data["is_sig"] & (data[coef_col] < 0), "direction"] = "Decreased"

    if label_mode == "interaction":
        data["plot_label"] = data[label_col].astype(str)
    elif label_mode == "cellpair":
        data["plot_label"] = (
            data[source_col].astype(str) + " → " + data[target_col].astype(str)
        )
    elif label_mode == "full":
        data["plot_label"] = (
            data[source_col].astype(str)
            + " → "
            + data[target_col].astype(str)
            + "\n"
            + data[label_col].astype(str)
        )
    else:
        raise ValueError("label_mode must be one of: 'interaction', 'cellpair', 'full'")

    if legend_position not in {"default", "top-left", "manual"}:
        raise ValueError(
            "legend_position must be one of: 'default', 'top-left', 'manual'"
        )

    if label_interactions is None:
        label_interactions = set()
    else:
        try:
            label_interactions = {
                tuple(str(part) for part in interaction)
                for interaction in label_interactions
            }
        except TypeError as exc:
            raise ValueError(
                "label_interactions must be an iterable of "
                "(source, target, interaction_name) tuples"
            ) from exc

        if any(len(interaction) != 3 for interaction in label_interactions):
            raise ValueError(
                "label_interactions must contain 3-item tuples of "
                "(source, target, interaction_name)"
            )

    data["label_score"] = data["neglog10_q"]

    if max_labels_per_side is None:
        max_labels_per_side = top_n

    coef_abs_max = float(np.nanmax(np.abs(data[coef_col].values)))
    x_lim = max(effect_thresh * 1.25, coef_abs_max * 1.08)

    y_max = float(np.nanmax(data["neglog10_q"].values))
    y_thresh_line = -np.log10(q_thresh)

    if min_label_separation_x is None:
        min_label_separation_x = 0.04 * x_lim

    if ymax is not None:
        ymax = float(ymax)
        if ymax <= 0:
            raise ValueError("ymax must be > 0 when provided")

    if y_break is None:
        use_break = False
    else:
        y_break = float(y_break)
        use_break = y_max > y_break

    colors = {
        "Not significant": "#C9C9C9",
        "Increased": "#C93C37",
        "Decreased": "#2B6CB0",
    }
    order = ["Not significant", "Decreased", "Increased"]

    if use_break:
        fig, (ax_top, ax) = plt.subplots(
            2,
            1,
            sharex=True,
            figsize=figsize,
            dpi=dpi,
            gridspec_kw={"height_ratios": [upper_height_ratio, lower_height_ratio]},
            constrained_layout=False,
        )
        fig.subplots_adjust(hspace=hspace)
        plot_axes = [ax, ax_top]
    else:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, constrained_layout=False)
        ax_top = None
        plot_axes = [ax]

    for axis in plot_axes:
        for group in order:
            sub = data[data["direction"] == group]
            if sub.empty:
                continue

            is_sig_group = group != "Not significant"

            axis.scatter(
                sub[coef_col],
                sub["neglog10_q"],
                s=point_size_sig if is_sig_group else point_size_ns,
                alpha=alpha_sig if is_sig_group else alpha_ns,
                c=colors[group],
                edgecolors="white" if is_sig_group else "none",
                linewidths=0.35 if is_sig_group else 0.0,
                label=group if axis is ax else "_nolegend_",
                zorder=3 if is_sig_group else 2,
            )

        axis.axvline(effect_thresh, linestyle="--", linewidth=1.0, color="0.45", zorder=1)
        axis.axvline(-effect_thresh, linestyle="--", linewidth=1.0, color="0.45", zorder=1)
        axis.axhline(y_thresh_line, linestyle="--", linewidth=1.0, color="0.45", zorder=1)

    ax.set_xlim(-x_lim, x_lim)

    if use_break:
        upper_min = y_break + break_gap
        upper_max = ymax if ymax is not None else y_max * 1.04

        if upper_max <= upper_min:
            upper_max = upper_min * 1.05

        ax.set_ylim(0, y_break)
        ax_top.set_ylim(upper_min, upper_max)
    else:
        y_lim = ymax if ymax is not None else max(y_thresh_line * 1.15, y_max * 1.06)
        ax.set_ylim(0, y_lim)

    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(r"$-\log_{10}(q)$", fontsize=11)

    ax.grid(axis="y", alpha=0.14, linewidth=0.8)
    ax.grid(axis="x", visible=False)

    if legend_position == "default":
        resolved_legend_loc = "upper right"
        resolved_legend_bbox_to_anchor = None
        resolved_legend_borderaxespad = (
            0.4 if legend_borderaxespad is None else legend_borderaxespad
        )
        stats_text_y = 0.99
    elif legend_position == "top-left":
        resolved_legend_loc = "upper left"
        resolved_legend_bbox_to_anchor = (
            (0.0, 1.0)
            if legend_bbox_to_anchor is None
            else legend_bbox_to_anchor
        )
        resolved_legend_borderaxespad = (
            0.0 if legend_borderaxespad is None else legend_borderaxespad
        )
        stats_text_y = 0.82
    else:
        resolved_legend_loc = "upper right" if legend_loc is None else legend_loc
        resolved_legend_bbox_to_anchor = legend_bbox_to_anchor
        resolved_legend_borderaxespad = (
            0.4 if legend_borderaxespad is None else legend_borderaxespad
        )
        stats_text_y = 0.99

    ax.text(
        0.01,
        stats_text_y,
        f"n = {len(data):,}\nsignificant = {int(data['is_sig'].sum()):,}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
    )

    legend_kwargs = {
        "frameon": False,
        "loc": resolved_legend_loc,
        "fontsize": 9,
        "handletextpad": 0.4,
        "borderaxespad": resolved_legend_borderaxespad,
    }
    if resolved_legend_bbox_to_anchor is not None:
        legend_kwargs["bbox_to_anchor"] = resolved_legend_bbox_to_anchor

    ax.legend(**legend_kwargs)

    if use_break:
        ax.spines["top"].set_visible(False)
        ax_top.spines["bottom"].set_visible(False)

        ax.tick_params(axis="x", top=False)
        ax_top.tick_params(axis="x", bottom=False, labelbottom=False)

        ax_top.grid(False)
        ax_top.set_title(title, fontsize=13, fontweight="bold", pad=10)

        if show_break_label:
            ax.text(
                0.99,
                1.03,
                f"y-axis truncated at {y_break:g}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=8,
                color="0.4",
            )

        sns.despine(ax=ax, offset=5, trim=True, top=True, right=True)
        sns.despine(ax=ax_top, offset=5, trim=True, bottom=True, right=True)
    else:
        ax.set_title(title, fontsize=13, fontweight="bold", pad=10)
        sns.despine(ax=ax, offset=5, trim=True)

    fig.canvas.draw()

    def _select_spaced(sub, x_col, score_col, max_n, min_dx):
        sub = sub.sort_values(score_col, ascending=False).copy()
        chosen = []
        chosen_x = []

        for _, row in sub.iterrows():
            x_value = row[x_col]
            if all(abs(x_value - prev_x) >= min_dx for prev_x in chosen_x):
                chosen.append(row)
                chosen_x.append(x_value)

            if len(chosen) >= max_n:
                break

        if not chosen:
            return sub.head(0).copy()

        return pd.DataFrame(chosen)

    texts_lower = []
    texts_upper = []

    sig = data[data["is_sig"]].copy()
    label_key = pd.Series(
        list(
            zip(
                data[source_col].astype(str),
                data[target_col].astype(str),
                data[label_col].astype(str),
            )
        ),
        index=data.index,
    )

    label_frames = []

    if label_interactions:
        manual_lab = data[label_key.isin(label_interactions)].copy()
        if not manual_lab.empty:
            label_frames.append(manual_lab)

    if not override_auto_labels and top_n > 0 and not sig.empty:
        pos = sig[sig[coef_col] > 0].copy()
        neg = sig[sig[coef_col] < 0].copy()

        pos_lab = _select_spaced(
            pos,
            x_col=coef_col,
            score_col="label_score",
            max_n=max_labels_per_side,
            min_dx=min_label_separation_x,
        )
        neg_lab = _select_spaced(
            neg,
            x_col=coef_col,
            score_col="label_score",
            max_n=max_labels_per_side,
            min_dx=min_label_separation_x,
        )

        auto_lab = pd.concat([neg_lab, pos_lab], axis=0)
        if not auto_lab.empty:
            label_frames.append(auto_lab)

    if label_frames:
        to_label = pd.concat(label_frames, axis=0)
        to_label = to_label.assign(
            _label_key=list(
                zip(
                    to_label[source_col].astype(str),
                    to_label[target_col].astype(str),
                    to_label[label_col].astype(str),
                )
            )
        )
        to_label = to_label.drop_duplicates(subset="_label_key", keep="first")

        for _, row in to_label.iterrows():
            if use_break and row["neglog10_q"] > y_break:
                target_ax = ax_top
                target_texts = texts_upper
            else:
                target_ax = ax
                target_texts = texts_lower

            text = target_ax.text(
                row[coef_col],
                row["neglog10_q"],
                row["plot_label"],
                fontsize=7.8,
                ha="center",
                va="bottom",
                zorder=4,
                bbox=dict(
                    boxstyle="round,pad=0.18",
                    fc="white",
                    ec="none",
                    alpha=0.2,
                ),
            )
            target_texts.append(text)

    fig.canvas.draw()

    if texts_lower:
        adjust_text(
            texts_lower,
            ax=ax,
            only_move={"points": "y", "text": "xy"},
            autoalign="y",
            force_text=(0.04, 0.20),
            force_points=(0.04, 0.15),
            expand_text=(1.03, 1.16),
            expand_points=(1.02, 1.08),
            lim=300,
            arrowprops=dict(
                arrowstyle="-",
                color="0.35",
                lw=0.6,
                shrinkA=2,
                shrinkB=2,
                alpha=0.75,
            ),
        )

    if use_break and texts_upper:
        adjust_text(
            texts_upper,
            ax=ax_top,
            only_move={"points": "y", "text": "xy"},
            autoalign="y",
            force_text=(0.04, 0.18),
            force_points=(0.03, 0.12),
            expand_text=(1.02, 1.10),
            expand_points=(1.01, 1.05),
            lim=150,
            arrowprops=dict(
                arrowstyle="-",
                color="0.35",
                lw=0.6,
                shrinkA=2,
                shrinkB=2,
                alpha=0.75,
            ),
        )

    fig.canvas.draw()

    if use_break:
        return fig, (ax, ax_top), data

    return fig, ax, data
