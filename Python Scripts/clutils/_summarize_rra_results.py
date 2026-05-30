"""Human-readable summaries for robust rank aggregation interaction results."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import importlib.util
from pathlib import Path
import re
from typing import Union

import numpy as np
import pandas as pd


PathLike = Union[str, Path]

_CELL_GROUP_CAR_CYTOTOXIC_MARKERS = ("cytotoxic", "cd8", "t cell")
_CELL_GROUP_TNK_MARKERS = (
    "cytotoxic",
    "cd8",
    "nk",
    "trm",
    "\u03b3\u03b4",
    "gammadelta",
    "treg",
    "th17",
    "cd4",
    "tfh",
    "naive/central",
    "naive/effector",
)
_CELL_GROUP_MYELOID_MARKERS = (
    "macrophage",
    "monocyte",
    "cdc",
    "dc1",
    "dc2",
    "granulocyte",
    "mast",
)
_CELL_GROUP_STROMAL_MARKERS = (
    "fibroblast",
    "telocyte",
    "myofibroblast",
    "pericyte",
)
_CELL_GROUP_ENDOTHELIAL_MARKERS = (
    "endothelium",
    "endothelial",
    "venous",
    "arterial",
    "arteriolar",
    "capillary",
    "lymphatic",
)
_CELL_GROUP_EPITHELIAL_MARKERS = (
    "ta cells",
    "crypt",
    "absorptive",
    "goblet",
    "colonocyte",
    "epithelium",
    "epithelial",
)
_CELL_GROUP_TNK_LABELS = {"T/NK-lineage", "CAR+ cytotoxic/T"}
_CYTOTOXIC_TARGET_PATTERN = r"cytotoxic|CD8|NK|Trm|CAR"


def _contains_any(text: str, tokens: Sequence[str]) -> bool:
    """Return whether a normalized label contains any marker token."""
    return any(token in text for token in tokens)


def _percent(numerator, denominator):
    """Compute percentages while avoiding divide-by-zero artifacts."""
    if np.isscalar(denominator):
        return 100 * numerator / denominator if denominator else np.nan

    safe_denominator = denominator.replace(0, np.nan)
    return 100 * numerator / safe_denominator


def _require_columns(df: pd.DataFrame, columns: Sequence[str]) -> None:
    """Raise a clear error when required columns are missing."""
    missing_columns = [column for column in columns if column not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")


def _coerce_numeric_columns(df: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    """Return a copy with known numeric summary columns coerced to numeric."""
    out = df.copy()
    for column in columns:
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce")
    return out


def _resolve_excel_writer_engine(engine: str | None) -> str:
    """Choose an installed pandas Excel writer engine for xlsx output."""
    if engine is not None:
        return engine

    for candidate in ("xlsxwriter", "openpyxl"):
        if importlib.util.find_spec(candidate) is not None:
            return candidate

    raise ImportError(
        "Writing .xlsx files requires either 'xlsxwriter' or 'openpyxl'. "
        "Install one of those packages or pass a supported engine explicitly."
    )


def _make_excel_sheet_name(name: str, used_names: set[str]) -> str:
    """Normalize sheet names so each summary key maps to a valid Excel tab."""
    cleaned_name = re.sub(r"[\\/*?:\[\]]", "_", str(name)).strip().strip("'")
    cleaned_name = cleaned_name[:31] or "Sheet"

    candidate = cleaned_name
    suffix = 1
    while candidate in used_names:
        suffix_text = f"_{suffix}"
        candidate = f"{cleaned_name[: 31 - len(suffix_text)]}{suffix_text}"
        suffix += 1

    used_names.add(candidate)
    return candidate


def _build_overall_summary(
    sig: pd.DataFrame,
    *,
    fdr_threshold: float,
    total_all: int,
    median_rank_col: str,
    mean_rank_col: str,
    n_samples_col: str,
    top10_col: str,
) -> pd.DataFrame:
    """Summarize overall signal prevalence and rank characteristics."""
    total_sig = len(sig)
    overall_rows = [
        ("Total RRA interactions in input", total_all, 100.0),
        (
            f"FDR-significant interactions, FDR < {fdr_threshold}",
            total_sig,
            _percent(total_sig, total_all),
        ),
    ]

    if total_sig == 0:
        return pd.DataFrame(overall_rows, columns=["metric", "value", "percent"])

    overall_rows.extend(
        [
            (
                "Median median-rank percentile among significant",
                sig[median_rank_col].median(),
                np.nan,
            ),
            (
                "Mean median-rank percentile among significant",
                sig[median_rank_col].mean(),
                np.nan,
            ),
            (
                "Median mean-rank percentile among significant",
                sig[mean_rank_col].median(),
                np.nan,
            ),
        ]
    )

    if n_samples_col in sig.columns:
        overall_rows.append(
            (
                "Median samples observed among significant",
                sig[n_samples_col].median(),
                np.nan,
            )
        )

    if top10_col in sig.columns:
        overall_rows.append(
            (
                "Median top-10 frequency among significant",
                sig[top10_col].median(),
                np.nan,
            )
        )

    overall_rows.extend(
        [
            (
                "T/NK-lineage target interactions",
                int(sig["is_tnk_target"].sum()),
                _percent(sig["is_tnk_target"].sum(), total_sig),
            ),
            (
                "Cytotoxic/CD8/NK-like target interactions",
                int(sig["is_cytotoxic_target"].sum()),
                _percent(sig["is_cytotoxic_target"].sum(), total_sig),
            ),
            (
                "CAR-positive target interactions",
                int(sig["is_carpos_target"].sum()),
                _percent(sig["is_carpos_target"].sum(), total_sig),
            ),
            (
                "CD74-axis interactions",
                int(sig["is_cd74_axis"].sum()),
                _percent(sig["is_cd74_axis"].sum(), total_sig),
            ),
            (
                "MIF-CD74-CXCR4 interactions",
                int(sig["is_mif_cd74_cxcr4"].sum()),
                _percent(sig["is_mif_cd74_cxcr4"].sum(), total_sig),
            ),
        ]
    )

    return pd.DataFrame(overall_rows, columns=["metric", "value", "percent"])


def _build_key_stats(sig: pd.DataFrame, *, total_all: int) -> pd.DataFrame:
    """Return compact numeric facts that are convenient for manuscript text."""
    total_sig = len(sig)

    return pd.DataFrame(
        {
            "statistic": [
                "n_significant",
                "percent_significant_of_input",
                "n_tnk_targets",
                "percent_tnk_targets",
                "n_cytotoxic_targets",
                "percent_cytotoxic_targets",
                "n_carpos_targets",
                "n_cd74_axis",
                "percent_cd74_axis",
                "n_mif_cd74_cxcr4",
                "percent_mif_cd74_cxcr4",
            ],
            "value": [
                total_sig,
                _percent(total_sig, total_all),
                int(sig["is_tnk_target"].sum()),
                _percent(sig["is_tnk_target"].sum(), total_sig),
                int(sig["is_cytotoxic_target"].sum()),
                _percent(sig["is_cytotoxic_target"].sum(), total_sig),
                int(sig["is_carpos_target"].sum()),
                int(sig["is_cd74_axis"].sum()),
                _percent(sig["is_cd74_axis"].sum(), total_sig),
                int(sig["is_mif_cd74_cxcr4"].sum()),
                _percent(sig["is_mif_cd74_cxcr4"].sum(), total_sig),
            ],
        }
    )


def _build_example_columns(
    sig: pd.DataFrame,
    *,
    source_col: str,
    target_col: str,
    interaction_col: str,
    pval_col: str,
    fdr_col: str,
    median_rank_col: str,
    mean_rank_col: str,
    n_samples_col: str,
    top10_col: str,
    logprob_col: str,
) -> list[str]:
    """Return a stable set of columns for example interaction tables."""
    columns = [
        source_col,
        target_col,
        interaction_col,
        "interaction_family",
        "source_group",
        "target_group",
        pval_col,
        fdr_col,
        median_rank_col,
        mean_rank_col,
    ]
    columns = [column for column in columns if column in sig.columns]

    for optional_column in (top10_col, n_samples_col, logprob_col):
        if optional_column in sig.columns:
            columns.append(optional_column)

    return columns


def classify_cell_group(cell: object) -> str:
    """Assign a coarse biological cell-group label used in RRA summaries."""
    if pd.isna(cell):
        return "Unknown"

    label = str(cell).lower()

    if "car" in label and _contains_any(label, _CELL_GROUP_CAR_CYTOTOXIC_MARKERS):
        return "CAR+ cytotoxic/T"

    if _contains_any(label, _CELL_GROUP_TNK_MARKERS):
        return "T/NK-lineage"

    if _contains_any(label, _CELL_GROUP_MYELOID_MARKERS):
        return "Myeloid/mast"

    if _contains_any(label, _CELL_GROUP_STROMAL_MARKERS):
        return "Stromal/fibroblast"

    if _contains_any(label, _CELL_GROUP_ENDOTHELIAL_MARKERS):
        return "Endothelial/vascular"

    if _contains_any(label, _CELL_GROUP_EPITHELIAL_MARKERS):
        return "Epithelial/crypt"

    if "glial" in label or "glia" in label:
        return "Glial"

    return "Other"


def classify_interaction_family(interaction: object) -> str:
    """Assign a coarse ligand-receptor family label used in high-level summaries."""
    if pd.isna(interaction):
        return "Unknown"

    label = str(interaction).upper()

    if "CD74" in label:
        if "MIF" in label and "CXCR4" in label:
            return "MIF-CD74-CXCR4"
        if "MIF" in label and "CD44" in label:
            return "MIF-CD74-CD44"
        if "APP" in label:
            return "APP-CD74"
        return "CD74-associated"

    if label.startswith("HLA") or "_HLA" in label or "B2M" in label:
        return "MHC-I / HLA"

    if label.startswith("COL") or "_COL" in label:
        return "Collagen/ECM"

    if label.startswith("LAM") or "_LAM" in label:
        return "Laminin/ECM"

    if _contains_any(label, ("FN1", "SPP1", "ICAM", "ITGA", "ITGB", "SDC", "CD44")):
        return "Adhesion / integrin / matrix"

    if (
        label.startswith("CXCL")
        or "_CXCL" in label
        or label.startswith("CCL")
        or "_CCL" in label
    ):
        return "Chemokine"

    if label.startswith("MIF") or "_MIF" in label:
        return "MIF-associated"

    if (
        label.startswith("MDK")
        or "_MDK" in label
        or label.startswith("MK")
        or "_MK" in label
    ):
        return "MDK/MK"

    if _contains_any(label, ("TNF", "TNFSF", "TNFRSF")):
        return "TNF/TNFR"

    if _contains_any(label, ("TGFB", "TGFBR", "ACVR")):
        return "TGF-beta"

    if "JAM" in label:
        return "JAM"

    if _contains_any(label, ("NECTIN", "TIGIT", "LILRB", "KLR")):
        return "Checkpoint / NK receptor"

    if _contains_any(label, ("PGE", "PTGER", "PTGES")):
        return "Prostaglandin"

    if _contains_any(label, ("VEGF", "FLT", "KDR")):
        return "VEGF"

    if _contains_any(label, ("IGFBP", "IGF")):
        return "IGF/IGFBP"

    return "Other"


def add_rra_annotations(
    df: pd.DataFrame,
    source_col: str = "source",
    target_col: str = "target",
    interaction_col: str = "interaction_name",
) -> pd.DataFrame:
    """Add coarse cell-group and interaction-family annotations used in summaries."""
    _require_columns(df, [source_col, target_col, interaction_col])

    out = df.copy()
    out["source_group"] = out[source_col].map(classify_cell_group)
    out["target_group"] = out[target_col].map(classify_cell_group)
    out["interaction_family"] = out[interaction_col].map(classify_interaction_family)

    target_text = out[target_col].fillna("").astype(str)
    interaction_upper = out[interaction_col].fillna("").astype(str).str.upper()

    out["is_tnk_target"] = out["target_group"].isin(_CELL_GROUP_TNK_LABELS)
    out["is_carpos_target"] = target_text.str.contains("CAR", case=False, na=False)
    out["is_cytotoxic_target"] = target_text.str.contains(
        _CYTOTOXIC_TARGET_PATTERN,
        case=False,
        na=False,
        regex=True,
    )
    out["is_mif_cd74_cxcr4"] = interaction_upper.eq("MIF_CD74_CXCR4")
    out["is_cd74_axis"] = interaction_upper.str.contains("CD74", na=False)

    return out


def summarize_group(
    df: pd.DataFrame,
    group_cols: str | Sequence[str],
    fdr_col: str = "rra_fdr",
    median_rank_col: str = "median_rank_pct",
    mean_rank_col: str = "mean_rank_pct",
    n_samples_col: str = "n_samples_observed",
    top10_col: str = "top10_freq",
    total_n: int | None = None,
    top_n: int | None = 30,
) -> pd.DataFrame:
    """Aggregate RRA summary metrics by one or more grouping columns."""
    group_keys = [group_cols] if isinstance(group_cols, str) else list(group_cols)
    if not group_keys:
        raise ValueError("group_cols must contain at least one column name.")

    if total_n is None:
        total_n = len(df)

    aggregation = {
        "n_interactions": (fdr_col, "size"),
        "min_fdr": (fdr_col, "min"),
        "median_fdr": (fdr_col, "median"),
        "median_rank_pct_median": (median_rank_col, "median"),
        "median_rank_pct_mean": (median_rank_col, "mean"),
        "mean_rank_pct_median": (mean_rank_col, "median"),
    }

    if n_samples_col in df.columns:
        aggregation["median_n_samples_observed"] = (n_samples_col, "median")
        aggregation["min_n_samples_observed"] = (n_samples_col, "min")

    if top10_col in df.columns:
        aggregation["median_top10_freq"] = (top10_col, "median")
        aggregation["mean_top10_freq"] = (top10_col, "mean")

    summary = (
        df.groupby(group_keys, dropna=False)
        .agg(**aggregation)
        .reset_index()
        .sort_values(["n_interactions", "median_rank_pct_median"], ascending=[False, False])
    )
    summary["percent_of_significant"] = _percent(summary["n_interactions"], total_n)

    if top_n is not None:
        summary = summary.head(top_n)

    return summary


def summarize_rra_results(
    rra_df: pd.DataFrame,
    fdr_threshold: float = 0.05,
    source_col: str = "source",
    target_col: str = "target",
    interaction_col: str = "interaction_name",
    fdr_col: str = "rra_fdr",
    pval_col: str = "rra_pval",
    median_rank_col: str = "median_rank_pct",
    mean_rank_col: str = "mean_rank_pct",
    n_samples_col: str = "n_samples_observed",
    top1_col: str = "top1_freq",
    top5_col: str = "top5_freq",
    top10_col: str = "top10_freq",
    logprob_col: str = "median_logprob",
    top_n: int = 30,
) -> dict[str, pd.DataFrame]:
    """Build a compact set of high-level tables from interaction-level RRA results.

    Args:
        rra_df: Interaction-level robust rank aggregation results.
        fdr_threshold: Significance cutoff applied to ``fdr_col``.
        source_col: Sender cell-type column.
        target_col: Receiver cell-type column.
        interaction_col: Ligand-receptor interaction label column.
        fdr_col: Multiple-testing adjusted RRA p-value column.
        pval_col: Nominal RRA p-value column.
        median_rank_col: Median within-sample rank percentile column.
        mean_rank_col: Mean within-sample rank percentile column.
        n_samples_col: Sample support column.
        top1_col: Frequency of appearing in the top 1 percent.
        top5_col: Frequency of appearing in the top 5 percent.
        top10_col: Frequency of appearing in the top 10 percent.
        logprob_col: Optional log-probability summary column.
        top_n: Maximum number of rows to retain for ranked summary tables.

    Returns:
        Dictionary of summary tables with stable keys for downstream notebooks,
        plotting helpers, and exports.

    Raises:
        ValueError: If required RRA columns are missing.
    """
    required_columns = [
        source_col,
        target_col,
        interaction_col,
        fdr_col,
        median_rank_col,
        mean_rank_col,
    ]
    _require_columns(rra_df, required_columns)

    numeric_columns = [
        fdr_col,
        pval_col,
        median_rank_col,
        mean_rank_col,
        n_samples_col,
        top1_col,
        top5_col,
        top10_col,
        logprob_col,
    ]
    df = _coerce_numeric_columns(rra_df, numeric_columns)
    df = add_rra_annotations(
        df,
        source_col=source_col,
        target_col=target_col,
        interaction_col=interaction_col,
    )

    sig = df.loc[df[fdr_col] < fdr_threshold].copy()
    total_all = len(df)
    total_sig = len(sig)

    source_summary = summarize_group(
        sig,
        source_col,
        fdr_col,
        median_rank_col,
        mean_rank_col,
        n_samples_col,
        top10_col,
        total_n=total_sig,
        top_n=top_n,
    )
    target_summary = summarize_group(
        sig,
        target_col,
        fdr_col,
        median_rank_col,
        mean_rank_col,
        n_samples_col,
        top10_col,
        total_n=total_sig,
        top_n=top_n,
    )
    source_group_summary = summarize_group(
        sig,
        "source_group",
        fdr_col,
        median_rank_col,
        mean_rank_col,
        n_samples_col,
        top10_col,
        total_n=total_sig,
        top_n=None,
    )
    target_group_summary = summarize_group(
        sig,
        "target_group",
        fdr_col,
        median_rank_col,
        mean_rank_col,
        n_samples_col,
        top10_col,
        total_n=total_sig,
        top_n=None,
    )
    source_target_summary = summarize_group(
        sig,
        [source_col, target_col],
        fdr_col,
        median_rank_col,
        mean_rank_col,
        n_samples_col,
        top10_col,
        total_n=total_sig,
        top_n=top_n,
    )
    source_group_target_group_summary = summarize_group(
        sig,
        ["source_group", "target_group"],
        fdr_col,
        median_rank_col,
        mean_rank_col,
        n_samples_col,
        top10_col,
        total_n=total_sig,
        top_n=top_n,
    )
    interaction_summary = summarize_group(
        sig,
        interaction_col,
        fdr_col,
        median_rank_col,
        mean_rank_col,
        n_samples_col,
        top10_col,
        total_n=total_sig,
        top_n=top_n,
    )
    interaction_family_summary = summarize_group(
        sig,
        "interaction_family",
        fdr_col,
        median_rank_col,
        mean_rank_col,
        n_samples_col,
        top10_col,
        total_n=total_sig,
        top_n=None,
    )

    cytotoxic_target_df = sig.loc[sig["is_cytotoxic_target"]].copy()
    carpos_target_df = sig.loc[sig["is_carpos_target"]].copy()
    mif_cd74_cxcr4_df = sig.loc[sig["is_mif_cd74_cxcr4"]].copy()

    cytotoxic_target_summary_by_interaction = summarize_group(
        cytotoxic_target_df,
        interaction_col,
        fdr_col,
        median_rank_col,
        mean_rank_col,
        n_samples_col,
        top10_col,
        total_n=len(cytotoxic_target_df),
        top_n=top_n,
    )
    cytotoxic_target_summary_by_source = summarize_group(
        cytotoxic_target_df,
        source_col,
        fdr_col,
        median_rank_col,
        mean_rank_col,
        n_samples_col,
        top10_col,
        total_n=len(cytotoxic_target_df),
        top_n=top_n,
    )
    carpos_target_summary_by_interaction = summarize_group(
        carpos_target_df,
        interaction_col,
        fdr_col,
        median_rank_col,
        mean_rank_col,
        n_samples_col,
        top10_col,
        total_n=len(carpos_target_df),
        top_n=top_n,
    )
    carpos_target_summary_by_source = summarize_group(
        carpos_target_df,
        source_col,
        fdr_col,
        median_rank_col,
        mean_rank_col,
        n_samples_col,
        top10_col,
        total_n=len(carpos_target_df),
        top_n=top_n,
    )
    mif_cd74_cxcr4_by_target = summarize_group(
        mif_cd74_cxcr4_df,
        target_col,
        fdr_col,
        median_rank_col,
        mean_rank_col,
        n_samples_col,
        top10_col,
        total_n=len(mif_cd74_cxcr4_df),
        top_n=top_n,
    )
    mif_cd74_cxcr4_by_source = summarize_group(
        mif_cd74_cxcr4_df,
        source_col,
        fdr_col,
        median_rank_col,
        mean_rank_col,
        n_samples_col,
        top10_col,
        total_n=len(mif_cd74_cxcr4_df),
        top_n=top_n,
    )
    mif_cd74_cxcr4_source_target = summarize_group(
        mif_cd74_cxcr4_df,
        [source_col, target_col],
        fdr_col,
        median_rank_col,
        mean_rank_col,
        n_samples_col,
        top10_col,
        total_n=len(mif_cd74_cxcr4_df),
        top_n=top_n,
    )

    example_columns = _build_example_columns(
        sig,
        source_col=source_col,
        target_col=target_col,
        interaction_col=interaction_col,
        pval_col=pval_col,
        fdr_col=fdr_col,
        median_rank_col=median_rank_col,
        mean_rank_col=mean_rank_col,
        n_samples_col=n_samples_col,
        top10_col=top10_col,
        logprob_col=logprob_col,
    )

    top_by_fdr = sig.sort_values(fdr_col, ascending=True).loc[:, example_columns].head(top_n)
    top_by_median_rank = (
        sig.sort_values(median_rank_col, ascending=False).loc[:, example_columns].head(top_n)
    )
    top_mif_cd74_cxcr4 = (
        mif_cd74_cxcr4_df.sort_values([median_rank_col, fdr_col], ascending=[False, True])
        .loc[:, example_columns]
        .head(top_n)
    )
    top_carpos_incoming = (
        carpos_target_df.sort_values([median_rank_col, fdr_col], ascending=[False, True])
        .loc[:, example_columns]
        .head(top_n)
    )

    return {
        "significant_interactions": sig,
        "overall": _build_overall_summary(
            sig,
            fdr_threshold=fdr_threshold,
            total_all=total_all,
            median_rank_col=median_rank_col,
            mean_rank_col=mean_rank_col,
            n_samples_col=n_samples_col,
            top10_col=top10_col,
        ),
        "key_stats": _build_key_stats(sig, total_all=total_all),
        "source_summary": source_summary,
        "target_summary": target_summary,
        "source_group_summary": source_group_summary,
        "target_group_summary": target_group_summary,
        "source_target_summary": source_target_summary,
        "source_group_target_group_summary": source_group_target_group_summary,
        "interaction_summary": interaction_summary,
        "interaction_family_summary": interaction_family_summary,
        "cytotoxic_target_summary_by_interaction": cytotoxic_target_summary_by_interaction,
        "cytotoxic_target_summary_by_source": cytotoxic_target_summary_by_source,
        "carpos_target_summary_by_interaction": carpos_target_summary_by_interaction,
        "carpos_target_summary_by_source": carpos_target_summary_by_source,
        "mif_cd74_cxcr4_by_target": mif_cd74_cxcr4_by_target,
        "mif_cd74_cxcr4_by_source": mif_cd74_cxcr4_by_source,
        "mif_cd74_cxcr4_source_target": mif_cd74_cxcr4_source_target,
        "top_by_fdr": top_by_fdr,
        "top_by_median_rank": top_by_median_rank,
        "top_mif_cd74_cxcr4": top_mif_cd74_cxcr4,
        "top_carpos_incoming": top_carpos_incoming,
    }


def print_rra_report(summaries: Mapping[str, pd.DataFrame], top_n: int = 10) -> None:
    """Print a compact console report from summarized RRA result tables."""
    sections = [
        ("OVERALL", "overall", None),
        ("KEY STATS", "key_stats", None),
        ("TOP TARGETS", "target_summary", top_n),
        ("TOP SOURCES", "source_summary", top_n),
        ("TOP SOURCE -> TARGET PAIRS", "source_target_summary", top_n),
        ("TOP INTERACTIONS", "interaction_summary", top_n),
        ("INTERACTION FAMILIES", "interaction_family_summary", None),
        (
            "CYTOTOXIC/CD8/NK TARGETS: TOP INTERACTIONS",
            "cytotoxic_target_summary_by_interaction",
            top_n,
        ),
        ("CAR+ TARGETS: TOP INCOMING SOURCES", "carpos_target_summary_by_source", top_n),
        (
            "MIF-CD74-CXCR4 SOURCE -> TARGET PAIRS",
            "mif_cd74_cxcr4_source_target",
            top_n,
        ),
    ]

    for title, key, limit in sections:
        print(f"\n=== {title} ===")
        table = summaries[key]
        if limit is not None:
            table = table.head(limit)
        print(table.to_string(index=False))


def export_rra_summary_to_excel(
    summaries: Mapping[str, pd.DataFrame],
    output_path: PathLike = "rra_summary_tables.xlsx",
    exclude_large: bool = True,
    include_index: bool = False,
    engine: str | None = None,
) -> Path:
    """Export RRA summary tables to a multi-sheet Excel workbook."""
    if not summaries:
        raise ValueError("summaries must contain at least one dataframe.")

    workbook_path = Path(output_path)
    workbook_path.parent.mkdir(parents=True, exist_ok=True)

    used_sheet_names: set[str] = set()
    writer_kwargs = {"engine": _resolve_excel_writer_engine(engine)}

    with pd.ExcelWriter(workbook_path, **writer_kwargs) as writer:
        for name, table in summaries.items():
            if exclude_large and name == "significant_interactions":
                continue

            if not isinstance(table, pd.DataFrame):
                raise TypeError(
                    "Each summary must be a pandas DataFrame. "
                    f"Key '{name}' has type {type(table).__name__}."
                )

            sheet_name = _make_excel_sheet_name(name, used_sheet_names)
            table.to_excel(writer, sheet_name=sheet_name, index=include_index)

    return workbook_path


__all__ = [
    "add_rra_annotations",
    "classify_cell_group",
    "classify_interaction_family",
    "export_rra_summary_to_excel",
    "print_rra_report",
    "summarize_group",
    "summarize_rra_results",
]