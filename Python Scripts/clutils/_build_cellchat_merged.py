"""Builders for merged CellChat input tables."""

from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path
from typing import Mapping, Optional, Sequence, Union
from zipfile import ZipFile

import numpy as np
import pandas as pd


PathLike = Union[str, Path]

DEFAULT_CATEGORICAL_COLUMNS = (
    "source",
    "target",
    "ligand",
    "receptor",
    "interaction_name",
    "interaction_name_2",
    "pathway_name",
    "annotation",
    "evidence",
    "sample",
    "condition",
    "chemistry",
)


def _iter_matching_members(zip_file: ZipFile, file_pattern: str) -> Sequence[str]:
    """Return archive members that match the requested CSV pattern.

    Args:
        zip_file: Open zip archive containing per-sample CellChat CSV files.
        file_pattern: Basename glob pattern used to select CSV members.

    Returns:
        Sorted archive member names that should be loaded.
    """
    members = []
    for info in zip_file.infolist():
        member_path = Path(info.filename)
        member_name = member_path.name

        if info.is_dir():
            continue
        if "__MACOSX" in member_path.parts or member_name.startswith("._"):
            continue
        if fnmatch(member_name, file_pattern):
            members.append(info.filename)

    return sorted(members)


def _sample_from_member(member_name: str) -> str:
    """Infer the sample identifier from an archive member path.

    Args:
        member_name: Archive member path for a per-sample CSV.

    Returns:
        Sample identifier parsed from the file name.
    """
    file_name = Path(member_name).name
    return file_name.removeprefix("cellchat_").removesuffix(".csv")


def _normalize_chemistry_values(
    chemistry_by_sample: Mapping[str, str],
) -> Mapping[str, str]:
    """Normalize chemistry labels to the compact values used downstream.

    Args:
        chemistry_by_sample: Sample-to-chemistry mapping.

    Returns:
        Mapping with normalized chemistry values.
    """
    replacements = {
        "3' V3": "3PV3",
        "5' V3": "5PV3",
    }
    return {
        str(sample): replacements.get(str(chemistry), str(chemistry))
        for sample, chemistry in chemistry_by_sample.items()
    }


def build_cellchat_merged(
    zip_file: PathLike = "/home/jupyter/sideproject/colitis/data/raw_v4/per_sample_csvs.zip",
    output_file: PathLike = "/home/jupyter/sideproject/colitis/data/intermediate_v4/cellchat_merged.parquet",
    file_pattern: str = "cellchat_*.csv",
    excluded_populations: Sequence[str] = ("Mixed", "Proliferating"),
    chemistry_by_sample: Optional[Mapping[str, str]] = None,
    default_chemistry: str = "unknown",
    categorical_columns: Sequence[str] = DEFAULT_CATEGORICAL_COLUMNS,
) -> pd.DataFrame:
    """Build the merged CellChat table directly from the raw_v4 zip archive.

    This loader reads all matching per-sample CSV files from the provided zip
    archive without extracting them to disk, concatenates them into a single
    dataframe, applies the same core cleanup as the legacy pipeline, and writes
    the merged result to parquet.

    Args:
        zip_file: Path to the zip archive containing per-sample CellChat CSV
            files.
        output_file: Destination parquet path for the merged dataset.
        file_pattern: Basename glob pattern used to select members from the zip
            archive.
        excluded_populations: Cell populations to remove when they appear in the
            source or target columns.
        chemistry_by_sample: Optional mapping from sample identifier to
            sequencing chemistry. When omitted, the function preserves an input
            chemistry column if present or fills a default value.
        default_chemistry: Fallback chemistry label used when no chemistry is
            available for a sample.
        categorical_columns: Columns to coerce to pandas categorical dtype when
            present after loading.

    Returns:
        The merged, filtered dataframe written to parquet.

    Raises:
        FileNotFoundError: If the zip archive does not exist or contains no
            matching CSV files.
        ValueError: If required columns are missing from the merged dataframe.
    """
    zip_path = Path(zip_file)
    output_path = Path(output_file)

    if not zip_path.exists():
        raise FileNotFoundError(f"Zip archive does not exist: {zip_path}")

    with ZipFile(zip_path) as archive:
        members = _iter_matching_members(archive, file_pattern)
        if not members:
            raise FileNotFoundError(
                f"No files found matching pattern {file_pattern!r} in {zip_path}"
            )

        frames = []
        for member_name in members:
            with archive.open(member_name) as handle:
                frame = pd.read_csv(handle)

            if "sample" not in frame.columns:
                frame["sample"] = _sample_from_member(member_name)

            frames.append(frame)

    df = pd.concat(frames, ignore_index=True)

    required_columns = {"source", "target", "prob", "sample"}
    missing_columns = sorted(required_columns.difference(df.columns))
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    if chemistry_by_sample is not None:
        chemistry_lookup = _normalize_chemistry_values(chemistry_by_sample)
        df["chemistry"] = df["sample"].astype(str).map(chemistry_lookup)
    elif "chemistry" not in df.columns:
        df["chemistry"] = default_chemistry

    df["chemistry"] = df["chemistry"].fillna(default_chemistry)

    if "condition" not in df.columns:
        raise ValueError(
            "Missing required columns: ['condition']. Provide condition in the input CSVs."
        )

    df["prob"] = pd.to_numeric(df["prob"], errors="coerce")
    df["logprob"] = np.log(df["prob"] + 1e-6)

    if excluded_populations:
        excluded = set(excluded_populations)
        mask = ~df["target"].isin(excluded) & ~df["source"].isin(excluded)
        df = df.loc[mask].copy()

    for column in categorical_columns:
        if column in df.columns:
            df[column] = df[column].astype(str).astype("category")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output_path)

    return df