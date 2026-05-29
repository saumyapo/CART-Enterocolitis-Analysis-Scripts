"""Helpers for writing dataframes into Excel workbooks."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Union

import pandas as pd


PathLike = Union[str, Path]


def _resolve_excel_writer_engine() -> str:
    """Return the engine required for replacing or appending worksheets."""
    if importlib.util.find_spec("openpyxl") is None:
        raise ImportError(
            "Updating .xlsx workbooks requires 'openpyxl'. "
            "Install it and retry."
        )

    return "openpyxl"


def write_dataframe_to_excel(
    df: pd.DataFrame,
    output_file: PathLike,
    sheet_name: str,
    include_index: bool = False,
) -> Path:
    """Write a dataframe to an Excel tab, creating or replacing as needed.

    If the workbook does not exist, it is created. If it exists and the target
    sheet already exists, that sheet is replaced. Otherwise the new sheet is
    appended to the workbook.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    writer_kwargs = {
        "engine": _resolve_excel_writer_engine(),
        "mode": "a" if output_path.exists() else "w",
    }
    if output_path.exists():
        writer_kwargs["if_sheet_exists"] = "replace"

    with pd.ExcelWriter(output_path, **writer_kwargs) as writer:
        df.to_excel(writer, sheet_name=str(sheet_name), index=include_index)

    return output_path