from typing import Optional, Sequence
import pandas as pd


def filter_by(
	df: pd.DataFrame,
	sources: Optional[Sequence[str]] = ("T cells", "Endothelial cells"),
	targets: Optional[Sequence[str]] = ("T cells", "Endothelial cells"),
) -> pd.DataFrame:
	"""Remove rows matching selected source and target labels.

	Rows are removed when the ``source`` value is present in ``sources`` or the
	``target`` value is present in ``targets``. If neither filter list is
	provided, a copy of ``df`` is returned unchanged.
	"""
	if not isinstance(df, pd.DataFrame):
		raise TypeError("df must be a pandas DataFrame.")

	missing_cols = [column for column in ("source", "target") if column not in df.columns]
	if missing_cols:
		raise KeyError(
			"df must contain 'source' and 'target' columns. "
			f"Missing: {', '.join(missing_cols)}"
		)

	if sources is None and targets is None:
		return df.copy()

	mask = pd.Series(False, index=df.index)

	if sources is not None:
		source_values = {str(value) for value in sources}
		mask |= df["source"].astype(str).isin(source_values)

	if targets is not None:
		target_values = {str(value) for value in targets}
		mask |= df["target"].astype(str).isin(target_values)

	return df.loc[~mask].copy()