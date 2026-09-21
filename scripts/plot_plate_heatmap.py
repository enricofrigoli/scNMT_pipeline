"""Lay per-cell QC metrics out on the 384-well plate they were sequenced from.

Spatial artefacts (edge evaporation, pipetting drift, thermal gradients) are
invisible in a scatter plot but obvious on the plate grid.
"""

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Cell names end in a well ID such as ..._A1 or ..._P24.
WELL_PATTERN = re.compile(r"([A-P])(\d{1,2})$")
PLATE_ROWS = [chr(code) for code in range(ord("A"), ord("P") + 1)]
PLATE_COLUMNS = list(range(1, 25))
MIN_PARSED_FRACTION = 0.5


def well_positions(samples: pd.Series) -> pd.DataFrame:
    """Split each cell name into plate coordinates, leaving misses as NA."""
    extracted = samples.str.extract(WELL_PATTERN)
    extracted.columns = ["row", "column"]
    extracted["column"] = pd.to_numeric(extracted["column"], errors="coerce")
    valid = extracted["row"].notna() & extracted["column"].between(1, 24)
    return extracted.where(valid)


def plot_metric(axis, table: pd.DataFrame, metric: str) -> None:
    """Draw one metric as a 16x24 grid with wells that have no value left blank."""
    grid = pd.DataFrame(index=PLATE_ROWS, columns=PLATE_COLUMNS, dtype=float)
    for row, column, value in zip(table["row"], table["column"], table[metric]):
        if pd.notna(row) and pd.notna(column):
            grid.loc[row, int(column)] = value
    image = axis.imshow(grid.to_numpy(dtype=float), cmap="viridis", aspect="auto")
    axis.set_title(metric, fontsize=9)
    axis.set_xticks(range(len(PLATE_COLUMNS)))
    axis.set_xticklabels(PLATE_COLUMNS, fontsize=5)
    axis.set_yticks(range(len(PLATE_ROWS)))
    axis.set_yticklabels(PLATE_ROWS, fontsize=5)
    axis.figure.colorbar(image, ax=axis, fraction=0.025, pad=0.02)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metrics", type=Path, help="Per-cell metrics TSV")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output image")
    parser.add_argument(
        "-m", "--metric", action="append", required=True,
        help="Column to plot; repeat for several panels",
    )
    parser.add_argument("--title", default="Plate QC")
    args = parser.parse_args()

    table = pd.read_csv(args.metrics, sep="\t", comment="#")
    table = pd.concat([table, well_positions(table["sample"].astype(str))], axis=1)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    parsed = table["row"].notna().mean() if len(table) else 0.0
    metrics = [metric for metric in args.metric if metric in table.columns]
    if parsed < MIN_PARSED_FRACTION or not metrics:
        # Sample names outside the plate convention are not an error: the run is
        # still valid, only this view does not apply to it.
        figure, axis = plt.subplots(figsize=(8, 3), layout="constrained")
        reason = (
            f"only {parsed:.0%} of cell names end in a well ID (A1-P24)"
            if parsed < MIN_PARSED_FRACTION
            else f"none of the requested metrics are present: {', '.join(args.metric)}"
        )
        axis.text(0.5, 0.5, f"No plate layout plotted:\n{reason}", ha="center", va="center")
        axis.set_axis_off()
        figure.savefig(args.output, dpi=200)
        plt.close(figure)
        return

    figure, axes = plt.subplots(
        len(metrics), 1, figsize=(9, 2.6 * len(metrics)), squeeze=False, layout="constrained"
    )
    for axis, metric in zip(axes[:, 0], metrics):
        plot_metric(axis, table, metric)
    figure.suptitle(args.title, fontsize=11)
    figure.savefig(args.output, dpi=200)
    plt.close(figure)


if __name__ == "__main__":
    main()
