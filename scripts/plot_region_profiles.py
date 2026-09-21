"""Plot per-cell CpG methylation and GpC accessibility around a set of regions.

Both marks share one panel per cell: if the GpC labelling worked, accessibility
rises over the region while CpG methylation falls, and the two curves separate.
A cell whose GpC trace stays flat did not get labelled.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd


PROFILE_COLUMNS = ["position", "cell_name", "meth_frac"]
# methscan names cells after the BED files it was given, so the mark comes back
# attached to the cell name.
DEFAULT_SUFFIXES = ("_HCG", "_GCH", ".hcg", ".gch")
MARK_STYLE = {"CpG": "tab:blue", "GpC": "tab:red"}


def load_profile(path: Path, column: str, bin_size: int, suffixes: tuple) -> pd.DataFrame:
    """Bin one methscan profile so 1-bp resolution becomes a drawable line."""
    table = pd.read_csv(path, usecols=PROFILE_COLUMNS)
    for suffix in suffixes:
        table["cell_name"] = table["cell_name"].str.removesuffix(suffix)
    table["position_binned"] = (table["position"] // bin_size) * bin_size
    binned = (
        table.groupby(["position_binned", "cell_name"], as_index=False)["meth_frac"]
        .mean()
        .rename(columns={"meth_frac": column})
    )
    return binned


def draw_page(cells: list, profile: pd.DataFrame, columns: int, rows: int, region: str):
    """Draw one page of small multiples, one cell per panel."""
    figure, axes = plt.subplots(
        rows, columns, figsize=(2.2 * columns, 1.9 * rows),
        squeeze=False, sharex=True, sharey=True, layout="constrained",
    )
    flat = axes.flatten()
    for axis, cell in zip(flat, cells):
        panel = profile[profile["cell_name"] == cell].sort_values("position_binned")
        position_kb = panel["position_binned"] / 1000
        for mark, color in MARK_STYLE.items():
            if mark in panel.columns:
                axis.plot(position_kb, panel[mark], color=color, linewidth=0.6)
        axis.set_title(cell, fontsize=5)
        axis.set_ylim(0, 1)
        axis.set_yticks([0, 0.5, 1])
        axis.set_yticklabels(["0%", "50%", "100%"], fontsize=5)
        axis.tick_params(axis="x", labelsize=5)
        axis.axvline(0, color="grey", linewidth=0.4, linestyle=":")
    for axis in flat[len(cells):]:
        axis.set_axis_off()
    figure.supxlabel(f"position relative to {region} [kb]", fontsize=8)
    figure.supylabel("methylation / accessibility", fontsize=8)
    handles = [
        plt.Line2D([], [], color=color, label=mark) for mark, color in MARK_STYLE.items()
    ]
    figure.legend(handles=handles, loc="outside upper right", fontsize=7, ncol=2)
    return figure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpg", type=Path, required=True, help="CpG profile CSV")
    parser.add_argument("--gpc", type=Path, required=True, help="GpC profile CSV")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output PDF")
    parser.add_argument("--region", default="TSS", help="Region name used in labels")
    parser.add_argument("--bin-size", type=int, default=10, help="Binning in bp")
    parser.add_argument("--columns", type=int, default=6, help="Panels per row")
    parser.add_argument("--rows", type=int, default=4, help="Panel rows per page")
    parser.add_argument(
        "--strip-suffix", action="append", default=None,
        help="Suffix removed from cell names; repeatable",
    )
    args = parser.parse_args()

    suffixes = tuple(args.strip_suffix) if args.strip_suffix else DEFAULT_SUFFIXES
    profile = load_profile(args.cpg, "CpG", args.bin_size, suffixes).merge(
        load_profile(args.gpc, "GpC", args.bin_size, suffixes),
        on=["position_binned", "cell_name"], how="outer",
    )
    cells = sorted(profile["cell_name"].dropna().unique())
    per_page = args.columns * args.rows
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(args.output) as pdf:
        for start in range(0, len(cells), per_page):
            figure = draw_page(
                cells[start:start + per_page], profile, args.columns, args.rows, args.region
            )
            pdf.savefig(figure)
            plt.close(figure)


if __name__ == "__main__":
    main()
