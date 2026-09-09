"""Plot observed CpG counts against global methylation from Methscan cell stats."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def plot_cell_stats(cell_stats_csv: Path, output: Path) -> None:
    """Save one point per cell, with methylation expressed as a percentage."""
    cell_stats = pd.read_csv(cell_stats_csv)
    required = {"global_meth_frac", "n_obs"}
    missing = required.difference(cell_stats.columns)
    if missing:
        raise ValueError(f"Missing columns in {cell_stats_csv}: {', '.join(sorted(missing))}")

    methylation_percent = pd.to_numeric(cell_stats["global_meth_frac"]) * 100
    observed_cpgs = pd.to_numeric(cell_stats["n_obs"])
    output.parent.mkdir(parents=True, exist_ok=True)
    with plt.style.context("ggplot"):
        fig, ax = plt.subplots(figsize=(6, 4.5), layout="constrained")
        try:
            ax.scatter(methylation_percent, observed_cpgs, color="black", s=12)
            ax.set_xlabel("global DNA methylation %")
            ax.set_ylabel("# of observed CpG sites")
            fig.savefig(output, dpi=300)
        finally:
            plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cell_stats_csv", type=Path, help="Methscan cell_stats.csv")
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("cell_stats.png"),
        help="Output image (default: cell_stats.png)",
    )
    args = parser.parse_args()
    plot_cell_stats(args.cell_stats_csv, args.output)


if __name__ == "__main__":
    main()
