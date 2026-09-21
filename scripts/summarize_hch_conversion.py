"""Derive per-cell bisulfite conversion from BISCUIT pileup methylation averages.

In NOMe-seq the GpC methyltransferase methylates GpC, so CpH-based conversion
estimates are contaminated by the accessibility signal. HCH is the only class
that is neither CpG nor GpC, which makes it the usable conversion proxy.
"""

import argparse
from pathlib import Path

import pandas as pd


AVERAGE_SUFFIX = "_meth_average.tsv"
# Written by "biscuit pileup -N -w". The H-prefixed columns only exist in
# NOMe-seq mode; their absence means the pileup ran without -N.
REQUIRED_COLUMNS = ("HCHn", "HCHb")
REPORTED = {"HCG": "cpg_meth_pct", "GC": "gpc_meth_pct"}


def weighted_beta(table: pd.DataFrame, count_column: str, beta_column: str) -> float:
    """Average a per-chromosome beta by its observation count."""
    counts = pd.to_numeric(table[count_column], errors="coerce").fillna(0)
    betas = pd.to_numeric(table[beta_column], errors="coerce").fillna(0)
    total = counts.sum()
    return float((betas * counts).sum() / total) if total else float("nan")


def parse_meth_average(path: Path) -> dict:
    """Reduce one per-chromosome table to the scalars used for QC."""
    table = pd.read_csv(path, sep="\t")
    missing = [column for column in REQUIRED_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(
            f"{path} lacks {', '.join(missing)}; run biscuit pileup with -N so it "
            f"reports NOMe-seq contexts. Found: {list(table.columns)}"
        )
    record = {"sample": path.name.removesuffix(AVERAGE_SUFFIX)}
    hch_beta = weighted_beta(table, "HCHn", "HCHb")
    record["hch_meth_pct"] = round(hch_beta, 4)
    record["conversion_rate_pct"] = round(100 - hch_beta, 4)
    record["hch_sites"] = int(pd.to_numeric(table["HCHn"], errors="coerce").fillna(0).sum())
    for prefix, column in REPORTED.items():
        if {f"{prefix}n", f"{prefix}b"} <= set(table.columns):
            record[column] = round(weighted_beta(table, f"{prefix}n", f"{prefix}b"), 4)
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tables", type=Path, nargs="+", help="*_meth_average.tsv files")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output TSV")
    args = parser.parse_args()

    table = pd.DataFrame([parse_meth_average(path) for path in args.tables])
    columns = ["sample", "conversion_rate_pct", "hch_meth_pct", "hch_sites", *REPORTED.values()]
    table = table.reindex(columns=columns)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.sort_values("sample").to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
