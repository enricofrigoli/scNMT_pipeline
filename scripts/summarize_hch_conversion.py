"""Derive per-cell bisulfite conversion from BISCUIT pileup methylation averages.

In NOMe-seq the GpC methyltransferase methylates GpC, so CpH-based conversion
estimates are contaminated by the accessibility signal. HCH is the only class
that is neither CpG nor GpC, which makes it the usable conversion proxy.
"""

import argparse
import csv
import math
from pathlib import Path

import pandas as pd


AVERAGE_SUFFIX = "_meth_average.tsv"
# Written by "biscuit pileup -N -w". The H-prefixed columns only exist in
# NOMe-seq mode; their absence means the pileup ran without -N.
REQUIRED_COLUMNS = ("sample", "chrm", "HCHn", "HCHb")
REPORTED = {"HCG": "cpg_meth_pct", "GC": "gpc_meth_pct"}


def parse_count_beta(count: str, beta: str, label: str) -> tuple[int, float]:
    """Read a site count and a percentage without turning missing data into zero."""
    count = count.strip()
    if not count.isascii() or not count.isdecimal():
        raise ValueError(f"{label}: expected a non-negative integer count, got {count!r}")
    sites = int(count)
    beta = beta.strip()
    if beta in ("", "."):
        if sites:
            raise ValueError(f"{label}: missing methylation percentage for {sites} sites")
        return sites, float("nan")
    try:
        # BISCUIT writes percentages, e.g. 2.000%; bare numeric percentages are
        # accepted too, always in the same 0-100 units.
        percentage = float(beta.removesuffix("%"))
    except ValueError as exc:
        raise ValueError(f"{label}: invalid methylation percentage {beta!r}") from exc
    if not math.isfinite(percentage) or not 0 <= percentage <= 100:
        raise ValueError(f"{label}: methylation percentage outside 0-100: {beta!r}")
    return sites, percentage


def read_meth_average(path: Path) -> pd.DataFrame:
    """Validate TSV rows before pandas can infer an index from extra columns."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        columns = next(reader, [])
        if not columns or len(columns) != len(set(columns)):
            raise ValueError(f"{path}: missing or duplicate column names")
        missing = [column for column in REQUIRED_COLUMNS if column not in columns]
        if missing:
            raise ValueError(
                f"{path} lacks {', '.join(missing)}; run biscuit pileup with -N so it "
                f"reports NOMe-seq contexts. Found: {columns}"
            )
        prefixes = ["HCH"]
        for prefix in REPORTED:
            pair = {f"{prefix}n", f"{prefix}b"}
            if pair.intersection(columns):
                if not pair.issubset(columns):
                    raise ValueError(f"{path}: both {prefix}n and {prefix}b are required")
                prefixes.append(prefix)

        rows = []
        chromosomes = set()
        samples = set()
        for values in reader:
            label = f"{path}:{reader.line_num}"
            if len(values) != len(columns):
                raise ValueError(
                    f"{label}: expected {len(columns)} fields, found {len(values)}; "
                    "regenerate malformed statistics with BISCUIT 1.10.2"
                )
            row = dict(zip(columns, values))
            if not row["sample"].strip() or not row["chrm"].strip():
                raise ValueError(f"{label}: sample and chromosome must be non-empty")
            samples.add(row["sample"])
            if len(samples) > 1:
                raise ValueError(f"{label}: expected statistics for one sample")
            if row["chrm"] in chromosomes:
                raise ValueError(f"{label}: duplicate chromosome {row['chrm']!r}")
            chromosomes.add(row["chrm"])
            for prefix in prefixes:
                count_column, beta_column = f"{prefix}n", f"{prefix}b"
                row[count_column], row[beta_column] = parse_count_beta(
                    row[count_column], row[beta_column], f"{label} {prefix}"
                )
            rows.append(row)
    return pd.DataFrame(rows, columns=columns)


def weighted_beta(table: pd.DataFrame, count_column: str, beta_column: str) -> float:
    """Average validated percentages, omitting sites with no observations."""
    observed = table.loc[table[count_column] > 0]
    counts = observed[count_column]
    betas = observed[beta_column]
    total = counts.sum()
    return float((betas * counts).sum() / total) if total else float("nan")


def parse_meth_average(path: Path) -> dict:
    """Reduce one per-chromosome table to the scalars used for QC."""
    table = read_meth_average(path)
    # The aggregate may include sites on chromosomes omitted from the report
    # because they have no HCG coverage. Use it alone, never sum it with rows.
    whole_genome = table.loc[table["chrm"] == "WholeGenome"]
    if not whole_genome.empty:
        table = whole_genome
    record = {"sample": path.name.removesuffix(AVERAGE_SUFFIX)}
    hch_beta = weighted_beta(table, "HCHn", "HCHb")
    record["hch_meth_pct"] = round(hch_beta, 4)
    record["conversion_rate_pct"] = round(100 - hch_beta, 4)
    record["hch_sites"] = int(table["HCHn"].sum())
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
