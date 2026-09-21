"""Merge per-cell QC tables into one table MultiQC renders as custom content."""

import argparse
from pathlib import Path

import pandas as pd


# MultiQC reads these commented keys from the top of the file; the remainder is
# an ordinary TSV, so the same file serves as the human-readable QC table.
MULTIQC_HEADER = """\
# id: '{section_id}'
# section_name: '{section_name}'
# description: 'One row per cell, merged from the per-tool QC outputs.'
# plot_type: 'table'
# pconfig:
#     id: '{section_id}_table'
#     namespace: 'scNMT'
#     col1_header: 'Cell'
"""


def merge_tables(paths: list[Path]) -> pd.DataFrame:
    """Outer-join per-tool tables so a missing tool leaves blanks, not gaps."""
    merged = None
    for path in paths:
        table = pd.read_csv(path, sep="\t", comment="#")
        if "sample" not in table.columns:
            raise ValueError(f"{path} has no 'sample' column; found {list(table.columns)}")
        table = table.drop_duplicates(subset="sample")
        merged = table if merged is None else merged.merge(table, on="sample", how="outer")
    if merged is None:
        raise ValueError("No input tables were provided")
    return merged.sort_values("sample")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tables", type=Path, nargs="+", help="Per-tool TSVs")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output TSV")
    parser.add_argument("--section-id", default="per_cell_metrics")
    parser.add_argument("--section-name", default="Per-cell metrics")
    args = parser.parse_args()

    merged = merge_tables(args.tables)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        handle.write(
            MULTIQC_HEADER.format(
                section_id=args.section_id, section_name=args.section_name
            )
        )
        merged.to_csv(handle, sep="\t", index=False)


if __name__ == "__main__":
    main()
