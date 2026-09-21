"""Collect per-cell mapping statistics from STAR Log.final.out files."""

import argparse
from pathlib import Path

import pandas as pd


# STAR writes "  label |\tvalue"; only these rows are worth tracking per cell.
STAR_FIELDS = {
    "Number of input reads": "input_reads",
    "Average input read length": "avg_input_read_length",
    "Uniquely mapped reads number": "uniquely_mapped",
    "Uniquely mapped reads %": "uniquely_mapped_pct",
    "Average mapped length": "avg_mapped_length",
    "Mismatch rate per base, %": "mismatch_rate_pct",
    "% of reads mapped to multiple loci": "multimapped_pct",
    "% of reads mapped to too many loci": "too_many_loci_pct",
    "% of reads unmapped: too short": "unmapped_too_short_pct",
    "% of reads unmapped: other": "unmapped_other_pct",
    "% of chimeric reads": "chimeric_pct",
}
LOG_SUFFIX = "_Log.final.out"


def parse_star_log(path: Path) -> dict:
    """Read one STAR summary into a flat record keyed by cell name."""
    record = {"sample": path.name.removesuffix(LOG_SUFFIX)}
    for line in path.read_text().splitlines():
        if "|" not in line:
            continue
        label, _, value = line.partition("|")
        column = STAR_FIELDS.get(label.strip())
        if column is None:
            continue
        record[column] = pd.to_numeric(value.strip().rstrip("%"), errors="coerce")
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", type=Path, nargs="+", help="STAR *_Log.final.out files")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output TSV")
    args = parser.parse_args()

    table = pd.DataFrame([parse_star_log(path) for path in args.logs])
    table = table.reindex(columns=["sample", *STAR_FIELDS.values()])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.sort_values("sample").to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
