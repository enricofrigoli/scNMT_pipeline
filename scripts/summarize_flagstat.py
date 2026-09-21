"""Collect per-cell alignment counts from samtools flagstat output."""

import argparse
import re
from pathlib import Path

import pandas as pd


# "N + M <label> (P% : Q%)" with the QC-failed count and percentages optional.
FLAGSTAT_PATTERN = re.compile(
    r"^(?P<passed>\d+)\s*\+\s*(?P<failed>\d+)\s+(?P<label>[^(]+?)\s*(?:\(.*\))?$"
)
FLAGSTAT_FIELDS = {
    "in total": "total_reads",
    "primary": "primary_reads",
    "secondary": "secondary_reads",
    "supplementary": "supplementary_reads",
    "duplicates": "duplicate_reads",
    "primary duplicates": "primary_duplicate_reads",
    "mapped": "mapped_reads",
    "primary mapped": "primary_mapped_reads",
    "properly paired": "properly_paired_reads",
    "singletons": "singleton_reads",
}
FLAGSTAT_SUFFIX = ".flagstat"


def parse_flagstat(path: Path) -> dict:
    """Read one flagstat report and derive the mapping percentage."""
    record = {"sample": path.name.removesuffix(FLAGSTAT_SUFFIX)}
    for line in path.read_text().splitlines():
        match = FLAGSTAT_PATTERN.match(line.strip())
        if match is None:
            continue
        # "in total (QC-passed reads + QC-failed reads)" keeps its parenthetical.
        label = match.group("label").split("(")[0].strip()
        column = FLAGSTAT_FIELDS.get(label)
        if column is not None:
            record[column] = int(match.group("passed"))
    total = record.get("primary_reads") or record.get("total_reads")
    mapped = record.get("primary_mapped_reads", record.get("mapped_reads"))
    record["mapped_pct"] = (
        round(100 * mapped / total, 4) if total and mapped is not None else float("nan")
    )
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", type=Path, nargs="+", help="samtools flagstat files")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output TSV")
    args = parser.parse_args()

    table = pd.DataFrame([parse_flagstat(path) for path in args.reports])
    table = table.reindex(columns=["sample", "mapped_pct", *FLAGSTAT_FIELDS.values()])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.sort_values("sample").to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
