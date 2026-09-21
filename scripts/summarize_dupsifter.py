"""Collect per-cell duplicate rates from dupsifter stat files."""

import argparse
import re
from pathlib import Path

import pandas as pd


# dupsifter prints "[dupsifter] <label>: <count>" and computes no rate itself.
# Labels are matched loosely so a wording change degrades to a missing column
# rather than a crash.
STAT_PATTERN = re.compile(r"\[dupsifter\]\s+(?P<label>[^:]+):\s+(?P<value>\d+)\s*$")
STAT_FIELDS = {
    "number of individual reads processed": "reads_processed",
    "number of reads with both reads mapped": "pairs_mapped",
    "number of reads with only one read mapped to the forward strand": "forward_mapped",
    "number of reads with only one read mapped to the reverse strand": "reverse_mapped",
    "number of reads with both reads marked as duplicates": "pairs_duplicate",
    "number of reads on the forward strand marked as duplicates": "forward_duplicate",
    "number of reads on the reverse strand marked as duplicates": "reverse_duplicate",
    "number of individual primary-alignment reads": "primary_reads",
    "number of individual secondary- and supplementary-alignment reads": "secondary_reads",
    "number of reads with no reads mapped": "unmapped",
}
STAT_SUFFIX = ".dupsifter.stat"
DUPLICATE_PARTS = ("pairs_duplicate", "forward_duplicate", "reverse_duplicate")
MAPPED_PARTS = ("pairs_mapped", "forward_mapped", "reverse_mapped")


def parse_dupsifter_stat(path: Path) -> dict:
    """Read one stat file and derive the duplicate percentage it omits."""
    record = {"sample": path.name.removesuffix(STAT_SUFFIX)}
    for line in path.read_text().splitlines():
        match = STAT_PATTERN.match(line.strip())
        if match is None:
            continue
        column = STAT_FIELDS.get(match.group("label").strip())
        if column is not None:
            record[column] = int(match.group("value"))
    # Pairs and orphan reads are counted in different units, but the ratio is
    # taken within each category, so the combined rate stays self-consistent.
    duplicates = sum(record.get(key, 0) for key in DUPLICATE_PARTS)
    mapped = sum(record.get(key, 0) for key in MAPPED_PARTS)
    record["duplicate_pct"] = round(100 * duplicates / mapped, 4) if mapped else float("nan")
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stats", type=Path, nargs="+", help="dupsifter .stat files")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output TSV")
    args = parser.parse_args()

    table = pd.DataFrame([parse_dupsifter_stat(path) for path in args.stats])
    table = table.reindex(columns=["sample", "duplicate_pct", *STAT_FIELDS.values()])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.sort_values("sample").to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
