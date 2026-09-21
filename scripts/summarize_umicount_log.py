"""Extract the per-cell counting summary UMITE writes into its log."""

import argparse
import re
from pathlib import Path

import pandas as pd


# umicount logs one INFO line per BAM. The rest of the file is a per-read
# warning stream, so the log is streamed and only these lines are kept.
SUMMARY_PATTERN = re.compile(
    r"\[INFO\]\s+(?P<sample>\S+):\s+"
    r"(?P<total_reads>\d+) reads,\s+"
    r"(?P<uncounted>\d+) uncounted reads \((?P<uncounted_pct>[\d.]+)%\),\s+"
    r"(?P<umi_reads>\d+) U-reads \((?P<umi_reads_pct>[\d.]+)%\),\s+"
    r"(?P<internal_reads>\d+) R-reads \((?P<internal_reads_pct>[\d.]+)%\),\s+"
    r"(?P<duplicate_reads>\d+) D-reads \((?P<duplicate_reads_pct>[\d.]+)%\),\s+"
    r"(?P<corrected_umis>\d+) counts in corrected UMIs"
)
NUMERIC_COLUMNS = [
    "total_reads", "uncounted", "uncounted_pct", "umi_reads", "umi_reads_pct",
    "internal_reads", "internal_reads_pct", "duplicate_reads",
    "duplicate_reads_pct", "corrected_umis",
]


def parse_umicount_log(path: Path, suffix: str) -> pd.DataFrame:
    """Stream the log so a multi-hundred-megabyte file is never held in memory."""
    records = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if "[INFO]" not in line:
                continue
            match = SUMMARY_PATTERN.search(line)
            if match is None:
                continue
            record = match.groupdict()
            record["sample"] = record["sample"].removesuffix(suffix)
            records.append(record)
    if not records:
        raise ValueError(f"No per-cell summary lines found in {path}")
    table = pd.DataFrame(records)
    table[NUMERIC_COLUMNS] = table[NUMERIC_COLUMNS].apply(pd.to_numeric)
    return table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path, help="umicount.log")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output TSV")
    parser.add_argument(
        "--samplename-suffix", default="_Aligned.qn_sorted.bam",
        help="BAM suffix stripped to recover the cell name",
    )
    args = parser.parse_args()

    table = parse_umicount_log(args.log, args.samplename_suffix)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.sort_values("sample").to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
