"""Regression checks against the BISCUIT NOMe methylation statistics format."""

import csv
import importlib.util
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/summarize_hch_conversion.py"
SPEC = importlib.util.spec_from_file_location("summarize_hch_conversion", SCRIPT)
CONVERSION = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CONVERSION)

BISCUIT_COLUMNS = (
    "sample", "chrm", "HCGn", "HCGb", "HCHGn", "HCHGb",
    "HCHHn", "HCHHb", "HCHn", "HCHb", "GCn", "GCb",
)
OUTPUT_COLUMNS = (
    "sample", "conversion_rate_pct", "hch_meth_pct", "hch_sites",
    "cpg_meth_pct", "gpc_meth_pct",
)


def chromosome_row(**overrides):
    row = {
        "sample": "cell.dedup_sorted",
        "chrm": "chr1",
        "HCGn": "100",
        "HCGb": "70.000%",
        "HCHGn": "25",
        "HCHGb": "2.000%",
        "HCHHn": "75",
        "HCHHb": "2.000%",
        "HCHn": "100",
        "HCHb": "2.000%",
        "GCn": "100",
        "GCb": "60.000%",
    }
    row.update(overrides)
    return row


class HchConversionTests(unittest.TestCase):
    def setUp(self):
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        self.root = Path(temporary_directory.name)

    def write_table(self, rows, *, sample="cell", columns=BISCUIT_COLUMNS):
        path = self.root / f"{sample}_meth_average.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(columns)
            for row in rows:
                writer.writerow([row[column] for column in columns])
        return path

    def test_literal_percentage_suffixes_and_filename_sample_name(self):
        path = self.write_table([chromosome_row()], sample="plate4_A1")

        record = CONVERSION.parse_meth_average(path)

        self.assertEqual(record, {
            "sample": "plate4_A1",
            "conversion_rate_pct": 98.0,
            "hch_meth_pct": 2.0,
            "hch_sites": 100,
            "cpg_meth_pct": 70.0,
            "gpc_meth_pct": 60.0,
        })

    def test_chromosomes_are_weighted_by_observed_sites(self):
        path = self.write_table([
            chromosome_row(),
            chromosome_row(chrm="chr2", HCHn="300", HCHb="4.000%"),
        ])

        record = CONVERSION.parse_meth_average(path)

        self.assertEqual(record["hch_sites"], 400)
        self.assertEqual(record["hch_meth_pct"], 3.5)
        self.assertEqual(record["conversion_rate_pct"], 96.5)

    def test_whole_genome_is_authoritative_and_not_counted_twice(self):
        path = self.write_table([
            chromosome_row(),
            chromosome_row(chrm="chr2", HCHn="300", HCHb="4.000%"),
            # BISCUIT can omit a chromosome row when that chromosome has no HCG.
            chromosome_row(
                chrm="WholeGenome", HCHn="500", HCHb="5.000%",
                HCGn="200", HCGb="75.000%", GCn="200", GCb="55.000%",
            ),
        ])

        record = CONVERSION.parse_meth_average(path)

        self.assertEqual(record["hch_sites"], 500)
        self.assertEqual(record["hch_meth_pct"], 5.0)
        self.assertEqual(record["conversion_rate_pct"], 95.0)
        self.assertEqual(record["cpg_meth_pct"], 75.0)
        self.assertEqual(record["gpc_meth_pct"], 55.0)

    def test_rejects_duplicate_whole_genome_rows(self):
        path = self.write_table([
            chromosome_row(chrm="WholeGenome"),
            chromosome_row(chrm="WholeGenome"),
        ])
        with self.assertRaises(ValueError):
            CONVERSION.parse_meth_average(path)

    def test_rejects_duplicate_chromosome_rows_even_with_whole_genome(self):
        path = self.write_table([
            chromosome_row(), chromosome_row(),
            chromosome_row(chrm="WholeGenome"),
        ])
        with self.assertRaises(ValueError):
            CONVERSION.parse_meth_average(path)

    def test_rejects_multiple_samples_in_one_table(self):
        path = self.write_table([
            chromosome_row(),
            chromosome_row(chrm="chr2", sample="another_cell.dedup_sorted"),
        ])
        with self.assertRaises(ValueError):
            CONVERSION.parse_meth_average(path)

    def test_rejects_missing_hch_columns(self):
        for missing in ("HCHn", "HCHb"):
            with self.subTest(missing=missing):
                columns = tuple(c for c in BISCUIT_COLUMNS if c != missing)
                path = self.write_table([chromosome_row()], columns=columns)
                with self.assertRaises(ValueError):
                    CONVERSION.parse_meth_average(path)

    def test_rejects_malformed_row_lengths_instead_of_inferring_an_index(self):
        row = chromosome_row()
        fields = [row[column] for column in BISCUIT_COLUMNS]
        for malformed in (fields[:-1], fields + ["2.000%"] * 4):
            with self.subTest(field_count=len(malformed)):
                path = self.root / "malformed_meth_average.tsv"
                path.write_text(
                    "\t".join(BISCUIT_COLUMNS) + "\n"
                    + "\t".join(malformed) + "\n",
                    encoding="utf-8",
                )
                with self.assertRaises(ValueError) as error:
                    CONVERSION.parse_meth_average(path)
                self.assertIn(path.name, str(error.exception))

    def test_rejects_invalid_counts_for_each_reported_context(self):
        for column in ("HCHn", "HCGn", "GCn"):
            for value in ("-1", "1.5", "nan", "inf", "wrong", "", "."):
                with self.subTest(column=column, value=value):
                    path = self.write_table([chromosome_row(**{column: value})])
                    with self.assertRaises(ValueError):
                        CONVERSION.parse_meth_average(path)

    def test_rejects_invalid_percentages_for_each_reported_context(self):
        for column in ("HCHb", "HCGb", "GCb"):
            for value in ("-0.1%", "100.001%", "nan", "inf%", "wrong"):
                with self.subTest(column=column, value=value):
                    path = self.write_table([chromosome_row(**{column: value})])
                    with self.assertRaises(ValueError):
                        CONVERSION.parse_meth_average(path)

    def test_rejects_missing_beta_with_positive_count(self):
        for column in ("HCHb", "HCGb", "GCb"):
            for value in (".", ""):
                with self.subTest(column=column, value=value):
                    path = self.write_table([chromosome_row(**{column: value})])
                    with self.assertRaises(ValueError):
                        CONVERSION.parse_meth_average(path)

    def test_zero_coverage_allows_missing_betas_and_reports_nan(self):
        for missing_beta in (".", ""):
            with self.subTest(missing_beta=missing_beta):
                row = chromosome_row()
                for column in BISCUIT_COLUMNS[2:]:
                    row[column] = "0" if column.endswith("n") else missing_beta
                record = CONVERSION.parse_meth_average(self.write_table([row]))

                self.assertEqual(record["hch_sites"], 0)
                for column in OUTPUT_COLUMNS[1:]:
                    if column != "hch_sites":
                        self.assertTrue(math.isnan(record[column]), column)

    def test_zero_coverage_row_does_not_change_weighted_beta(self):
        path = self.write_table([
            chromosome_row(chrm="chr1", HCHn="0", HCHb="."),
            chromosome_row(chrm="chr2", HCHn="100", HCHb="4.000%"),
        ])

        record = CONVERSION.parse_meth_average(path)

        self.assertEqual(record["hch_sites"], 100)
        self.assertEqual(record["conversion_rate_pct"], 96.0)

    def test_header_only_table_reports_nan_instead_of_perfect_conversion(self):
        record = CONVERSION.parse_meth_average(self.write_table([]))

        self.assertEqual(record["hch_sites"], 0)
        for column in OUTPUT_COLUMNS[1:]:
            if column != "hch_sites":
                self.assertTrue(math.isnan(record[column]), column)

    def test_numeric_percentages_without_suffix_remain_supported(self):
        path = self.write_table([
            chromosome_row(HCHb="2.0", HCGb="70.0", GCb="60.0"),
        ])

        record = CONVERSION.parse_meth_average(path)

        self.assertEqual(record["conversion_rate_pct"], 98.0)
        self.assertEqual(record["cpg_meth_pct"], 70.0)
        self.assertEqual(record["gpc_meth_pct"], 60.0)

    def test_cli_keeps_output_columns_and_sorts_cell_names(self):
        cell_b = self.write_table([chromosome_row()], sample="cell_B")
        cell_a = self.write_table(
            [chromosome_row(HCHb="4.000%")], sample="cell_A",
        )
        output = self.root / "qc/tables/conversion_metrics.tsv"

        result = subprocess.run(
            [sys.executable, str(SCRIPT), str(cell_b), str(cell_a), "-o", str(output)],
            capture_output=True, text=True, check=False,
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        with output.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(reader.fieldnames, list(OUTPUT_COLUMNS))
            rows = list(reader)
        self.assertEqual([row["sample"] for row in rows], ["cell_A", "cell_B"])
        self.assertEqual([float(row["conversion_rate_pct"]) for row in rows], [96, 98])
        self.assertEqual([int(row["hch_sites"]) for row in rows], [100, 100])


if __name__ == "__main__":
    unittest.main()
