"""Exercise Snakemake's real DAG builder without running the aligners."""

import gzip
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml


REPO = Path(__file__).resolve().parents[1]


FIXTURE_READ_LENGTH = 50


def write_fastq(path, read_length=FIXTURE_READ_LENGTH):
    """Placeholder reads must be readable: --sjdbOverhang is derived from them."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write(f"@fixture\n{'A' * read_length}\n+\n{'I' * read_length}\n")


class WorkflowDryRunTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="scnmt-dag-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for filename in ("Snakefile", "defaultconfig.yaml", "workflow_config.py"):
            shutil.copy2(REPO / filename, self.root / filename)
        for directory in ("modules", "envs", "scripts"):
            shutil.copytree(REPO / directory, self.root / directory)
        references = self.root / "references"
        references.mkdir()
        self.reference = {}
        for key, filename in (
            ("genome", "genome.fa"),
            ("genes", "genes.gtf"),
        ):
            path = references / filename
            path.touch()
            self.reference[key] = str(path)
        for modality, specific_sample in (("cDNA", "rna_only"), ("gDNA", "dna_only")):
            layer = self.root / "data" / modality
            fastq = layer / "fastq"
            fastq.mkdir(parents=True)
            rows = [("shared-cell", "same-id"), (specific_sample, "run1")]
            if modality == "gDNA":
                rows.append((specific_sample, "run2"))
            (layer / "metadata.tsv").write_text(
                "Sample Name\tUnique ID / Lane\n"
                + "".join(f"{sample}\t{fqid}\n" for sample, fqid in rows)
            )
            for _, fqid in rows:
                directory = fastq if modality == "cDNA" else fastq / fqid / "fastq"
                directory.mkdir(parents=True, exist_ok=True)
                for read in (1, 2):
                    write_fastq(directory / f"{fqid}_R{read}.fastq.gz")
            (fastq / "qc.csv").write_text("facility sidecar,not sample metadata\n")

    def dry_run(self, modality, pipeline, reference=None, cli=(), **overrides):
        config = {
            "dataset": "fixture",
            "modality": modality,
            "pipeline": pipeline,
            "reference": self.reference if reference is None else reference,
        }
        config.update(overrides)
        (self.root / "snakeconfig.yaml").write_text(yaml.safe_dump(config))
        completed = subprocess.run(
            [sys.executable, "-m", "snakemake", "-n", "--cores", "2", "--printshellcmds", *cli],
            cwd=self.root, capture_output=True, text=True, timeout=60,
        )
        output = completed.stdout + completed.stderr
        self.assertEqual(completed.returncode, 0, output)
        return output

    def test_all_modality_combinations(self):
        for modality, pipeline in (
            ("cDNA", "star_umite"),
            ("gDNA", "biscuit_methscan"),
            ("both", "star_umite"),
        ):
            with self.subTest(modality=modality, pipeline=pipeline):
                output = self.dry_run(modality, pipeline)
                if modality in ("cDNA", "both"):
                    self.assertIn(f"results/cDNA/fixture.{pipeline}.h5ad", output)
                    self.assertIn("data/cDNA/fastq/same-id_R1.fastq.gz", output)
                    self.assertIn("rna_only", output)
                else:
                    self.assertNotIn("cDNA_", output)
                if modality in ("gDNA", "both"):
                    self.assertIn("results/gDNA/fixture.biscuit_methscan.h5ad", output)
                    self.assertIn("data/gDNA/fastq/same-id/fastq/same-id_R1.fastq.gz", output)
                    self.assertIn("dna_only", output)
                    self.assertIn("gDNA_merge_trimmed_reads", output)
                    self.assertIn("results/gDNA/qc_plots/cell_stats.png", output)
                    self.assertIn("rule gDNA_plot_methscan_cell_stats:", output)
                    plot_job = output.split("rule gDNA_plot_methscan_cell_stats:", 1)[1]
                    plot_job = plot_job.split("\n\n", 1)[0]
                    self.assertIn("methscan/compact_data/cell_stats.csv", plot_job)
                    self.assertNotIn("methscan/filtered_data", plot_job)
                else:
                    self.assertNotIn("gDNA_", output)
                    self.assertNotIn("cell_stats.png", output)
                # Every declared generated file and log must live below results/,
                # except the shared STAR index, which is generated beside the genome.
                index_dir = Path(self.reference["genome"]).parent / "star_index"
                for line in output.splitlines():
                    if line.startswith(("    output: ", "    log: ")):
                        paths = line.split(": ", 1)[1].split(", ")
                        for path in paths:
                            self.assertTrue(
                                Path(path).is_relative_to(self.root / "results")
                                or Path(path).is_relative_to(index_dir),
                                line,
                            )

    def test_modern_facility_metadata_builds_gdna_and_both_dags(self):
        metadata = self.root / "data" / "gDNA" / "metadata.tsv"
        rows = [("shared-cell", "same-id"), ("dna_only", "run1"), ("dna_only", "run2")]
        metadata.write_text(
            "SAMPLE_NAME\tFASTQ_FILE\tREAD\n"
            + "".join(
                f"{sample}\t{fqid}_R{read}.fastq.gz\t{read}\n"
                for sample, fqid in rows for read in (1, 2)
            )
            + "\tUndetermined_S0_L001_R1_001.fastq.gz\t1\n"
        )
        original = metadata.read_bytes()
        for modality in ("gDNA", "both"):
            with self.subTest(modality=modality):
                output = self.dry_run(modality, "star_umite")
                self.assertIn("results/gDNA/fixture.biscuit_methscan.h5ad", output)
                self.assertIn("gDNA_merge_trimmed_reads", output)
                self.assertIn("data/gDNA/fastq/run2/fastq/run2_R2.fastq.gz", output)
                self.assertNotIn("Undetermined", output)
                if modality == "both":
                    self.assertIn("results/cDNA/fixture.star_umite.h5ad", output)
                self.assertEqual(metadata.read_bytes(), original)

    def methscan_command(self, output, step):
        commands = []
        for line in output.splitlines():
            command = line.removeprefix("Shell command: ").strip()
            if command.startswith(f"methscan {step} "):
                commands.append(shlex.split(command))
        self.assertEqual(len(commands), 1, output)
        return commands[0]

    def test_methscan_commands_follow_configured_defaults(self):
        output = self.dry_run("gDNA", "biscuit_methscan")
        defaults = yaml.safe_load((self.root / "defaultconfig.yaml").read_text())
        directory = self.root / "results" / "gDNA" / "methscan"
        filtered = str(directory / "filtered_data")
        vmrs = str(directory / "fixture_VMRs.bed")
        self.assertEqual(
            self.methscan_command(output, "filter"),
            [
                "methscan", "filter", *shlex.split(defaults["methscan_filter_args"]),
                str(directory / "compact_data"), filtered,
                "2>", str(directory / "methscan_filter.log"),
            ],
        )
        self.assertEqual(
            self.methscan_command(output, "smooth"),
            [
                "methscan", "smooth", *shlex.split(defaults["methscan_smooth_args"]),
                filtered, "2>>", str(directory / "methscan_scan.log"),
            ],
        )
        self.assertEqual(
            self.methscan_command(output, "scan"),
            [
                "methscan", "scan", *shlex.split(defaults["methscan_scan_args"]),
                "--threads", "2", filtered, vmrs,
                "2>>", str(directory / "methscan_scan.log"),
            ],
        )
        self.assertEqual(
            self.methscan_command(output, "matrix"),
            [
                "methscan", "matrix", *shlex.split(defaults["methscan_matrix_args"]),
                "--threads", "2", vmrs, filtered,
                str(directory / "VMR_matrices"), "2>",
                str(directory / "methscan_matrix.log"),
            ],
        )

    def test_methscan_arguments_reach_their_commands_in_both_mode(self):
        arguments = {
            "filter": "--min-sites 250 --min-meth 5 --max-meth 95",
            "smooth": "--bandwidth 750 --use-weights",
            "scan": "--min-cells 3 --stepsize 50",
            "matrix": "--threads 1",
        }
        output = self.dry_run(
            "both", "star_umite",
            **{f"methscan_{step}_args": value for step, value in arguments.items()},
        )
        for step, value in arguments.items():
            with self.subTest(step=step):
                command = self.methscan_command(output, step)
                expected_args = shlex.split(value)
                self.assertTrue(
                    any(
                        command[index:index + len(expected_args)] == expected_args
                        for index in range(len(command))
                    ),
                    command,
                )
        matrix = self.methscan_command(output, "matrix")
        thread_values = [
            matrix[index + 1] for index, token in enumerate(matrix) if token == "--threads"
        ]
        self.assertEqual(thread_values, ["1", "2"])
        self.assertNotIn("10000", self.methscan_command(output, "filter"))
        self.assertIn("cDNA_align_to_ref", output)

    def test_empty_filter_arguments_remove_default_thresholds(self):
        output = self.dry_run("gDNA", "biscuit_methscan", methscan_filter_args="")
        command = self.methscan_command(output, "filter")
        self.assertNotIn("--min-sites", command)
        self.assertNotIn("--min-meth", command)
        self.assertNotIn("--max-meth", command)
        self.assertTrue(command[2].endswith("/methscan/compact_data"))

    def test_filter_cell_list_is_a_dependency_without_extra_positional_arguments(self):
        cell_names = self.root / "selected_cells.txt"
        cell_names.write_text("shared-cell\n", encoding="utf-8")
        directory = self.root / "results" / "gDNA" / "methscan"
        for separator in (" ", "="):
            with self.subTest(separator=separator):
                arguments = f"--cell-names{separator}{cell_names} --keep"
                output = self.dry_run(
                    "gDNA", "biscuit_methscan", methscan_filter_args=arguments
                )
                filter_job = output.split("rule gDNA_filter_methscan_data:", 1)[1]
                input_line = next(
                    line for line in filter_job.splitlines() if line.startswith("    input: ")
                )
                self.assertIn(str(cell_names), input_line)
                self.assertIn(str(directory / "compact_data"), input_line)
                self.assertEqual(
                    self.methscan_command(output, "filter"),
                    ["methscan", "filter", *shlex.split(arguments),
                     str(directory / "compact_data"), str(directory / "filtered_data"),
                     "2>", str(directory / "methscan_filter.log")],
                )

    def test_gdna_only_needs_genome_and_its_own_inputs(self):
        shutil.rmtree(self.root / "data" / "cDNA")
        output = self.dry_run("gDNA", "star_umite", {"genome": self.reference["genome"]})
        self.assertIn("gDNA_align_to_ref", output)
        self.assertNotIn("cDNA_", output)

    def test_cdna_only_does_not_need_gdna_inputs(self):
        shutil.rmtree(self.root / "data" / "gDNA")
        output = self.dry_run("cDNA", "star_umite")
        self.assertIn("cDNA_align_to_ref", output)
        self.assertNotIn("gDNA_", output)

    def test_umicount_strand_mode_reaches_the_gtf_dump(self):
        # umicount refuses a dump parsed under a different --stranded mode, so
        # the mode has to name the dump and be passed to both rules.
        shutil.rmtree(self.root / "data" / "gDNA")
        for mode, arguments in (("no", "--UMI_correct"),
                                ("reverse", "--UMI_correct --stranded reverse")):
            with self.subTest(mode=mode):
                output = self.dry_run("cDNA", "star_umite", umicount_args=arguments)
                dump = f"results/cDNA/reference/umicount_GTF_dump.{mode}.pkl"
                self.assertIn(f"--GTF_dump {self.root / dump} --stranded {mode}", output)
                self.assertIn(f"--GTF_skip_parse {self.root / dump}", output)
                if mode != "no":
                    self.assertIn(f"--stranded {mode} \\\n            --combine_unspliced", output)

    def test_cli_modality_overrides_config_file(self):
        output = self.dry_run("cDNA", "star_umite", cli=("--config", "modality=both"))
        self.assertIn("cDNA_align_to_ref", output)
        self.assertIn("gDNA_align_to_ref", output)


if __name__ == "__main__":
    unittest.main()
