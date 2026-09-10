"""Batch isolation, reusable reference indexes, and incremental DAG regressions."""

import copy
import csv
import gzip
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from workflow_config import (
    normalize_config,
    normalize_run_config,
    prepare_batch_configs,
    prepare_modality_config,
)


REPO = Path(__file__).resolve().parents[1]
STAR_COMPONENTS = ("Genome", "SA", "SAindex", "genomeParameters.txt")
BISCUIT_SUFFIXES = (
    ".bis.amb", ".bis.ann", ".bis.pac", ".dau.bwt", ".dau.sa", ".par.bwt", ".par.sa"
)


FIXTURE_READ_LENGTH = 50


def write_fastq(path, read_length=FIXTURE_READ_LENGTH):
    """Placeholder reads must be readable: --sjdbOverhang is derived from them."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write(f"@fixture\n{'A' * read_length}\n+\n{'I' * read_length}\n")


class BatchFixtures:
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="scnmt-batches-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        reference_dir = self.root / "reference"
        reference_dir.mkdir()
        self.genome = reference_dir / "genome.fa"
        self.genome.write_text(">chr1\nACGT\n")
        self.genes = reference_dir / "genes.gtf"
        self.genes.touch()
        self.star_index = reference_dir / "star"
        self.star_index.mkdir()
        for filename in STAR_COMPONENTS:
            (self.star_index / filename).touch()
        self.biscuit_index = reference_dir / "biscuit" / "genome"
        self.biscuit_index.parent.mkdir()
        for suffix in BISCUIT_SUFFIXES:
            Path(str(self.biscuit_index) + suffix).touch()
        self.reference = {"genome": str(self.genome), "genes": str(self.genes)}
        self.external_reference = self.reference | {
            "star_index": str(self.star_index),
            "biscuit_index": str(self.biscuit_index),
        }

    def write_batch(self, batch, modalities, rows=(("shared_cell", "same_id"),)):
        for modality in modalities:
            layer = self.root / "data" / batch / modality
            fastq = layer / "fastq"
            fastq.mkdir(parents=True, exist_ok=True)
            self.write_table(layer / "metadata.csv", rows)
            for _, fqid in rows:
                for read in (1, 2):
                    write_fastq(fastq / f"{fqid}_R{read}.fastq.gz")

    def write_table(self, path, rows):
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(("Sample Name", "Unique ID / Lane"))
            writer.writerows(rows)

    def write_facility_table(self, path, rows):
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(("SAMPLE_NAME", "FASTQ_FILE", "READ"))
            writer.writerows(
                (sample, f"{fqid}_R{read}.fastq.gz", read)
                for sample, fqid in rows for read in (1, 2)
            )

    def run_config(self, batches, **overrides):
        return {"reference": copy.deepcopy(self.reference), "batches": batches} | overrides

    def prepared(self, batches, **overrides):
        config = normalize_run_config(self.run_config(batches, **overrides), self.root)
        return prepare_batch_configs(config, self.root)


class BatchConfigTests(BatchFixtures, unittest.TestCase):
    def test_three_mixed_batches_keep_shared_cells_and_fastq_ids_independent(self):
        selections = {"batch_01": "both", "batch_02": "cDNA", "batch_03": "gDNA"}
        for batch, modality in selections.items():
            self.write_batch(batch, ("cDNA", "gDNA") if modality == "both" else (modality,))
        batches = self.prepared({key: {"modality": value} for key, value in selections.items()})
        self.assertEqual(set(batches), set(selections))
        for batch, config in batches.items():
            self.assertEqual(config["dataset"], batch)
            self.assertEqual(config["datadir"], str(self.root / "data" / batch))
            self.assertEqual(config["outdir"], str(self.root / "results" / batch))
            modalities = ("cDNA", "gDNA") if selections[batch] == "both" else (selections[batch],)
            for modality in modalities:
                layer = prepare_modality_config(config, modality)
                self.assertEqual(layer["sample_to_fqid"], {"shared_cell": ["same_id"]})
                self.assertEqual(layer["outdir"], str(self.root / "results" / batch / modality))
                self.assertEqual(
                    layer["fqid_to_reads"]["same_id"]["read1"],
                    str(self.root / "data" / batch / modality / "fastq" / "same_id_R1.fastq.gz"),
                )

    def test_run_helpers_preserve_legacy_single_batch_roots(self):
        raw = {"dataset": "legacy", "modality": "gDNA", "reference": self.reference}
        normalized = normalize_run_config(raw, self.root)
        batches = prepare_batch_configs(normalized, self.root)
        self.assertEqual(set(batches), {None})
        self.assertEqual(batches[None]["dataset"], "legacy")
        self.assertEqual(batches[None]["datadir"], str(self.root / "data"))
        self.assertEqual(batches[None]["outdir"], str(self.root / "results"))

    def test_batch_resequencing_across_tables_retains_runs_once(self):
        self.write_batch("batch_01", ("gDNA",), (("cell", "first"), ("cell", "second")))
        layer_dir = self.root / "data" / "batch_01" / "gDNA"
        for modern in (False, True):
            with self.subTest(modern=modern):
                write_table = self.write_facility_table if modern else self.write_table
                write_table(layer_dir / "metadata.csv", [("cell", "first")])
                write_table(layer_dir / "repeat.csv", [("cell", "second"), ("cell", "first")])
                config = self.prepared({"batch_01": {"modality": "gDNA"}})["batch_01"]
                self.assertEqual(
                    prepare_modality_config(config, "gDNA")["sample_to_fqid"],
                    {"cell": ["first", "second"]},
                )

    def test_batch_overrides_merge_references_and_keep_processing_settings_local(self):
        alternate_genes = self.root / "alternate.gtf"
        alternate_genes.touch()
        for batch in ("batch_01", "batch_02"):
            self.write_batch(batch, ("cDNA",), (("keep", "first"), ("other", "second")))
        config = self.prepared(
            {"batch_01": {"modality": "cDNA", "reference": {"genes": str(alternate_genes)},
                          "samples": ["keep"], "umicount_args": "--UMI_correct"},
             "batch_02": {"modality": "cDNA"}},
            umicount_args="--mm_count_primary",
        )
        self.assertEqual(config["batch_01"]["reference"]["genome"], str(self.genome))
        self.assertEqual(config["batch_01"]["reference"]["genes"], str(alternate_genes))
        self.assertEqual(config["batch_02"]["reference"]["genes"], str(self.genes))
        self.assertEqual(config["batch_01"]["umicount_args"], "--UMI_correct")
        self.assertEqual(config["batch_02"]["umicount_args"], "--mm_count_primary")
        self.assertEqual(
            set(prepare_modality_config(config["batch_01"], "cDNA")["sample_to_fqid"]),
            {"keep"},
        )
        self.assertEqual(
            set(prepare_modality_config(config["batch_02"], "cDNA")["sample_to_fqid"]),
            {"keep", "other"},
        )

    def test_batch_metadata_override_is_independent(self):
        self.write_batch("delivery", ("gDNA",))
        delivery = self.root / "data" / "delivery" / "gDNA"
        config = self.prepared({"batch_01": {
            "modality": "gDNA",
            "ilse_info": {"gDNA": {"metadata": str(delivery / "metadata.csv"),
                                    "fastqdir": str(delivery / "fastq")}},
        }})["batch_01"]
        layer = prepare_modality_config(config, "gDNA")
        self.assertEqual(layer["sample_to_fqid"], {"shared_cell": ["same_id"]})
        self.assertEqual(layer["outdir"], str(self.root / "results" / "batch_01" / "gDNA"))

    def test_batch_keys_and_modality_are_required_and_unambiguous(self):
        invalid = [
            {}, [], {"bad-name": {"modality": "cDNA"}}, {"../bad": {"modality": "gDNA"}},
            {"batch_01": {}}, {"batch_01": {"modality": "RNA"}}, {"batch_01": "gDNA"},
            {"batch_01": {"modality": "cDNA", "dataset": "other"}},
            {"batch_01": {"modality": "cDNA", "datadir": "/other"}},
            {"batch_01": {"modality": "cDNA", "pipeline": "star_umite"}},
        ]
        for batches in invalid:
            with self.subTest(batches=batches):
                with self.assertRaises(ValueError):
                    self.prepared(batches)
        with self.assertRaisesRegex(ValueError, "ilse_info"):
            self.prepared({"batch_01": {"modality": "gDNA"}}, ilse_info={"metadata": "/x.csv"})

    def test_external_index_files_are_passed_to_the_selected_layer(self):
        self.write_batch("batch_01", ("cDNA", "gDNA"))
        config = self.prepared(
            {"batch_01": {"modality": "both"}}, reference=self.external_reference
        )["batch_01"]
        cdna = prepare_modality_config(config, "cDNA")
        gdna = prepare_modality_config(config, "gDNA")
        self.assertEqual(cdna["star_index_dir"], str(self.star_index))
        self.assertEqual(set(cdna["star_index_inputs"]), {
            str(self.star_index / name) for name in STAR_COMPONENTS
        })
        self.assertEqual(gdna["biscuit_index_prefix"], str(self.biscuit_index))
        self.assertEqual(set(gdna["biscuit_index_inputs"]), {
            str(self.biscuit_index) + suffix for suffix in BISCUIT_SUFFIXES
        })

    def test_external_star_index_replaces_cdna_genome_fasta_requirement(self):
        self.write_batch("batch_01", ("cDNA",))
        config = self.prepared(
            {"batch_01": {"modality": "cDNA"}},
            reference={"star_index": str(self.star_index), "genes": str(self.genes)},
        )["batch_01"]
        self.assertEqual(prepare_modality_config(config, "cDNA")["pipeline"], "star_umite")
        with self.assertRaisesRegex(ValueError, "genome"):
            self.prepared(
                {"batch_01": {"modality": "gDNA"}},
                reference={"biscuit_index": str(self.biscuit_index)},
            )

    def test_missing_external_index_components_fail_only_for_selected_modalities(self):
        for key, missing, active, inactive in (
            ("star_index", self.star_index / "SAindex", "cDNA", "gDNA"),
            ("biscuit_index", Path(str(self.biscuit_index) + ".dau.sa"), "gDNA", "cDNA"),
        ):
            with self.subTest(key=key):
                missing.unlink()
                try:
                    raw = {"dataset": "fixture", "reference": self.external_reference}
                    with self.assertRaisesRegex(ValueError, key):
                        normalize_config(raw | {"modality": active}, self.root)
                    with self.assertRaisesRegex(ValueError, key):
                        normalize_config(raw | {"modality": "both"}, self.root)
                    normalize_config(raw | {"modality": inactive}, self.root)
                finally:
                    missing.touch()

    def test_external_index_paths_must_be_absolute(self):
        for key, modality in (("star_index", "cDNA"), ("biscuit_index", "gDNA")):
            with self.subTest(key=key):
                with self.assertRaisesRegex(ValueError, key):
                    normalize_config(
                        {"dataset": "fixture", "modality": modality,
                         "reference": self.reference | {key: "relative/index"}}, self.root
                    )


class BatchWorkflowTests(BatchFixtures, unittest.TestCase):
    def setUp(self):
        super().setUp()
        for filename in ("Snakefile", "defaultconfig.yaml", "workflow_config.py"):
            shutil.copy2(REPO / filename, self.root / filename)
        for directory in ("modules", "envs", "scripts"):
            shutil.copytree(REPO / directory, self.root / directory)

    def run_snakemake(self, config, *arguments):
        (self.root / "snakeconfig.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
        result = subprocess.run(
            [sys.executable, "-m", "snakemake", "--cores", "2", "--printshellcmds", *arguments],
            cwd=self.root, capture_output=True, text=True, timeout=60,
        )
        output = result.stdout + result.stderr
        self.assertEqual(result.returncode, 0, output)
        return output

    def dry_run(self, config):
        return self.run_snakemake(config, "-n")

    def job_blocks(self, output, rule):
        return [block.split("\n\n", 1)[0] for block in output.split(f"rule {rule}:\n")[1:]]

    def complete_fixture(self, config):
        initial_dag = self.dry_run(config)
        # Materialize declared outputs, then record real rule/input/param metadata.
        # Snakemake touch mode invokes no aligners or downstream analysis tools.
        completed_outputs = []
        directories = {"compact_data", "filtered_data", "VMR_matrices"}
        for line in initial_dag.splitlines():
            if line.startswith("    output: "):
                for filename in line.split(": ", 1)[1].split(", "):
                    path = Path(filename)
                    self.assertTrue(path.is_relative_to(self.root / "results"))
                    path.parent.mkdir(parents=True, exist_ok=True)
                    if path.name in directories:
                        path.mkdir(exist_ok=True)
                    else:
                        path.touch()
                    completed_outputs.append(path)
        self.run_snakemake(config, "--touch", "--forceall", "--notemp")
        self.run_snakemake(config, "--delete-temp-output")
        return [path for path in completed_outputs if path.exists()]

    def test_mixed_batches_use_shared_external_indexes_without_rebuilding(self):
        batches = {"batch_01": {"modality": "both"}, "batch_02": {"modality": "cDNA"},
                   "batch_03": {"modality": "gDNA"}}
        for batch, settings in batches.items():
            modes = ("cDNA", "gDNA") if settings["modality"] == "both" else (settings["modality"],)
            self.write_batch(batch, modes)
        output = self.dry_run(self.run_config(batches, reference=self.external_reference))
        for batch, modality in (("batch_01", "cDNA"), ("batch_01", "gDNA"),
                                ("batch_02", "cDNA"), ("batch_03", "gDNA")):
            pipeline = "star_umite" if modality == "cDNA" else "biscuit_methscan"
            self.assertIn(f"results/{batch}/{modality}/{batch}.{pipeline}.h5ad", output)
            self.assertIn(f"data/{batch}/{modality}/fastq/same_id_R1.fastq.gz", output)
            jobs = self.job_blocks(output, f"{batch}_{modality}_align_to_ref")
            self.assertEqual(len(jobs), 1, output)
            components = ([str(self.star_index / name) for name in STAR_COMPONENTS]
                          if modality == "cDNA" else
                          [str(self.biscuit_index) + suffix for suffix in BISCUIT_SUFFIXES])
            input_line = next(line for line in jobs[0].splitlines() if line.startswith("    input:"))
            for component in components:
                self.assertIn(component, input_line)
        self.assertNotIn("prepare_star_indices", output)
        self.assertNotIn("prepare_biscuit_index", output)
        self.assertIn("batch_01_gDNA_prepare_biscuit_reference", output)
        cleanup = self.job_blocks(output, "unload_star_genome")
        self.assertEqual(len(cleanup), 1, output)
        for batch in ("batch_01", "batch_02"):
            self.assertIn(f"results/{batch}/cDNA/star_alignments/shared_cell/", cleanup[0])
            count_job = self.job_blocks(output, f"{batch}_cDNA_count_umis")[0]
            self.assertNotIn("star_cleanup", count_job)
        self.assertIn("results/reference/star_cleanup/", cleanup[0])

    def test_generated_and_external_index_aliases_share_one_cleanup(self):
        for batch in ("batch_01", "batch_02"):
            self.write_batch(batch, ("cDNA",))
        real_results = self.root / "actual_results"
        real_results.mkdir()
        alias_results = self.root / "results_alias"
        alias_results.symlink_to(real_results, target_is_directory=True)
        generated_index = self.genome.parent / "star_index"
        generated_index.mkdir(parents=True)
        for filename in STAR_COMPONENTS:
            (generated_index / filename).touch()
        index_alias = self.root / "star_index_alias"
        index_alias.symlink_to(generated_index, target_is_directory=True)
        config = self.run_config(
            {"batch_01": {"modality": "cDNA"},
             "batch_02": {"modality": "cDNA", "reference": {"star_index": str(index_alias)}}},
            outdir=str(alias_results),
        )
        output = self.dry_run(config)
        cleanup = self.job_blocks(output, "unload_star_genome")
        self.assertEqual(len(cleanup), 1, output)
        for batch in ("batch_01", "batch_02"):
            self.assertIn(f"/{batch}/cDNA/star_alignments/shared_cell/", cleanup[0])
        self.assertIn(str(generated_index), cleanup[0])

    def test_batches_without_external_indexes_keep_index_builders(self):
        self.write_batch("batch_01", ("cDNA", "gDNA"))
        output = self.dry_run(self.run_config({"batch_01": {"modality": "both"}}))
        # The STAR index is generated once beside the genome, not per batch.
        self.assertIn("rule generate_star_index:", output)
        self.assertNotIn("batch_01_cDNA_prepare_star_indices", output)
        self.assertIn(str(self.genome.parent / "star_index" / "Genome"), output)
        self.assertNotIn("results/batch_01/cDNA/reference/star_index", output)
        self.assertIn(f"--sjdbOverhang {FIXTURE_READ_LENGTH - 1}", output)
        self.assertIn("batch_01_gDNA_prepare_biscuit_index", output)
        self.assertIn("results/batch_01/gDNA/reference/biscuit_index", output)

    def test_two_cdna_batches_generate_one_shared_star_index(self):
        for batch in ("batch_01", "batch_02"):
            self.write_batch(batch, ("cDNA",))
        output = self.dry_run(
            self.run_config({batch: {"modality": "cDNA"} for batch in ("batch_01", "batch_02")})
        )
        self.assertEqual(len(self.job_blocks(output, "generate_star_index")), 1, output)
        self.assertEqual(len(self.job_blocks(output, "unload_star_genome")), 1, output)

    def test_longest_sampled_read_sets_the_shared_overhang(self):
        self.write_batch("batch_01", ("cDNA",), rows=(("cell_a", "short_id"),))
        self.write_batch("batch_02", ("cDNA",), rows=(("cell_b", "long_id"),))
        for read in (1, 2):
            write_fastq(
                self.root / "data" / "batch_02" / "cDNA" / "fastq" / f"long_id_R{read}.fastq.gz",
                read_length=FIXTURE_READ_LENGTH + 25,
            )
        output = self.dry_run(
            self.run_config({batch: {"modality": "cDNA"} for batch in ("batch_01", "batch_02")})
        )
        self.assertIn(f"--sjdbOverhang {FIXTURE_READ_LENGTH + 24}", output)

    def test_appending_gdna_batch_leaves_completed_batch_untouched(self):
        self.write_batch("batch_01", ("gDNA",))
        config = self.run_config(
            {"batch_01": {"modality": "gDNA"}}, reference=self.external_reference
        )
        outputs = self.complete_fixture(config)
        stable = self.dry_run(config)
        self.assertNotIn("rule batch_01_gDNA_", stable)
        mtimes = {path: path.stat().st_mtime_ns for path in outputs}
        self.write_batch("batch_02", ("gDNA",))
        config["batches"]["batch_02"] = {"modality": "gDNA"}
        appended = self.dry_run(config)
        self.assertIn("rule batch_02_gDNA_align_to_ref:", appended)
        self.assertIn("rule batch_02_gDNA_construct_methscan_matrix:", appended)
        self.assertNotIn("rule batch_01_gDNA_", appended)
        self.assertEqual(mtimes, {path: path.stat().st_mtime_ns for path in outputs})

    def test_resequencing_updates_one_gdna_cell_and_aggregates_only_in_its_batch(self):
        self.write_batch("batch_01", ("gDNA",), (("cell_a", "first"), ("cell_b", "other")))
        self.write_batch("batch_02", ("gDNA",), (("cell_a", "first"),))
        config = self.run_config(
            {"batch_01": {"modality": "gDNA"}, "batch_02": {"modality": "gDNA"}},
            reference=self.external_reference,
        )
        outputs = self.complete_fixture(config)
        unchanged_batch = [p for p in outputs if p.is_relative_to(self.root / "results" / "batch_02")]
        previous_mtimes = {p: p.stat().st_mtime_ns for p in unchanged_batch}
        layer_dir = self.root / "data" / "batch_01" / "gDNA"
        self.write_facility_table(layer_dir / "resequenced.csv", [("cell_a", "second")])
        for read in (1, 2):
            write_fastq(layer_dir / "fastq" / f"second_R{read}.fastq.gz")
        updated = self.dry_run(config)
        for rule in ("merge_trimmed_reads", "align_to_ref", "deduplicate_and_sort_bam"):
            jobs = self.job_blocks(updated, f"batch_01_gDNA_{rule}")
            self.assertEqual(len(jobs), 1, updated)
            self.assertIn("cell_a", jobs[0])
            self.assertNotIn("cell_b", jobs[0])
        self.assertIn("second_R1.fastq.gz", updated)
        self.assertIn("rule batch_01_gDNA_prepare_methscan_data_and_rename_columns:", updated)
        self.assertIn("rule batch_01_gDNA_construct_methscan_matrix:", updated)
        self.assertIn("rule batch_01_gDNA_build_meth_anndata:", updated)
        self.assertNotIn("rule batch_02_gDNA_", updated)
        self.assertEqual(previous_mtimes, {p: p.stat().st_mtime_ns for p in unchanged_batch})

    def test_changing_strand_mode_recounts_only_the_affected_batch(self):
        for batch in ("batch_01", "batch_02"):
            self.write_batch(batch, ("cDNA",))
        config = self.run_config(
            {"batch_01": {"modality": "cDNA"},
             "batch_02": {"modality": "cDNA", "umicount_args": "--stranded no"}},
            reference=self.external_reference,
        )
        self.complete_fixture(config)
        config["batches"]["batch_01"]["umicount_args"] = "--UMI_correct --stranded reverse"
        changed = self.dry_run(config)
        for rule in ("parse_dump_GTF", "count_umis", "build_trsc_anndata"):
            self.assertIn(f"rule batch_01_cDNA_{rule}:", changed)
            self.assertNotIn(f"rule batch_02_cDNA_{rule}:", changed)
        self.assertIn("umicount_GTF_dump.reverse.pkl", changed)
        for batch in ("batch_01", "batch_02"):
            self.assertNotIn(f"rule {batch}_cDNA_extract_umi:", changed)
            self.assertNotIn(f"rule {batch}_cDNA_align_to_ref:", changed)

    def test_appending_a_batch_preserves_completed_jobs_with_shared_star_cleanup(self):
        self.write_batch("batch_01", ("cDNA",))
        config = self.run_config(
            {"batch_01": {"modality": "cDNA"}}, reference=self.external_reference
        )
        completed_outputs = self.complete_fixture(config)
        stable = self.dry_run(config)
        self.assertNotIn("rule batch_01_cDNA_count_umis:", stable)
        self.assertNotIn("rule batch_01_cDNA_align_to_ref:", stable)
        previous_mtimes = {path: path.stat().st_mtime_ns for path in completed_outputs}
        self.write_batch("batch_02", ("cDNA",))
        config["batches"]["batch_02"] = {"modality": "cDNA"}
        appended = self.dry_run(config)
        self.assertIn("rule batch_02_cDNA_align_to_ref:", appended)
        self.assertIn("rule batch_02_cDNA_count_umis:", appended)
        self.assertNotIn("rule batch_01_cDNA_align_to_ref:", appended)
        self.assertNotIn("rule batch_01_cDNA_count_umis:", appended)
        self.assertNotIn("rule batch_01_cDNA_build_trsc_anndata:", appended)
        self.assertEqual(len(self.job_blocks(appended, "unload_star_genome")), 1)
        self.assertEqual(previous_mtimes, {path: path.stat().st_mtime_ns for path in completed_outputs})


if __name__ == "__main__":
    unittest.main()
