"""Regression checks for independent, already demultiplexed input modalities."""

import copy
import csv
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from workflow_config import normalize_config, prepare_modality_config, read_metadata


class WorkflowConfigTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.genome = self.root / "genome.fa"
        self.genome.write_text(">chr1\nACGT\n", encoding="utf-8")
        self.genes = self.root / "genes.gtf"
        self.genes.touch()
        self.config = {
            "dataset": "test_dataset",
            "reference": {"genome": str(self.genome), "genes": str(self.genes)},
        }

    def write_metadata(
        self, path, rows, *, delimiter=",", columns=("Sample Name", "Unique ID / Lane")
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter=delimiter)
            writer.writerow(columns)
            writer.writerows(rows)
        return path

    def write_pair(self, directory, fqid):
        directory.mkdir(parents=True, exist_ok=True)
        paths = {}
        for read in (1, 2):
            path = directory / f"{fqid}_R{read}.fastq.gz"
            path.touch()
            paths[f"read{read}"] = str(path)
        return paths

    def create_layer(self, modality, rows):
        directory = self.root / "data" / modality
        metadata = self.write_metadata(directory / "metadata.csv", rows)
        fastq = directory / "fastq"
        fastq.mkdir(exist_ok=True)
        for _, fqid in rows:
            if fqid:
                self.write_pair(fastq, fqid)
        return metadata, fastq

    def write_modern_metadata(self, path, rows):
        return self.write_metadata(
            path, rows, columns=("SAMPLE_NAME", "FASTQ_FILE", "READ")
        )

    def create_modern_layer(self, modality, rows):
        metadata, fastq = self.create_layer(modality, rows)
        self.write_modern_metadata(
            metadata,
            [(sample, f"{fqid}_R{read}.fastq.gz", read)
             for sample, fqid in rows for read in (1, 2)],
        )
        return metadata, fastq

    def normalized(self, **overrides):
        config = copy.deepcopy(self.config)
        config.update(overrides)
        return normalize_config(config, self.root)

    def test_modern_metadata_mate_rows_resolve_to_one_read_pair(self):
        metadata, fastq = self.create_modern_layer("gDNA", [("cell_A", "facility_001")])
        original = metadata.read_bytes()
        layer = prepare_modality_config(self.normalized(modality="gDNA"), "gDNA")
        self.assertEqual(layer["sample_to_fqid"], {"cell_A": ["facility_001"]})
        self.assertEqual(
            layer["fqid_to_reads"],
            {"facility_001": {
                "read1": str(fastq / "facility_001_R1.fastq.gz"),
                "read2": str(fastq / "facility_001_R2.fastq.gz"),
            }},
        )
        self.assertEqual(metadata.read_bytes(), original)

    def test_modern_metadata_keeps_resequenced_stems_and_deduplicates_mate_rows(self):
        self.create_modern_layer(
            "gDNA", [("cell_A", "run_1"), ("cell_A", "run_2"), ("cell_A", "run_1")]
        )
        layer = prepare_modality_config(self.normalized(modality="gDNA"), "gDNA")
        self.assertEqual(layer["sample_to_fqid"], {"cell_A": ["run_1", "run_2"]})
        self.assertEqual(set(layer["fqid_to_reads"]), {"run_1", "run_2"})

    def test_modern_metadata_rejects_mates_assigned_to_different_samples(self):
        metadata, _ = self.create_modern_layer("gDNA", [("cell_A", "id")])
        self.write_modern_metadata(
            metadata, [("cell_A", "id_R1.fastq.gz", 1), ("cell_B", "id_R2.fastq.gz", 2)]
        )
        with self.assertRaisesRegex(ValueError, "assigned to both"):
            prepare_modality_config(self.normalized(modality="gDNA"), "gDNA")

    def test_modern_metadata_rejects_read_column_mismatch(self):
        metadata = self.root / "modern.csv"
        for read in (2, 3, "R1", ""):
            with self.subTest(read=read):
                self.write_modern_metadata(metadata, [("cell", "id_R1.fastq.gz", read)])
                with self.assertRaisesRegex(ValueError, "READ"):
                    read_metadata(str(metadata))

    def test_modern_metadata_rejects_malformed_fastq_filename(self):
        metadata = self.root / "modern.csv"
        for filename in ("id_R1_001.fastq.gz", "id.fastq.gz", "id_R1.fastq", "/id_R1.fastq.gz"):
            with self.subTest(filename=filename):
                self.write_modern_metadata(metadata, [("cell", filename, 1)])
                with self.assertRaisesRegex(ValueError, "FASTQ_FILE"):
                    read_metadata(str(metadata))

    def test_modern_metadata_requires_all_three_columns(self):
        columns = ("SAMPLE_NAME", "FASTQ_FILE", "READ")
        metadata = self.root / "modern.csv"
        for missing in columns:
            with self.subTest(missing=missing):
                self.write_metadata(metadata, [], columns=[c for c in columns if c != missing])
                with self.assertRaisesRegex(ValueError, missing):
                    read_metadata(str(metadata))

    def test_modern_metadata_skips_undetermined_reads_with_blank_sample(self):
        metadata, _ = self.create_modern_layer("gDNA", [("cell", "id")])
        self.write_modern_metadata(
            metadata,
            [("cell", "id_R1.fastq.gz", 1), ("cell", "id_R2.fastq.gz", 2),
             ("", "Undetermined_S0_L001_R1_001.fastq.gz", 1),
             ("", "Undetermined_S0_L001_R2_001.fastq.gz", 2)],
        )
        layer = prepare_modality_config(self.normalized(modality="gDNA"), "gDNA")
        self.assertEqual(layer["sample_to_fqid"], {"cell": ["id"]})

    def test_modern_metadata_rejects_blank_sample_for_named_fastq(self):
        metadata = self.write_modern_metadata(
            self.root / "modern.csv", [("", "id_R1.fastq.gz", 1)]
        )
        with self.assertRaisesRegex(ValueError, "SAMPLE_NAME|sample name"):
            read_metadata(str(metadata))

    def test_modern_metadata_preserves_sample_name_whitespace(self):
        sample = " sample with spaces "
        metadata, _ = self.create_modern_layer("gDNA", [(sample, "id")])
        table = read_metadata(str(metadata))
        self.assertEqual(table["Sample Name"].tolist(), [sample, sample])
        with self.assertRaisesRegex(ValueError, "sample name"):
            prepare_modality_config(self.normalized(modality="gDNA"), "gDNA")

    def test_literal_na_sample_names_survive_both_metadata_schemas(self):
        for modern in (False, True):
            with self.subTest(modern=modern):
                create_layer = self.create_modern_layer if modern else self.create_layer
                create_layer("gDNA", [("NA", "NA")])
                layer = prepare_modality_config(self.normalized(modality="gDNA"), "gDNA")
                self.assertEqual(layer["sample_to_fqid"], {"NA": ["NA"]})

    def test_legacy_metadata_still_skips_blank_fastq_ids(self):
        self.create_layer("gDNA", [("cell", "id"), ("unused", ""), ("", "")])
        layer = prepare_modality_config(self.normalized(modality="gDNA"), "gDNA")
        self.assertEqual(layer["sample_to_fqid"], {"cell": ["id"]})

    def test_shared_sample_and_fastq_ids_stay_in_separate_layers(self):
        for modality in ("cDNA", "gDNA"):
            self.create_layer(modality, [("cell_1", "facility_1")])
        config = self.normalized(modality="both")
        layers = {
            modality: prepare_modality_config(config, modality)
            for modality in ("cDNA", "gDNA")
        }
        self.assertEqual(config["datadir"], str(self.root / "data"))
        self.assertEqual(config["outdir"], str(self.root / "results"))
        for modality, layer in layers.items():
            self.assertEqual(layer["sample_to_fqid"], {"cell_1": ["facility_1"]})
            self.assertEqual(layer["outdir"], str(self.root / "results" / modality))
            self.assertEqual(
                layer["fqid_to_reads"]["facility_1"]["read1"],
                str(self.root / "data" / modality / "fastq" / "facility_1_R1.fastq.gz"),
            )
        self.assertEqual(layers["cDNA"]["pipeline"], "star_umite")
        self.assertEqual(layers["gDNA"]["pipeline"], "biscuit_methscan")
        self.assertNotIn("sample_to_fqid", config)
        self.assertEqual(config["modality"], "both")

    def test_both_accepts_different_cell_sets_without_matching(self):
        self.create_layer("cDNA", [("rna_cell", "rna_id")])
        self.create_layer("gDNA", [("dna_cell", "dna_id"), ("other", "other_id")])
        config = self.normalized(modality="both")
        self.assertEqual(
            set(prepare_modality_config(config, "cDNA")["sample_to_fqid"]),
            {"rna_cell"},
        )
        self.assertEqual(
            set(prepare_modality_config(config, "gDNA")["sample_to_fqid"]),
            {"dna_cell", "other"},
        )

    def test_gdna_does_not_require_cdna_inputs_or_genes(self):
        self.create_layer("gDNA", [("dna_cell", "dna_id")])
        config = self.normalized(
            modality="gDNA",
            reference={"genome": str(self.genome)},
            ilse_info={"cDNA": {"metadata": "/missing/rna.csv"}},
        )
        layer = prepare_modality_config(config, "gDNA")
        self.assertEqual(layer["sample_to_fqid"], {"dna_cell": ["dna_id"]})
        self.assertFalse((self.root / "data" / "cDNA").exists())

    def test_cdna_does_not_validate_unselected_gdna_override(self):
        self.create_layer("cDNA", [("cell", "id")])
        config = self.normalized(ilse_info={"gDNA": {"metadata": "/missing/dna.csv"}})
        self.assertEqual(config["modality"], "cDNA")
        self.assertEqual(
            prepare_modality_config(config, "cDNA")["sample_to_fqid"],
            {"cell": ["id"]},
        )

    def test_both_reports_missing_requested_layer(self):
        self.create_layer("cDNA", [("cell", "id")])
        config = self.normalized(modality="both")
        prepare_modality_config(config, "cDNA")
        with self.assertRaisesRegex(ValueError, "gDNA: input directory does not exist"):
            prepare_modality_config(config, "gDNA")

    def test_metadata_discovery_ignores_fastq_sidecars(self):
        metadata, fastq = self.create_layer("cDNA", [("cell", "id")])
        (fastq / "facility_metrics.csv").write_text("not,a,sample,table\n")
        (fastq / "id_R1.fastq.gz.md5").write_text("checksum\n")
        (metadata.parent / "README.txt").write_text("Sequencing delivery details\n")
        layer = prepare_modality_config(self.normalized(), "cDNA")
        self.assertEqual(layer["ilse_info"]["metadata"], [str(metadata)])
        self.assertEqual(layer["sample_to_fqid"], {"cell": ["id"]})

    def test_facility_reads_can_be_flat_or_nested(self):
        directory = self.root / "data" / "cDNA"
        rows = [("cell_a", "flat"), ("cell_b", "nested"), ("cell_c", "deep")]
        self.write_metadata(directory / "metadata.csv", rows)
        fastq = directory / "fastq"
        expected = {
            "flat": self.write_pair(fastq, "flat"),
            "nested": self.write_pair(fastq / "nested", "nested"),
            "deep": self.write_pair(fastq / "deep" / "fastq", "deep"),
        }
        layer = prepare_modality_config(self.normalized(), "cDNA")
        self.assertEqual(layer["fqid_to_reads"], expected)

    def test_selected_sample_missing_read_two_fails(self):
        _, fastq = self.create_layer("cDNA", [("cell", "id")])
        (fastq / "id_R2.fastq.gz").unlink()
        with self.assertRaisesRegex(ValueError, "cDNA: could not find a complete FASTQ pair"):
            prepare_modality_config(self.normalized(), "cDNA")

    def test_excluded_sample_does_not_require_fastqs(self):
        _, fastq = self.create_layer("cDNA", [("keep", "present"), ("skip", "missing")])
        (fastq / "missing_R2.fastq.gz").unlink()
        layer = prepare_modality_config(self.normalized(samples=["keep"]), "cDNA")
        self.assertEqual(layer["sample_to_fqid"], {"keep": ["present"]})

    def test_layer_filters_override_global_filter_independently(self):
        rows = [("cell_1", "id_1"), ("cell_2", "id_2"), ("cell_10", "id_10")]
        for modality in ("cDNA", "gDNA"):
            self.create_layer(modality, rows)
        config = self.normalized(
            modality="both",
            samples=["would_match_neither_layer"],
            ilse_info={"cDNA": {"samples": ["cell_1"]}, "gDNA": {"samples": "cell_2"}},
        )
        self.assertEqual(
            prepare_modality_config(config, "cDNA")["sample_to_fqid"],
            {"cell_1": ["id_1"]},
        )
        self.assertEqual(
            prepare_modality_config(config, "gDNA")["sample_to_fqid"],
            {"cell_2": ["id_2"]},
        )

    def test_empty_layer_filter_includes_all_despite_global_filter(self):
        self.create_layer("cDNA", [("one", "id_1"), ("two", "id_2")])
        config = self.normalized(samples=["one"], ilse_info={"cDNA": {"samples": []}})
        layer = prepare_modality_config(config, "cDNA")
        self.assertEqual(set(layer["sample_to_fqid"]), {"one", "two"})

    def test_legacy_flat_config_infers_gdna(self):
        metadata = self.write_metadata(self.root / "old_delivery" / "metadata.csv", [("cell", "id")])
        fastq = self.root / "old_delivery" / "reads"
        self.write_pair(fastq, "id")
        config = self.normalized(
            pipeline="biscuit_methscan",
            reference={"genome": str(self.genome)},
            ilse_info={"metadata": str(metadata), "fastqdir": str(fastq)},
        )
        self.assertEqual(config["modality"], "gDNA")
        self.assertEqual(config["ilse_info"]["gDNA"]["metadata"], str(metadata))
        self.assertEqual(
            prepare_modality_config(config, "gDNA")["sample_to_fqid"],
            {"cell": ["id"]},
        )

    def test_both_rejects_ambiguous_flat_metadata_config(self):
        metadata, fastq = self.create_layer("cDNA", [("cell", "id")])
        with self.assertRaisesRegex(ValueError, "requires separate ilse_info"):
            self.normalized(
                modality="both", ilse_info={"metadata": str(metadata), "fastqdir": str(fastq)}
            )

    def test_resequenced_sample_keeps_multiple_ids_and_deduplicates_rows(self):
        metadata, fastq = self.create_layer(
            "cDNA", [("cell", "run_1"), ("cell", "run_1"), ("cell", "run_2")]
        )
        extra = self.write_metadata(metadata.parent / "second.csv", [("cell", "run_2")])
        config = self.normalized(
            ilse_info={"cDNA": {"metadata": f"{metadata} {extra}", "fastqdir": str(fastq)}}
        )
        layer = prepare_modality_config(config, "cDNA")
        self.assertEqual(layer["sample_to_fqid"], {"cell": ["run_1", "run_2"]})
        self.assertEqual(set(layer["fqid_to_reads"]), {"run_1", "run_2"})

    def test_one_fastq_id_cannot_belong_to_two_samples_in_a_layer(self):
        self.create_layer("cDNA", [("one", "id"), ("two", "id")])
        with self.assertRaisesRegex(ValueError, "assigned to both"):
            prepare_modality_config(self.normalized(), "cDNA")

    def test_tsv_preserves_leading_zero_sample_and_fastq_ids(self):
        directory = self.root / "data" / "gDNA"
        self.write_metadata(directory / "metadata.tsv", [("001", "0007")], delimiter="\t")
        self.write_pair(directory / "fastq", "0007")
        config = self.normalized(modality="gDNA", reference={"genome": str(self.genome)})
        self.assertEqual(
            prepare_modality_config(config, "gDNA")["sample_to_fqid"],
            {"001": ["0007"]},
        )

    def test_no_sample_matches_reports_clear_error(self):
        self.create_layer("cDNA", [("one", "id")])
        with self.assertRaisesRegex(ValueError, "no samples with FASTQ IDs were selected"):
            prepare_modality_config(self.normalized(samples=["other"]), "cDNA")

    def test_invalid_metadata_headers_report_required_columns(self):
        metadata, _ = self.create_layer("cDNA", [("one", "id")])
        metadata.write_text("Name,Identifier\none,id\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Missing column 'Sample Name'"):
            prepare_modality_config(self.normalized(), "cDNA")

    def test_cdna_requires_gene_annotations(self):
        for modality in ("cDNA", "both"):
            with self.subTest(modality=modality):
                with self.assertRaisesRegex(ValueError, "reference.genes"):
                    self.normalized(modality=modality, reference={"genome": str(self.genome)})
        self.normalized(pipeline="star_umite")

    def test_removed_mapper_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unexpected pipeline"):
            self.normalized(pipeline="salmon")

    def test_active_gdna_requires_methscan_arguments_to_be_strings(self):
        for modality in ("gDNA", "both"):
            for step in ("prepare", "filter", "smooth", "scan", "matrix"):
                key = f"methscan_{step}_args"
                for value in (None, 1, False, [], {}):
                    with self.subTest(modality=modality, key=key, value=value):
                        with self.assertRaisesRegex(ValueError, key):
                            self.normalized(modality=modality, **{key: value})

    def test_cdna_ignores_unused_methscan_argument_types(self):
        self.create_layer("cDNA", [("cell", "id")])
        unused_arguments = {
            f"methscan_{step}_args": {"unused": True}
            for step in ("prepare", "filter", "smooth", "scan", "matrix")
        }
        config = self.normalized(modality="cDNA", **unused_arguments)
        layer = prepare_modality_config(config, "cDNA")
        self.assertEqual(layer["sample_to_fqid"], {"cell": ["id"]})

    def test_methscan_argument_strings_survive_layer_preparation(self):
        self.create_layer("gDNA", [("cell", "id")])
        arguments = {
            "methscan_prepare_args": "--input-format biscuit_short --chunksize 1000",
            "methscan_filter_args": "--min-sites 250",
            "methscan_smooth_args": "",
            "methscan_scan_args": "--min-cells 3",
            "methscan_matrix_args": "--threads 1",
        }
        for modality in ("gDNA", "both"):
            with self.subTest(modality=modality):
                layer = prepare_modality_config(
                    self.normalized(modality=modality, **arguments), "gDNA"
                )
                for key, value in arguments.items():
                    self.assertEqual(layer[key], value)

    def test_active_gdna_rejects_incompatible_output_options(self):
        incompatible = {
            "methscan_matrix_args": "--sparse",
            "methscan_scan_args": "--write-header",
        }
        for modality in ("gDNA", "both"):
            for key, arguments in incompatible.items():
                with self.subTest(modality=modality, key=key):
                    with self.assertRaisesRegex(ValueError, key):
                        self.normalized(modality=modality, **{key: arguments})
        self.normalized(modality="cDNA", **incompatible)

    def test_active_gdna_rejects_malformed_argument_quoting(self):
        for modality in ("gDNA", "both"):
            for step in ("prepare", "filter", "smooth", "scan", "matrix"):
                key = f"methscan_{step}_args"
                with self.subTest(modality=modality, key=key):
                    with self.assertRaisesRegex(ValueError, key):
                        self.normalized(modality=modality, **{key: "--option 'unfinished"})

    def test_filter_cell_name_file_is_resolved_from_both_option_forms(self):
        cell_names = self.root / "selected_cells.txt"
        cell_names.write_text("cell_1\n", encoding="utf-8")
        for modality in ("gDNA", "both"):
            for separator in (" ", "="):
                with self.subTest(modality=modality, separator=separator):
                    arguments = f"--cell-names{separator}{cell_names} --discard"
                    config = self.normalized(
                        modality=modality, methscan_filter_args=arguments
                    )
                    self.assertEqual(config["methscan_filter_cell_names"], [str(cell_names)])
                    self.assertEqual(config["methscan_filter_args"], arguments)
        self.assertEqual(self.normalized(modality="gDNA")["methscan_filter_cell_names"], [])

    def test_filter_cell_name_option_requires_an_existing_absolute_file(self):
        invalid = (
            "--cell-names",
            "--cell-names=",
            "--cell-names --discard",
            "--cell-names relative.txt",
            f"--cell-names={self.root / 'missing.txt'}",
            f"--cell-names {self.root}",
        )
        for modality in ("gDNA", "both"):
            for arguments in invalid:
                with self.subTest(modality=modality, arguments=arguments):
                    with self.assertRaisesRegex(ValueError, "cell.names"):
                        self.normalized(modality=modality, methscan_filter_args=arguments)
        self.normalized(modality="cDNA", methscan_filter_args="--cell-names /missing.txt")

    def test_normalization_preserves_unknown_keys_and_original_config(self):
        original = copy.deepcopy(self.config)
        original["custom_setting"] = {"values": [1, 2]}
        before = copy.deepcopy(original)
        normalized = normalize_config(original, self.root)
        self.assertEqual(original, before)
        self.assertEqual(normalized["custom_setting"], before["custom_setting"])
        normalized["custom_setting"]["values"].append(3)
        self.assertEqual(original, before)

    def test_invalid_config_values_fail_with_value_error(self):
        invalid = [
            {"dataset": "../outside"},
            {"dataset": []},
            {"modality": "RNA"},
            {"modality": ["cDNA", "gDNA"]},
            {"pipeline": "unknown"},
            {"pipeline": "biscuit_methscan", "modality": "both"},
            {"datadir": "relative/data"},
            {"outdir": "relative/results"},
            {"reference": []},
            {"reference": {"genome": str(self.root), "genes": str(self.genes)}},
            {"ilse_info": ["invalid"]},
        ]
        for overrides in invalid:
            with self.subTest(overrides=overrides):
                with self.assertRaises(ValueError):
                    self.normalized(**overrides)

    def test_invalid_selected_layer_overrides_and_filters_fail(self):
        self.create_layer("cDNA", [("one", "id")])
        invalid = [
            {"cDNA": ["invalid"]},
            {"cDNA": {"metadata": 42}},
            {"cDNA": {"fastqdir": "relative/fastq"}},
            {"cDNA": {"samples": [3]}},
            {"cDNA": {"samples": "["}},
        ]
        for ilse_info in invalid:
            with self.subTest(ilse_info=ilse_info):
                with self.assertRaises(ValueError):
                    prepare_modality_config(self.normalized(ilse_info=ilse_info), "cDNA")


if __name__ == "__main__":
    unittest.main()
