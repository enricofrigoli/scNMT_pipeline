"""Validate configuration and resolve independent batches and demultiplexed layers."""

from collections import defaultdict
from copy import deepcopy
from pathlib import Path
import gzip
import os
import re
import shlex

import pandas as pd


MODALITIES = ("cDNA", "gDNA")
PIPELINES = ("star_umite", "biscuit_methscan")
METHSCAN_STEPS = ("prepare", "filter", "smooth", "scan", "matrix", "profile")
# umicount --stranded modes; "umi" strand-resolves only the UMI-containing
# readpairs, since SmartSeq3 internal fragments are not strand-specific.
UMICOUNT_STRAND_MODES = ("no", "umi", "yes", "reverse")
METADATA_EXTENSIONS = (".xls", ".xlsx", ".csv", ".tsv")
STAR_INDEX_FILES = ("Genome", "SA", "SAindex", "genomeParameters.txt")
BISCUIT_INDEX_SUFFIXES = (
    ".bis.amb", ".bis.ann", ".bis.pac", ".dau.bwt", ".dau.sa", ".par.bwt", ".par.sa",
)
STAR_INDEX_DIRNAME = "star_index"
READ_LENGTH_SAMPLE_FILES = 8
READ_LENGTH_SAMPLE_RECORDS = 1000
BATCH_OPTIONS = {
    "modality", "samples", "reference", "ilse_info",
    "umiextract_args", "umicount_args", "star_index_args", "star_align_args",
    "biscuit_index_alg", "biscuit_align_args", "biscuit_pileup_args",
    "generate_QC_plots",
    *(f"methscan_{step}_args" for step in METHSCAN_STEPS),
}


def absolute_path(value: str, label: str) -> Path:
    """Validate a configured path without requiring an output to exist yet."""
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty absolute path")
    path = Path(value)
    if not path.is_absolute() or any(char.isspace() for char in value):
        raise ValueError(f"{label} must be an absolute path without whitespace: {value}")
    return path


def validate_identifier(value: str, label: str) -> str:
    """Keep sample and dataset names safe as path components and shell arguments."""
    if (
        not isinstance(value, str)
        or not re.fullmatch(r"[\w.-]+", value)
        or value in (".", "..")
    ):
        raise ValueError(
            f"{label} must contain only letters, digits, underscores, dots or hyphens"
        )
    return value


def normalize_config(raw_config: dict, base_dir: str | Path) -> dict:
    """Resolve defaults and legacy settings, validating active references."""
    config = deepcopy(raw_config)
    config["dataset"] = validate_identifier(config.get("dataset"), "dataset")
    pipeline = config.get("pipeline", "star_umite")
    if pipeline not in PIPELINES:
        raise ValueError(f"Unexpected pipeline {pipeline!r}; choose star_umite or biscuit_methscan")
    config["pipeline"] = pipeline
    modality = config.get(
        "modality", "gDNA" if pipeline == "biscuit_methscan" else "cDNA"
    )
    if modality not in (*MODALITIES, "both"):
        raise ValueError("modality must be cDNA, gDNA or both")
    config["modality"] = modality
    selected = list(MODALITIES) if modality == "both" else [modality]
    if "cDNA" in selected and pipeline != "star_umite":
        raise ValueError(
            "cDNA requires pipeline: star_umite "
            "(also for modality: both)"
        )

    generate_qc_plots = config.get("generate_QC_plots", False)
    if not isinstance(generate_qc_plots, bool):
        raise ValueError("generate_QC_plots must be true or false")
    config["generate_QC_plots"] = generate_qc_plots

    base_dir = Path(base_dir).resolve()
    for key, default in (("datadir", "data"), ("outdir", "results")):
        config[key] = str(absolute_path(config.get(key, str(base_dir / default)), key))

    reference = config.get("reference")
    if not isinstance(reference, dict):
        raise ValueError('No "reference" mapping found in config file')
    required = []
    if "gDNA" in selected or ("cDNA" in selected and not reference.get("star_index")):
        required.append("genome")
    if "cDNA" in selected:
        required.append("genes")
    for key in required:
        path = absolute_path(reference.get(key), f"reference.{key}")
        if not path.is_file():
            raise ValueError(f"reference.{key} file does not exist: {path}")
    if generate_qc_plots and "gDNA" in selected:
        reference["profile_regions"] = validate_profile_regions(reference)

    if "cDNA" in selected and "star_index" in reference:
        index_dir = absolute_path(reference["star_index"], "reference.star_index").resolve()
        if not index_dir.is_dir():
            raise ValueError(f"reference.star_index directory does not exist: {index_dir}")
        validate_index_files([index_dir / name for name in STAR_INDEX_FILES], "reference.star_index")
        reference["star_index"] = str(index_dir)
    if "gDNA" in selected and "biscuit_index" in reference:
        prefix = absolute_path(reference["biscuit_index"], "reference.biscuit_index")
        validate_index_files(
            [Path(str(prefix) + suffix) for suffix in BISCUIT_INDEX_SUFFIXES],
            "reference.biscuit_index (file prefix)",
        )
        reference["biscuit_index"] = str(prefix)

    if "gDNA" in selected:
        config["methscan_filter_cell_names"] = validate_methscan_args(config)
    if "cDNA" in selected:
        config["umicount_stranded"] = validate_umicount_args(config)

    ilse_info = config.get("ilse_info") or {}
    if not isinstance(ilse_info, dict):
        raise ValueError("ilse_info must be a mapping keyed by cDNA and/or gDNA")
    if "metadata" in ilse_info or "fastqdir" in ilse_info:
        if modality == "both":
            raise ValueError(
                "modality: both requires separate ilse_info.cDNA and "
                "ilse_info.gDNA overrides; remove flat ilse_info to use data/ defaults"
            )
        if any(layer in ilse_info for layer in MODALITIES):
            raise ValueError("Do not mix flat and modality-specific ilse_info settings")
        ilse_info = {modality: ilse_info}
    config["ilse_info"] = ilse_info
    return config


def validate_index_files(paths: list[Path], label: str) -> None:
    """Check required index components without rebuilding or modifying them."""
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise ValueError(f"{label} is missing index files: {', '.join(missing)}")


def validate_profile_regions(reference: dict) -> dict:
    """Resolve the optional region files profiled for QC, keyed by plot name."""
    regions = reference.get("profile_regions") or {}
    if not isinstance(regions, dict):
        raise ValueError(
            "reference.profile_regions must be a mapping of names to BED files"
        )
    resolved = {}
    for name, value in regions.items():
        # The name becomes a wildcard, a path component and a plot title.
        label = validate_identifier(str(name), "reference.profile_regions name")
        path = absolute_path(value, f"reference.profile_regions.{label}")
        if not path.is_file():
            raise ValueError(
                f"reference.profile_regions.{label} file does not exist: {path}"
            )
        resolved[label] = str(path)
    return resolved


def qc_targets(layer: dict, modality: str) -> list[str]:
    """List the QC files a prepared layer produces when QC is enabled."""
    qc_dir = Path(layer["outdir"]) / "qc"
    targets = [
        str(qc_dir / "multiqc/multiqc_report.html"),
        str(qc_dir / f"tables/{layer['dataset']}.{modality}.per_cell_metrics_mqc.tsv"),
    ]
    if modality != "gDNA":
        return targets
    targets.append(str(qc_dir / "cell_stats.png"))
    targets.append(str(qc_dir / "plate/plate_qc.png"))
    regions = layer["reference"].get("profile_regions") or {}
    targets.extend(
        str(qc_dir / f"profiles/{region}_profiles.pdf") for region in regions
    )
    return targets


def shared_star_index_dir(genome: str) -> str:
    """Generate one index per genome instead of one per batch."""
    genome_dir = Path(genome).resolve().parent
    if not os.access(genome_dir, os.W_OK):
        raise ValueError(
            "cannot generate a STAR index beside a read-only reference genome: "
            f"{genome_dir}. Supply reference.star_index to reuse an existing index."
        )
    return str(genome_dir / STAR_INDEX_DIRNAME)


def detect_read_length(fastq: str, max_records: int = READ_LENGTH_SAMPLE_RECORDS) -> int:
    """Measure the longest read among the leading records of one gzipped FASTQ."""
    longest = 0
    with gzip.open(fastq, "rt") as handle:
        for index, line in enumerate(handle):
            if index >= max_records * 4:
                break
            if index % 4 == 1:
                longest = max(longest, len(line.strip()))
    return longest


def detect_sjdb_overhang(fqid_to_reads: dict, modality: str) -> int:
    """Apply STAR's max(read length) - 1 guidance by sampling a few read pairs."""
    fqids = sorted(fqid_to_reads)
    step = max(1, len(fqids) // READ_LENGTH_SAMPLE_FILES)
    longest = 0
    for fqid in fqids[::step][:READ_LENGTH_SAMPLE_FILES]:
        for read in ("read1", "read2"):
            longest = max(longest, detect_read_length(fqid_to_reads[fqid][read]))
    if longest < 2:
        raise ValueError(
            f"{modality}: could not measure a read length from the FASTQ files, so "
            "--sjdbOverhang is unknown; supply reference.star_index to skip generation"
        )
    return longest - 1


def normalize_run_config(raw_config: dict, base_dir: str | Path) -> dict:
    """Validate a batch manifest, or preserve the legacy single-dataset mode."""
    if "batches" not in raw_config:
        return normalize_config(raw_config, base_dir)
    config = deepcopy(raw_config)
    batches = config["batches"]
    if not isinstance(batches, dict) or not batches:
        raise ValueError("batches must be a non-empty mapping of batch names to settings")
    if config.get("ilse_info"):
        raise ValueError(
            "With batches, put ilse_info overrides inside each batch; "
            "global ilse_info would reuse inputs across batches"
        )
    if "modality" in config:
        raise ValueError("With batches, set modality inside each batch")
    if config.get("pipeline", "star_umite") not in PIPELINES:
        raise ValueError("Unexpected pipeline; choose star_umite or biscuit_methscan")
    reference = config.get("reference", {})
    if not isinstance(reference, dict):
        raise ValueError("reference must be a mapping")
    base_dir = Path(base_dir).resolve()
    for key, default in (("datadir", "data"), ("outdir", "results")):
        config[key] = str(absolute_path(config.get(key, str(base_dir / default)), key))
    for batch, settings in batches.items():
        if not isinstance(batch, str) or not re.fullmatch(r"[A-Za-z0-9_]+", batch):
            raise ValueError("Batch names must contain only ASCII letters, digits and underscores")
        if not isinstance(settings, dict):
            raise ValueError(f"Batch {batch}: settings must be a mapping")
        if settings.get("modality") not in (*MODALITIES, "both"):
            raise ValueError(f"Batch {batch}: modality must be cDNA, gDNA or both")
        unknown = settings.keys() - BATCH_OPTIONS
        if unknown:
            raise ValueError(
                f"Batch {batch}: unsupported settings: "
                f"{', '.join(sorted(str(key) for key in unknown))}"
            )
        if "reference" in settings and not isinstance(settings["reference"], dict):
            raise ValueError(f"Batch {batch}: reference must be a mapping")
    return config


def prepare_batch_configs(config: dict, base_dir: str | Path) -> dict:
    """Resolve independent input roots and settings before preparing modalities."""
    config = normalize_run_config(config, base_dir)
    if "batches" not in config:
        return {None: config}
    common = {key: value for key, value in config.items() if key != "batches"}
    batch_configs = {}
    for batch, settings in config["batches"].items():
        batch_config = deepcopy(common) | deepcopy(settings)
        batch_config["dataset"] = batch
        batch_config["datadir"] = str(Path(config["datadir"]) / batch)
        batch_config["outdir"] = str(Path(config["outdir"]) / batch)
        batch_config["pipeline"] = (
            "biscuit_methscan" if settings["modality"] == "gDNA" else "star_umite"
        )
        batch_config["reference"] = (
            deepcopy(config.get("reference", {}))
            | deepcopy(settings.get("reference", {}))
        )
        batch_config["ilse_info"] = deepcopy(settings.get("ilse_info", {}))
        try:
            batch_configs[batch] = normalize_config(batch_config, base_dir)
        except ValueError as exc:
            raise ValueError(f"Batch {batch}: {exc}") from exc
    return batch_configs


def validate_umicount_args(config: dict) -> str:
    """Validate counting options and resolve the strand mode shared with the GTF dump."""
    value = config.get("umicount_args", "")
    if not isinstance(value, str):
        raise ValueError(
            'umicount_args must be a string of command-line options; '
            'use "" for no extra options'
        )
    try:
        tokens = shlex.split(value)
    except ValueError as exc:
        raise ValueError(f"Invalid quoting in umicount_args: {exc}") from exc

    # The workflow supplies these itself, from the rule's inputs and outputs.
    reserved = {
        "-f": "--bams is built from the aligned BAMs",
        "--bams": "--bams is built from the aligned BAMs",
        "-d": "the output directory is fixed by the rule",
        "--output_dir": "the output directory is fixed by the rule",
        "-c": "cores come from the Snakemake thread allocation",
        "--cores": "cores come from the Snakemake thread allocation",
        "-l": "the log path is fixed by the rule",
        "--logfile": "the log path is fixed by the rule",
        "-g": "the GTF comes from reference.genes",
        "--gtf": "the GTF comes from reference.genes",
        "--GTF_dump": "the dump is written by the parse_dump_GTF rule",
        "--GTF_skip_parse": "the dump is read from the parse_dump_GTF rule",
        "--combine_unspliced": "the rule always combines spliced and unspliced UMIs",
        "--no_dedup": "the workflow requires deduplicated UMI and duplicate-count matrices",
    }
    stranded = "no"
    for index, token in enumerate(tokens):
        flag, separator, inline = token.partition("=")
        if flag in reserved:
            raise ValueError(f"umicount_args cannot use {flag}: {reserved[flag]}")
        if flag != "--stranded":
            continue
        if separator:
            stranded = inline
        elif index + 1 == len(tokens) or tokens[index + 1].startswith("-"):
            raise ValueError(
                "umicount_args --stranded requires a mode: "
                f"{', '.join(UMICOUNT_STRAND_MODES)}"
            )
        else:
            stranded = tokens[index + 1]
        if stranded not in UMICOUNT_STRAND_MODES:
            raise ValueError(
                f"umicount_args --stranded must be one of "
                f"{', '.join(UMICOUNT_STRAND_MODES)}: {stranded}"
            )
    return stranded


def validate_methscan_args(config: dict) -> list[str]:
    """Validate stage options and identify external inputs used by filtering."""
    tokens_by_step = {}
    for step in METHSCAN_STEPS:
        key = f"methscan_{step}_args"
        value = config.get(key, "")
        if not isinstance(value, str):
            raise ValueError(
                f'{key} must be a string of command-line options; '
                'use "" for no extra options'
            )
        try:
            tokens_by_step[step] = shlex.split(value)
        except ValueError as exc:
            raise ValueError(f"Invalid quoting in {key}: {exc}") from exc

    incompatible = (
        ("matrix", "--sparse", "H5AD conversion requires the default CSV matrices"),
        ("scan", "--write-header", "the downstream matrix step requires a BED without a header"),
    )
    for step, flag, reason in incompatible:
        if any(token.split("=", 1)[0] == flag for token in tokens_by_step[step]):
            raise ValueError(f"methscan_{step}_args cannot use {flag}: {reason}")

    cell_names = None
    tokens = tokens_by_step["filter"]
    for index, token in enumerate(tokens):
        if token == "--cell-names":
            if index + 1 == len(tokens) or tokens[index + 1].startswith("--"):
                raise ValueError("methscan_filter_args --cell-names requires a file path")
            cell_names = tokens[index + 1]
        elif token.startswith("--cell-names="):
            cell_names = token.split("=", 1)[1]
    if cell_names is None:
        return []
    path = absolute_path(cell_names, "methscan_filter_args --cell-names")
    if not path.is_file():
        raise ValueError(f"methscan_filter_args --cell-names file does not exist: {path}")
    return [str(path)]


def path_list(value: str | list, label: str, *, directory: bool = False) -> list[str]:
    """Accept existing whitespace-separated paths and YAML lists."""
    paths = value.split() if isinstance(value, str) else value
    if not isinstance(paths, list) or not paths:
        raise ValueError(f"No paths provided for {label}")
    result = []
    for value in paths:
        path = absolute_path(value, label)
        if not (path.is_dir() if directory else path.is_file()):
            kind = "directory" if directory else "file"
            raise ValueError(f"{label} {kind} does not exist: {path}")
        result.append(str(path))
    return result


def sample_patterns(value: str | list | None, modality: str) -> list:
    """Compile optional full-match sample filters independently for each layer."""
    value = value or []
    if isinstance(value, str):
        value = [value]
    if not isinstance(value, list) or any(not isinstance(p, str) for p in value):
        raise ValueError(
            f"{modality}: samples must be a regex string or a list of regex strings"
        )
    try:
        return [re.compile(pattern) for pattern in value]
    except re.error as exc:
        raise ValueError(f"{modality}: invalid sample filter regex: {exc}") from exc


def read_metadata(filename: str) -> pd.DataFrame:
    """Normalize facility per-read tables and legacy per-pair sample tables."""
    suffix = Path(filename).suffix.lower()
    if suffix in (".xls", ".xlsx"):
        table = pd.read_excel(filename, dtype=str, keep_default_na=False)
    elif suffix in (".csv", ".tsv"):
        table = pd.read_csv(
            filename, sep="\t" if suffix == ".tsv" else ",",
            dtype=str, keep_default_na=False,
        )
    else:
        raise ValueError(f"Unexpected metadata extension {suffix!r}: {filename}")
    table.columns = table.columns.str.strip()
    if table.columns.duplicated().any():
        raise ValueError(f"Duplicate metadata column names in {filename}")

    facility_format = "SAMPLE_NAME" in table.columns or "FASTQ_FILE" in table.columns
    required = (
        ("SAMPLE_NAME", "FASTQ_FILE", "READ") if facility_format
        else ("Sample Name", "Unique ID / Lane")
    )
    for column in required:
        if column not in table.columns:
            raise ValueError(
                f"Missing column {column!r} in {filename}. "
                f"Found: {list(table.columns)}. Expected SAMPLE_NAME, FASTQ_FILE "
                "and READ, or legacy Sample Name and Unique ID / Lane."
            )
    if not facility_format:
        return table

    # The facility includes an unnamed pair of reads that could not be assigned
    # to a cell. These are not cell libraries and must not enter the workflow.
    unnamed = table["SAMPLE_NAME"].str.strip().eq("")
    undetermined = table["FASTQ_FILE"].str.startswith("Undetermined")
    table = table.loc[~(unnamed & undetermined)].copy()
    if table["SAMPLE_NAME"].str.strip().eq("").any():
        raise ValueError(f"Missing SAMPLE_NAME for a named FASTQ file in {filename}")

    paired_ids = []
    for _, row in table.iterrows():
        match = re.fullmatch(r"([\w.-]+)_R([12])\.fastq\.gz", row["FASTQ_FILE"])
        if match is None:
            raise ValueError(
                f"Unexpected FASTQ_FILE {row['FASTQ_FILE']!r} in {filename}; "
                "expected <FASTQ ID>_R1.fastq.gz or <FASTQ ID>_R2.fastq.gz"
            )
        fqid, mate = match.groups()
        if row["READ"].strip() != mate:
            raise ValueError(
                f"READ {row['READ']!r} does not match FASTQ_FILE "
                f"{row['FASTQ_FILE']!r} in {filename}"
            )
        paired_ids.append(fqid)

    normalized = {
        "Sample Name": table["SAMPLE_NAME"],
        "Unique ID / Lane": pd.Series(paired_ids, index=table.index, dtype=str),
    }
    for column, values in normalized.items():
        if column in table.columns and not table[column].equals(values):
            raise ValueError(f"Conflicting facility and legacy column {column!r} in {filename}")
        table[column] = values
    return table


def find_fastq_pair(fqid: str, searchdirs: list[str], modality: str) -> dict[str, str]:
    """Resolve one exact read pair, ignoring facility sidecar files."""
    for searchdir in searchdirs:
        root = Path(searchdir)
        for directory in (root / fqid / "fastq", root / fqid, root):
            reads = {
                f"read{read}": directory / f"{fqid}_R{read}.fastq.gz"
                for read in (1, 2)
            }
            if all(path.is_file() for path in reads.values()):
                return {key: str(path) for key, path in reads.items()}
    raise ValueError(
        f"{modality}: could not find a complete FASTQ pair for ID {fqid!r}; "
        f"expected {fqid}_R1.fastq.gz and {fqid}_R2.fastq.gz in {searchdirs} "
        f"(directly, under {fqid}/, or under {fqid}/fastq/)"
    )


def prepare_modality_config(config: dict, modality: str) -> dict:
    """Build one isolated workflow's config and sample-to-read lookup tables."""
    layer_config = deepcopy(config)
    layer_config["modality"] = modality
    layer_config["pipeline"] = (
        config["pipeline"] if modality == "cDNA" else "biscuit_methscan"
    )
    layer_config["outdir"] = str(Path(config["outdir"]) / modality)
    reference_dir = Path(layer_config["outdir"]) / "reference"
    if modality == "cDNA":
        external_index = config["reference"].get("star_index")
        index_dir = external_index or shared_star_index_dir(config["reference"]["genome"])
        layer_config["star_index_dir"] = index_dir
        layer_config["star_index_generated"] = external_index is None
        layer_config["star_index_inputs"] = [
            str(Path(index_dir) / name) for name in STAR_INDEX_FILES
        ]
    else:
        external_index = config["reference"].get("biscuit_index")
        index_dir = reference_dir / "biscuit_index"
        prefix = external_index or str(index_dir / "genome.fa")
        layer_config["biscuit_index_prefix"] = prefix
        layer_config["biscuit_index_inputs"] = (
            [prefix + suffix for suffix in BISCUIT_INDEX_SUFFIXES]
            if external_index else [str(index_dir)]
        )
    layer_dir = Path(config["datadir"]) / modality
    overrides = config["ilse_info"].get(modality) or {}
    if not isinstance(overrides, dict):
        raise ValueError(f"ilse_info.{modality} must be a mapping")
    metadata = overrides.get("metadata")
    if metadata is None:
        if not layer_dir.is_dir():
            raise ValueError(f"{modality}: input directory does not exist: {layer_dir}")
        metadata = [
            str(path) for path in sorted(layer_dir.iterdir())
            if path.is_file() and path.suffix.lower() in METADATA_EXTENSIONS
        ]
        if not metadata:
            raise ValueError(
                f"{modality}: no metadata tables found directly in {layer_dir}; "
                "expected .xls, .xlsx, .csv or .tsv"
            )
    metadata = path_list(metadata, f"{modality} metadata")
    fastqdirs = path_list(
        overrides.get("fastqdir", [str(layer_dir / "fastq")]),
        f"{modality} FASTQ", directory=True,
    )
    layer_config["ilse_info"] = {"metadata": metadata, "fastqdir": fastqdirs}
    layer_config["samples"] = overrides.get("samples", config.get("samples", []))
    patterns = sample_patterns(layer_config["samples"], modality)

    sample_to_fqid = defaultdict(list)
    fqid_to_reads = {}
    fqid_to_sample = {}
    for filename in metadata:
        for _, row in read_metadata(filename).iterrows():
            if pd.isna(row["Unique ID / Lane"]) or not row["Unique ID / Lane"].strip():
                continue
            sample = row["Sample Name"]
            if patterns and not any(p.fullmatch(str(sample)) for p in patterns):
                continue
            sample = validate_identifier(sample, f"{modality} sample name in {filename}")
            fqid = validate_identifier(
                row["Unique ID / Lane"].strip(), f"{modality} FASTQ ID in {filename}"
            )
            if fqid in fqid_to_sample and fqid_to_sample[fqid] != sample:
                raise ValueError(
                    f"{modality}: FASTQ ID {fqid!r} is assigned to both "
                    f"{fqid_to_sample[fqid]!r} and {sample!r}"
                )
            if fqid not in fqid_to_reads:
                fqid_to_reads[fqid] = find_fastq_pair(fqid, fastqdirs, modality)
                fqid_to_sample[fqid] = sample
            if fqid not in sample_to_fqid[sample]:
                sample_to_fqid[sample].append(fqid)
    if not sample_to_fqid:
        raise ValueError(
            f"{modality}: no samples with FASTQ IDs were selected. Check metadata "
            "and samples filters; omit samples or set samples: [] to include all samples."
        )
    layer_config["sample_to_fqid"] = dict(sample_to_fqid)
    layer_config["fqid_to_reads"] = fqid_to_reads
    if layer_config.get("star_index_generated"):
        layer_config["star_sjdb_overhang"] = detect_sjdb_overhang(fqid_to_reads, modality)
    return layer_config
