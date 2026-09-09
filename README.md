This repo contains a collection of pipelines to process scNMT sequencing data, intended for DKFZ A290 internal use.



# Important notes

- A lot of steps of the pipeline have issues with the Windows filesystem (e.g. the O drive). As a result, it is recommended to have **the pipeline folder, the output directory and the reference files all on a Linux volume** (e.g. scratch).
- Paths in the config file (see below) must be **absolute** paths and must not contain whitespaces.
- STAR requires a significant amount of RAM (up to 10 times the size of the unzipped genome) at indexing step. **Make sure your machine has enough RAM left.**

# Preparation
## Environment setup

The pipeline assumes avaiability of the following commands.

- `conda`
- `snakemake` >=9.19 (required for batch modules)

One way to get `conda` is to install Miniforge, following instructions on [their GitHub repo](https://github.com/conda-forge/miniforge), but feel free to choose any alternative.

Once `conda` is available, create a conda environment called "snakemake" and install necessary packages, which are `snakemake`, `pandas`, `xlrd` and `openpyxl`. To ensure portability, you can use the the provided environment file:

```bash
conda env create -n snakemake -f envs/snakemake.yaml
```

## Reference genome and annotation

For most popular species, these files can be downloaded from [EMSELBL](https://www.ensembl.org/index.html).

For mouse GRCm39 release 115, use the following commands:

```bash
curl -O https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz
curl -O https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz
```

Similarly, for human GRCh38 release 115, use:

```bash
curl -O https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
curl -O https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
```

Decompress the downloaded references before configuring the pipeline, for example with `gunzip` on each `.fa.gz` and `.gtf.gz` file. Use the resulting `.fa` and `.gtf` paths in `snakeconfig.yaml`.

## Input layout

Organize data by batch, with separate modality folders inside each batch. For example, `batch_01` has both modalities, `batch_02` only cDNA, and `batch_03` only gDNA:

```bash
mkdir -p data/batch_01/cDNA/fastq data/batch_01/gDNA/fastq data/batch_02/cDNA/fastq data/batch_03/gDNA/fastq results
```

```text
scNMT_pipeline/
├── data/
│   ├── batch_01/
│   │   ├── cDNA/
│   │   │   ├── run1_meta.tsv
│   │   │   ├── run2_meta.tsv
│   │   │   └── fastq/
│   │   │       ├── <fqid>_R1.fastq.gz
│   │   │       ├── <fqid>_R2.fastq.gz
│   │   │       └── ... facility metadata and sidecar files
│   │   └── gDNA/
│   │       ├── metadata.tsv
│   │       └── fastq/
│   ├── batch_02/
│   │   └── cDNA/
│   │       ├── metadata.tsv
│   │       └── fastq/
│   └── batch_03/
│       └── gDNA/
│           ├── metadata.tsv
│           └── fastq/
└── results/
    ├── batch_01/
    │   ├── cDNA/
    │   └── gDNA/
    ├── batch_02/
    │   └── cDNA/
    └── batch_03/
        └── gDNA/
```

Only the modalities selected for a batch need data. Selecting `both` requires metadata and samples with FASTQ IDs in both modality folders; an empty or missing selected modality produces an error.

The pipeline discovers metadata tables (`.xls`, `.xlsx`, `.csv`, or `.tsv`) directly inside each selected `data/<batch>/<modality>/` folder. Keep unrelated tables inside `fastq/` so they are not interpreted as sample metadata. Facility metadata and sidecar files inside `fastq/` can remain in place.

Current facility tables, such as `48006_meta.tsv`, use `SAMPLE_NAME`, `FASTQ_FILE`, and `READ`, with one row for each mate:

| SAMPLE_NAME | FASTQ_FILE | READ |
| --- | --- | --- |
| cell_A | run1_cell_A_R1.fastq.gz | 1 |
| cell_A | run1_cell_A_R2.fastq.gz | 2 |

The pipeline derives the FASTQ ID (`run1_cell_A` in this example) by removing `_R1.fastq.gz` or `_R2.fastq.gz` and verifies that `READ` is `1` or `2` and matches the filename. The mate rows form one read pair. Legacy tables with `Sample Name` and `Unique ID / Lane` remain supported. The reader normalizes these formats internally without changing the metadata files. Rows with a blank `SAMPLE_NAME` are skipped only for `Undetermined*` FASTQs; other rows missing a sample name cause an error.

FASTQs are searched under the modality's `fastq/` folder, using either the flat layout shown above or the facility layout `<fqid>/fastq/<fqid>_R1.fastq.gz` and `<fqid>/fastq/<fqid>_R2.fastq.gz`.

Within one batch and modality, exactly matching `SAMPLE_NAME` values (or legacy `Sample Name` values) across metadata tables identify the same cell sequenced multiple times. Distinct FASTQ IDs are combined for that cell; repeated rows for the same FASTQ ID, including the two mate rows, contribute only one read pair. Keep sample names consistent across runs and use `samples` to select cells belonging to the dataset. Sample names, FASTQ IDs, and dataset names may contain only letters, digits, underscores, dots, and hyphens; sample names containing spaces remain unsupported. Different batches and modalities handle sample names and FASTQ IDs independently: identical names in different batches stay separate, and the pipeline does not intersect, match, or merge cells between modalities. Match cells manually using the generated outputs.

## Configuration

Copy `snakeconfig.example.yaml` to `snakeconfig.yaml`, set the shared reference paths, and list the batches to process:

```yaml
reference:
  genome: /absolute/path/to/genome.fa
  genes: /absolute/path/to/annotation.gtf
  star_index: /absolute/path/to/star_index
  biscuit_index: /absolute/path/to/biscuit/genome  # File prefix, not a directory.

batches:
  batch_01:
    modality: both
  batch_02:
    modality: cDNA
  batch_03:
    modality: gDNA

# Optional absolute overrides; defaults are data/ and results/ beside Snakefile.
# datadir: /absolute/path/to/data
# outdir: /absolute/path/to/results
```

Each batch must specify `modality: cDNA`, `gDNA`, or `both`. Batch IDs allow only letters, digits, and underscores; they name the input folders, output folders, and final matrices. cDNA uses STAR with UMITE (`star_umite`); gDNA uses BISCUIT with Methscan (`biscuit_methscan`).

`datadir` is the root containing the batch folders, and `outdir` is the output root. They default to `data/` and `results/` beside the `Snakefile`. All explicitly configured paths must be absolute and contain no whitespace.

### Shared references and existing indexes

The top-level `reference` settings apply to every batch. A batch may override individual reference fields through its own `reference` section. `reference.genes` is required for cDNA. `reference.genome` is required for gDNA and when building a STAR index; a cDNA-only batch using an existing STAR index can omit it. Use uncompressed genome and annotation files.

Supply either or both index settings to reuse existing indexes and skip the corresponding index builders:

| Setting | Value and required files |
| --- | --- |
| `reference.star_index` | Absolute directory containing `Genome`, `SA`, `SAindex`, and `genomeParameters.txt`. |
| `reference.biscuit_index` | Absolute prefix whose files have suffixes `.bis.amb`, `.bis.ann`, `.bis.pac`, `.dau.bwt`, `.dau.sa`, `.par.bwt`, and `.par.sa`. For `/refs/biscuit/genome`, one required file is `/refs/biscuit/genome.bis.amb`. |

Use indexes compatible with the aligner version and reference genome/annotation for the batch. Validation checks that the required files exist; it cannot verify that compatibility. If an index setting is omitted, the workflow builds that index under the corresponding batch's `results/<batch>/<modality>/reference/` directory. Multiple batches may share supplied indexes. STAR shared-memory cleanup waits for all alignments using the same index and does not make an existing batch's counting depend on a newly added batch.

### Batch overrides and resequencing

A batch can override `samples`, individual `reference` fields, processing option strings such as `star_align_args` or `methscan_filter_args`, and modality-specific metadata/FASTQ locations through `ilse_info`:

```yaml
batches:
  batch_01:
    modality: both
    samples: ['plate1_.*']
    ilse_info:
      cDNA:
        metadata:
          - /absolute/path/to/batch_01/cDNA/run1_meta.tsv
          - /absolute/path/to/batch_01/cDNA/run2_meta.tsv
        fastqdir: /absolute/path/to/batch_01/cDNA/fastq
      gDNA:
        samples: []
```

Each metadata or FASTQ directory override accepts a YAML list of absolute paths or a whitespace-separated string. Omitted fields use the batch's standard input layout. In batch mode, `ilse_info` must be inside a batch; a top-level `ilse_info` is rejected to avoid accidentally reusing inputs across batches.

The optional `samples` patterns must match the complete `SAMPLE_NAME` (or legacy `Sample Name`). For example, `['plate1_.*', 'plate2_.*']` selects those two plates. Omit the filter or use `samples: []` to include all samples. A modality's `ilse_info.<modality>.samples` overrides the batch filter; an explicit empty list includes all its samples. A filter selecting no samples with FASTQ IDs causes an error.

For resequencing, add the new metadata table directly to the same batch/modality folder and its FASTQs to that folder's `fastq/` directory. Exactly matching sample names automatically combine distinct FASTQ pairs across those tables. No merge-group configuration is needed. Add a new batch entry and input folder to process an independent batch; existing outputs are reused when their inputs and settings are unchanged.

Run the same Snakemake command after adding deliveries, keeping the existing `results/` and `.snakemake/` directories. Adding resequencing to an existing batch realigns affected cells using all their FASTQs and regenerates that batch's downstream aggregate outputs. Independent batches remain complete unless their own inputs, processing settings, or workflow rules change.

### Existing configurations without batches

The earlier configuration remains supported when `batches` is absent: use top-level `dataset`, `modality`, `reference`, and optional `ilse_info`, with inputs in `data/cDNA/` and `data/gDNA/`. Outputs retain their original paths, `results/cDNA/<dataset>.star_umite.h5ad` and `results/gDNA/<dataset>.biscuit_methscan.h5ad`. Existing input files do not need to move unless switching to the batch layout.

In this mode, if `modality` is omitted, `pipeline: star_umite` selects cDNA and `pipeline: biscuit_methscan` selects gDNA. Flat `ilse_info.metadata` and `ilse_info.fastqdir` fields remain valid for a single modality; `both` uses nested `ilse_info.cDNA` and `ilse_info.gDNA` sections. Set `pipeline: star_umite` for cDNA or both. Shared index settings are supported in this mode too.

## Methscan configuration

Methscan options for gDNA are configured through five option strings in `snakeconfig.yaml`. Set them at the top level for all batches or inside a batch to override its settings. Copy and edit the defaults below to override `defaultconfig.yaml`:

```yaml
methscan_prepare_args: >-
  --input-format biscuit_short
methscan_filter_args: >-
  --min-sites 50000
  --min-meth 40
  --max-meth 90
methscan_smooth_args: ""
methscan_scan_args: ""
methscan_matrix_args: ""
```

Each string is passed as command-line options to the corresponding `methscan prepare`, `filter`, `smooth`, `scan`, or `matrix` command. Empty strings use the tool defaults for the remaining options. The workflow controls threads and input/output paths; leave those out of these strings. Keep `--input-format biscuit_short` for the BISCUIT BED files generated by this pipeline. The preparation and filtering values shown above preserve the pipeline's existing defaults.

`methscan_matrix_args` cannot include `--sparse`, because H5AD conversion requires CSV output. `methscan_scan_args` cannot include `--write-header`, because the downstream matrix step requires a BED file without a header. To select cells from a file, include `--cell-names /absolute/path/to/cells.txt` in `methscan_filter_args`. This file must exist; the workflow tracks it as an input so changes trigger filtering again.

## Outputs

Each batch and modality writes its intermediate files, logs, and final matrix to its own directory beneath `results/` (or the configured `outdir`):

- cDNA: `results/<batch>/cDNA/<batch>.star_umite.h5ad`.
- gDNA: `results/<batch>/gDNA/<batch>.biscuit_methscan.h5ad`.

Generated reference artifacts are stored in `results/<batch>/<modality>/reference/`; supplied indexes are reused at their configured locations. A batch with both modalities produces two independent matrices for subsequent manual cell matching.

The `star_umite` branch pools exonic and intronic UMI evidence for each gene before UMI correction and deduplication (`--combine_unspliced`). The final H5AD uses `umicount/umite.U.tsv`, containing combined deduplicated UMI counts. Non-UMI internal fragment counts (`umite.R.tsv`) and duplicate counts (`umite.D.tsv`) remain separate TSV outputs.

### gDNA cell statistics plot

Batches with `modality: gDNA` or `modality: both` also produce `results/<batch>/gDNA/qc_plots/cell_stats.png` at 300 dpi. The scatter plot uses all cells from the unfiltered `methscan/compact_data/cell_stats.csv`: `global_meth_frac * 100` on the x-axis ("global DNA methylation %") and `n_obs` on the y-axis ("# of observed CpG sites"). To generate this plot and its prerequisites for `batch_01` with gDNA enabled:

```bash
snakemake --cores 1 --sdm conda batch_01_gDNA_plot_methscan_cell_stats
```

You can also plot an existing CSV directly in a Python environment with `pandas` and `matplotlib` (the workflow uses `envs/plotting.yaml`):

```bash
python scripts/plot_methscan_cell_stats.py results/batch_01/gDNA/methscan/compact_data/cell_stats.csv -o results/batch_01/gDNA/qc_plots/cell_stats.png
```

Omit `-o`/`--output` to save `cell_stats.png` in the current directory. Adjust the batch name and paths as needed. Without `batches`, the rule remains `gDNA_plot_methscan_cell_stats` and the plot remains under `results/gDNA/qc_plots/`.

# Running the pipeline

Change to the directory containing the script and activate the snakemake environment.

```bash
cd /path/to/pipeline/directory
conda activate snakemake
```

Run all batches and their selected modalities from `snakeconfig.yaml`:

```bash
snakemake --cores 24 --sdm conda --keep-incomplete
```

Change the number of cores to match your system's capabilities. For a dry run, add `-n`. Choose each batch's modalities in its configuration entry. For an existing configuration without `batches`, command-line selection remains available through `--config modality=cDNA`, `modality=gDNA`, or `modality=both`; cDNA and both require `pipeline=star_umite`.

# Known issues

If STAR is interrupted from outside (e.g. command line interruption or killed due to memory shortage), the loaded genome might stay in memory and become locked, preventing the next STAR instance to load any genome at all.
When this happens, use `STAR --genomeDir [/path/to/STAR/indices] --genomeLoad Remove`.
