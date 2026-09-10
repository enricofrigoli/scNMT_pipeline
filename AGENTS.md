# AGENTS.md

Repository guide for agentic coding in this project.

## Project Shape
- Snakemake pipeline for scNMT processing.
- Main workflow entrypoint: `Snakefile`.
- Default pipeline config: `defaultconfig.yaml`.
- User config: `snakeconfig.yaml`; copyable template: `snakeconfig.example.yaml`.
- Reference docs: `README.md`.

## Build / Lint / Test
- Create the pipeline env: `conda env create -n snakemake -f envs/snakemake.yaml`.
- Activate it before pipeline work: `conda activate snakemake`.
- Require Snakemake >=9.19 for dynamically named batch modules.
- Run the pipeline: `snakemake --cores 24 --sdm conda --keep-incomplete`.
- Dry run the workflow: `snakemake -n`.
- Run a smaller workflow slice: `snakemake --cores 1 --until <rule>`.
- Validate a config-driven run without executing jobs: `snakemake -n --configfile snakeconfig.yaml`.

## Single Test / Narrow Check
- Run regression tests with `python -m unittest discover -s tests -v`.
- For a narrow workflow check, run one rule path with `snakemake --cores 1 --until <rule>`.
- For config validation, prefer `snakemake -n` after editing `snakeconfig.yaml`.

## Config Expectations
- `snakeconfig.yaml` is the editable config file.
- It should follow the nested structure described in `README.md`.
- The primary configuration uses top-level `batches`, shared `reference`, optional `datadir`/`outdir`, and shared processing options.
- Each `batches.<batch>` requires `modality` and may override `samples`, individual `reference` fields, `ilse_info`, and processing option strings.
- Batch IDs allow only letters, digits, and underscores; do not introduce cross-batch merge groups.
- Without `batches`, preserve legacy `dataset`, `modality`, `reference`, `ilse_info`, `samples`, and `pipeline` behavior and paths.
- `reference.genes` is required for cDNA. `reference.genome` is required for gDNA and STAR index building; cDNA with a supplied STAR index can omit it. Use uncompressed references.
- Optional `reference.star_index` is an absolute directory containing `Genome`, `SA`, `SAindex`, and `genomeParameters.txt`.
- Optional `reference.biscuit_index` is an absolute file prefix with `.bis.amb`, `.bis.ann`, `.bis.pac`, `.dau.bwt`, `.dau.sa`, `.par.bwt`, and `.par.sa` files.
- Supplied indexes skip builders and can be shared between batches.
- Generate the STAR index once per genome in `star_index` beside `reference.genome`, from the top-level `Snakefile`; never declare it inside the per-batch module.
- Derive `--sjdbOverhang` from the sampled FASTQ read length minus one, taking the longest read across every batch sharing a generated index.
- Generate an omitted `reference.biscuit_index` per batch under its result directory.
- `datadir` and `outdir` default to `data/` and `results/` beside `Snakefile`.
- Batch `modality` selects `cDNA`, `gDNA`, or `both` and is mandatory; only legacy configurations infer omitted modality from the pipeline selector.
- Discover metadata tables in `data/<batch>/<modality>/` and reads under its `fastq/` folder; legacy configurations retain `data/<modality>/`.
- Accept current facility metadata columns `SAMPLE_NAME`, `FASTQ_FILE`, and `READ`; retain support for legacy `Sample Name` and `Unique ID / Lane`.
- Normalize metadata only in memory. Derive FASTQ IDs by stripping `_R1.fastq.gz` or `_R2.fastq.gz` and require `READ` to match the filename's mate (`1` or `2`).
- Group resequencing by exact `SAMPLE_NAME` (or legacy `Sample Name`) within each batch and modality; mate rows for an ID represent one read pair.
- Skip blank `SAMPLE_NAME` rows only for `Undetermined*` FASTQs; reject other rows with missing sample names.
- Optional batch `ilse_info.cDNA` and `ilse_info.gDNA` sections override `metadata`, `fastqdir`, and `samples` independently.
- Reject top-level `ilse_info` in batch mode to prevent accidental reuse of inputs across batches.
- Metadata and FASTQ directory overrides accept whitespace-separated strings or YAML lists.
- Legacy flat `ilse_info` inputs are valid only for a single modality.
- All explicit config paths must be absolute and contain no whitespace.
- Validate input paths for selected modalities only; both modalities must have data when `modality: both`.
- Never match or intersect cells or FASTQ IDs across batches or modalities; outputs are independent under `results/<batch>/<modality>/`.
- Final batch matrices are `<batch>.star_umite.h5ad` and `<batch>.biscuit_methscan.h5ad`; retain legacy output names when `batches` is absent.
- Deduplicate repeated FASTQ IDs within a cell.
- Configure Methscan through string fields `methscan_prepare_args`, `methscan_filter_args`, `methscan_smooth_args`, `methscan_scan_args`, and `methscan_matrix_args`.
- Configure BISCUIT pileup through `biscuit_pileup_args`, defaulting to `-p` to keep improper pairs; the workflow always adds `-N` for NOMe-seq.
- Install UMITE from `github.com/enricofrigoli/umite` at the commit pinned in `envs/umite.yaml`, not from PyPI; that fork provides CIGAR-aware lookup and `umicount --stranded`.
- Configure UMITE through the verbatim strings `umiextract_args` and `umicount_args`; reject options the workflow supplies itself (inputs, outputs, logs, cores, GTF options, `--combine_unspliced`).
- Resolve `umicount --stranded` (`no`, `umi`, `yes`, `reverse`; default `umi`) during config normalization, and pass the same mode to the GTF dump and to counting; the mode names the dump because UMITE refuses a mismatched one.
- Keep Methscan input/output paths and threads controlled by the workflow; preparation uses `--input-format biscuit_short` for BISCUIT BED input.
- Legacy dataset names, sample names, and FASTQ IDs use only letters, digits, underscores, dots, and hyphens; batch IDs are stricter (letters, digits, underscores only). Sample names containing spaces remain unsupported.

## Code Style
- Use 4-space indentation.
- Prefer explicit imports over wildcard imports.
- Group imports as stdlib, third-party, local.
- Use `yaml.safe_load` / `yaml.safe_dump`; avoid unsafe YAML APIs.
- Keep config parsing and validation close to the workflow rules in `Snakefile`.
- Use f-strings for readable error messages.
- Keep functions small and single-purpose.
- Prefer clear names over clever abbreviations.
- Use `PascalCase` for classes, `snake_case` for functions and variables, `UPPER_SNAKE_CASE` for constants.
- Annotate public helper functions when practical.

## Imports
- Standard library first, then third-party packages, then local modules.
- Avoid re-importing the same module under multiple names.
- Do not add unused imports.
- Prefer `pathlib.Path` or `os.path` consistently within a file.

## Formatting
- Keep lines reasonably short and readable.
- Prefer straightforward control flow over deeply nested conditionals.
- Use blank lines to separate logical sections.
- Preserve YAML key order when writing configs.
- Avoid unnecessary comments; add comments only where the logic is not obvious.

## Types / Data Handling
- Treat config values as strings until validation proves otherwise.
- Convert integers only after successful validation.
- Represent multi-path YAML fields as whitespace-separated strings or lists of paths.
- When loading config, handle missing sections defensively.
- Preserve unknown top-level keys when editing a config unless you are intentionally removing them.

## Naming Conventions
- Match existing pipeline terminology from the README and Snakefile.
- Use `batches`, `dataset`, `datadir`, `outdir`, `modality`, `reference`, `ilse_info`, `samples`, and `pipeline` consistently.
- Use descriptive names such as `load_config`, `validate_config`, and `build_config`.

## Error Handling
- Raise `ValueError` or a custom exception for invalid user input.
- Surface configuration validation errors instead of silently fixing them.
- Fail fast on missing mandatory config sections.
- Check absolute paths and existence where the Snakefile does so.
- Do not swallow file-write errors; report them clearly.

## Workflow Rules
- Keep `Snakefile` semantics in mind when editing config logic.
- For selected cDNA, require `pipeline: star_umite`; gDNA always uses `biscuit_methscan`.
- Require `reference.genes` only when cDNA is selected.
- Keep `snakeconfig.yaml` overwriting intentional and explicit.
- Coordinate STAR shared-memory cleanup once per index after all alignments using it; do not make counting depend on cleanup. Appending a batch must not force previous batches to recount.

## Repository Notes
- The repo currently has no `.cursor/rules/` files.
- The repo currently has no `.cursorrules` file.
- The repo currently has no `.github/copilot-instructions.md` file.
- If those files are added later, follow them in addition to this guide.

## Validation Checklist
- Verify `snakeconfig.yaml` remains valid YAML after any edit.
- Validate batch IDs and required batch modalities, or verify `dataset` is non-empty for legacy configurations.
- Verify explicit `datadir`, `outdir`, and reference paths are absolute.
- Verify selected modality metadata and FASTQ directory overrides contain absolute paths.
- Verify `pipeline` is one of the supported values before launching Snakemake.
- Prefer `snakemake -n --configfile snakeconfig.yaml` after config edits.

## Config Writing Rules
- Load an existing `snakeconfig.yaml` and preserve unknown top-level keys when reasonable.
- Write nested `batches`, `reference`, and batch `ilse_info` sections explicitly, not as flattened placeholders.
- Preserve supported multi-path values as whitespace-separated strings or YAML lists.
- Remove optional keys only when intentionally clearing them.
- Validate before writing configuration files.

## Agent Behavior
- Make focused changes.
- Do not touch unrelated files.
- Preserve existing user changes.
- Prefer minimal, reviewable diffs.
- When in doubt, align with `README.md`, `Snakefile`, and `defaultconfig.yaml`.
