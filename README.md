This repo contains a collection of pipelines to process scNMT sequencing data, intended for DKFZ A290 internal use.



# Important notes

- A lot of steps of the pipeline have issues with the Windows filesystem (e.g. the O drive). As a result, it is recommended to have **the pipeline folder, the output directory and the reference files all on a Linux volume** (e.g. scratch).
- Paths in the config file (see below) must be **absolute** paths and must not contain whitespaces.
- Both STAR and Salmon require a significant amount of RAM (up to 10 times the size of the unzipped genome) at indexing step. **Make sure your machine has enough RAM left.**

# Preparation
## Environment setup

The pipeline assumes avaiability of the following commands.

- `conda`
- `snakemake`

One way to get `conda` is to install Miniforge, following instructions on [their GitHub repo](https://github.com/conda-forge/miniforge), but feel free to choose any alternative.

Once `conda` is available, create a conda environment called "snakemake" and install necessary packages, which are `snakemake`, `pandas`, `xlrd` and `openpyxl`. To ensure portability, you can use the the provided environment file:

```bash
conda env create -n snakemake -f envs/snakemake.yaml
```

## Reference genome transcriptome and annotation

For most popular species, these files can be downloaded from [EMSELBL](https://www.ensembl.org/index.html).

For mouse GRCm39 release 115, use the following commands:

```bash
curl -O https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz
curl -O https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
curl -O https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz
```

Similarly, for human GRCh38 release 115, use:

```bash
curl -O https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
curl -O https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
curl -O https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
```

## Configuration

The config file `snakeconfig.yaml` is where you configure the pipeline parameters.

```yaml
dataset: [dataset_name]
outdir: [path_to_output_directory]

samples: []  # Include all samples; optionally list full-match regex patterns.

reference:
  genome: [path_to_genome_fasta]
  transcriptome: [path_to_transcriptome_fasta]
  genes: [path_to_annotation_gtf_gff]

ilse_info:
  metadata: [paths_to_metadata_files]
  fastqdir: [paths_to_fastq_directories]

pipeline: [pipeline_name]
```

- `[dataset_name]`: name of the dataset, used to name the output files.
- `[path_to_output_directory]`: the directory to store all outputs and log files of the pipeline.
- `[path_to_genome_fasta]`, `[path_to_transcriptome_fasta]` (optional), `[path_to_annotation_gtf_gff]`: path to the reference genome, reference transcriptome and gene annotation files. `[path_to_transcriptome_fasta]` is required only when `[pipeline_name]` is "salmon" (see below). 
- `[paths_to_metadata_files]`: metadata files (space_separated) downloaded from ILSe website, typically named `[ILSe ID]-result.xls`. The pipeline will identify samples by the "Sample Name" column in these metadata files, so if the same sample name appears in multiple metadata files, the pipeline will treat them as **one sample being sequenced multiple times**. Therefore, make sure that sample names are consistent across sequencing runs and remember to delete rows that do not belong to the dataset being mapped.
- `[paths_to_fastq_directories]`: directories (space_separated) to search for the FASTQ files, usually on the NGS drive. The pipeline will search for all FASTQ IDs specified in the "Unique ID / Lane" columns of the metadata files in these directories.
- `[pipeline_name]`: the pipeline used to map the sequences, possible options are "star_umite", "salmon" or "biscuit_methscan".
- `samples`: optional regular expression pattern(s) used to select samples from the "Sample Name" column. For example, `samples: ['plate1_.*', 'plate2_.*']` selects those two plates. Each pattern must match the complete sample name. Omit this key or use `samples: []` to include all samples. A filter that selects no samples with FASTQ IDs causes an explicit error.

> You can add any number of ILSe IDs to the `ilse_info` section, each with a metadata field and a fastq field. The pipeline will identify samples by the "Sample Name" column in the metadata file, so if the same sample name appears in two different metadata files, the pipeline will treat them as one sample being sequenced twice. Therefore, make sure that sample names are consistent across sequencing runs.

Finally, there is an experimental `tui.py` that provides a user interface to edit the file. This will require packages `textual` and `textual-dev` to run but will likely make it easier to create/edit config. Its sample-filter field accepts one regex per line; leave it blank to include all samples.

The `star_umite` branch pools exonic and intronic UMI evidence for each gene before UMI correction and deduplication (`--combine_unspliced`). The final H5AD uses `umicount/umite.U.tsv`, containing combined deduplicated UMI counts. Non-UMI internal fragment counts (`umite.R.tsv`) and duplicate counts (`umite.D.tsv`) remain separate TSV outputs.

# Running the pipeline

Change to the directory containing the script and activate the snakemake environment.

```bash
cd /path/to/pipeline/directory
conda activate snakemake
```

Use the following command to run the pipeline.

```bash
snakemake --cores 24 --sdm conda --keep-incomplete
```

Change the number of cores to match your system's capabilities. For a dry run, add the `-n` flag to the end of the command.


# Known issues

If STAR is interrupted from outside (e.g. command line interruption or killed due to memory shortage), the loaded genome might stay in memory and become locked, preventing the next STAR instance to load any genome at all.
When this happens, use `STAR --genomeDir [/path/to/STAR/indices] --genomeLoad Remove`.
