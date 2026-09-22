#!/bin/bash
#
# Submit the scNMT pipeline to LSF:
#
#     bsub < run_pipeline.sh
#
# The submitted job is only the coordinator: it stays alive for the whole run
# and submits one LSF job per rule instance, requesting the threads, mem_mb and
# walltime that rule declares. Arguments are forwarded to snakemake, so
# `./run_pipeline.sh -n` dry-runs the same configuration on the login node.
#
# Copying this pipeline elsewhere means editing the two #BSUB log paths and
# PIPELINE_DIR below; #BSUB lines are read by bsub before any shell variable
# exists, so they cannot be written in terms of PIPELINE_DIR.

#BSUB -J scnmt_coordinator
#BSUB -o /path/to/coordinator.log
#BSUB -e /path/to/coordinator.err
#BSUB -n 2
#BSUB -R "rusage[mem=8000]"
#BSUB -W 80:00

set -eo pipefail
ulimit -n 65536

module load Mamba/24.11.2-1

# Read each rule's mem_mb as the memory of the whole job. Without this the LSF
# plugin divides it by the job's threads and requests that much per core.
export SNAKEMAKE_LSF_MEMFMT=perjob

source "$HOME/.bashrc"
conda deactivate
conda activate /path/to/conda/env/with/snamekake

PIPELINE_DIR=/path/to/scNMT_pipeline
cd "$PIPELINE_DIR"

# The Snakefile reads snakeconfig.yaml relative to the working directory.
if [ ! -f snakeconfig.yaml ]; then
    echo "No snakeconfig.yaml in $PIPELINE_DIR: copy snakeconfig.example.yaml and edit it." >&2
    exit 1
fi

# Create the conda environments once, before any job is submitted, so that
# hundreds of jobs do not race to build the same environment:
#     snakemake --sdm conda --conda-create-envs-only --cores 1
snakemake -p \
  --executor lsf \
  --default-resources lsf_project=project_name mem_mb=4000 walltime=60 \
  --jobs 100 \
  --cores 16 \
  --local-cores 2 \
  --sdm conda \
  --latency-wait 60 \
  --retries 2 \
  --rerun-incomplete \
  --keep-going \
  -s "$PIPELINE_DIR/Snakefile" \
  "$@"

# --jobs limits how many LSF jobs run at once; --cores is the ceiling on a
# single job's threads and must stay at or above the largest threads: in the
# workflow (16, generate_star_index), or Snakemake silently requests fewer.
