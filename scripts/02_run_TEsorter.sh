#!/bin/bash
#SBATCH --job-name=TEsorter_LTR
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err

set -euo pipefail

WORKDIR="/data/users/yliu2/Organization_and_annotation"
OUTDIR="$WORKDIR/results/EDTA_annotation"
RAW="$OUTDIR/assembly.p_ctg.fa.mod.EDTA.raw"
CONTAINER_TS="/data/courses/assembly-annotation-course/CDS_annotation/containers/TEsorter_1.3.0.sif"


ls -lh "$RAW/assembly.p_ctg.fa.mod.LTR.raw.fa"

cd "$OUTDIR"
apptainer exec \
  --bind "$WORKDIR","/data/courses/assembly-annotation-course" \
  "$CONTAINER_TS" \
  TEsorter "$RAW/assembly.p_ctg.fa.mod.LTR.raw.fa" -db rexdb-plant -p "$SLURM_CPUS_PER_TASK"