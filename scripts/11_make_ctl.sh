#!/bin/bash
#SBATCH --job-name=maker_ctl
#SBATCH --partition=pshort_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/scripts/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/scripts/logs/%x_%j.err

set -euo pipefail


WORKDIR="/data/users/yliu2/Organization_and_annotation/gene_annotation"
mkdir -p "$WORKDIR"
cd "$WORKDIR"


module load Apptainer/1.2.4 2>/dev/null || module load Singularity/3.8.5 2>/dev/null || echo "No Apptainer/Singularity module found"


if command -v apptainer &>/dev/null; then
  echo "Running MAKER via Apptainer..."
  apptainer exec --bind "$WORKDIR" \
  /data/courses/assembly-annotation-course/CDS_annotation/containers/MAKER_3.01.03.sif \
  maker -CTL
elif command -v singularity &>/dev/null; then
  echo "Running MAKER via Singularity..."
  singularity exec --bind "$WORKDIR" \
  /data/courses/assembly-annotation-course/CDS_annotation/containers/MAKER_3.01.03.sif \
  maker -CTL
else
  echo "ERROR: Neither apptainer nor singularity found on this node." >&2
  exit 1
fi