#!/bin/bash
#SBATCH --job-name=EDTA_annotation
#SBATCH --partition=pitz_el8
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err

WORKDIR="/data/users/yliu2/Organization_and_annotation"
GENOME="/data/users/yliu2/assembly_annotation_course/assemblies/hifiasm/assembly.p_ctg.fa"
CONTAINER="/data/courses/assembly-annotation-course/CDS_annotation/containers/EDTA2.2.sif"
CDS="/data/courses/assembly-annotation-course/CDS_annotation/data/TAIR10_cds_20110103_representative_gene_model_updated"
OUTDIR="$WORKDIR/results/EDTA_annotation"

mkdir -p "$OUTDIR"
cd "$OUTDIR"


apptainer exec \
  --bind "$WORKDIR","/data/users/yliu2/assembly_annotation_course","/data/courses/assembly-annotation-course" \
  "$CONTAINER" \
  EDTA.pl \
  --genome "$GENOME" \
  --species others \
  --step all \
  --sensitive 1 \
  --cds "$CDS" \
  --anno 1 \
  --threads $SLURM_CPUS_PER_TASK \
  --overwrite 1