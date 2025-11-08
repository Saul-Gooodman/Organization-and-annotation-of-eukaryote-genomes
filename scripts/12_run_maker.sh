#!/bin/bash
#SBATCH --job-name=maker_run
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00        
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err


set -euo pipefail


WORKDIR="/data/users/yliu2/Organization_and_annotation"
ANNDIR="$WORKDIR/gene_annotation"
LOGDIR="$WORKDIR/scripts/logs"
mkdir -p "$LOGDIR"

cd "$ANNDIR"


SIF="/data/courses/assembly-annotation-course/CDS_annotation/containers/MAKER_3.01.03.sif"


if command -v apptainer >/dev/null 2>&1; then
  RUN="apptainer exec --bind /data/users,/data/courses $SIF maker"
elif command -v maker >/dev/null 2>&1; then
  RUN="maker"
else
  echo "ERROR: neither apptainer nor maker found in PATH." >&2
  exit 2
fi

echo "Using MAKER command: $RUN"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"


$RUN -base assembly_p_ctg \
  maker_opts.ctl maker_bopts.ctl maker_exe.ctl | tee -a "$LOGDIR/maker_run.log"


gff3_merge -d assembly_p_ctg.maker.output/assembly_p_ctg_master_datastore_index.log \
  > assembly_p_ctg.all.maker.gff

fasta_merge -d assembly_p_ctg.maker.output/assembly_p_ctg_master_datastore_index.log

echo "=============================="
echo "MAKER annotation completed."
echo "GFF3:        $ANNDIR/assembly_p_ctg.all.maker.gff"
echo "Transcripts: $ANNDIR/assembly_p_ctg.all.maker.transcripts.fasta"
echo "Proteins:    $ANNDIR/assembly_p_ctg.all.maker.proteins.fasta"
echo "=============================="