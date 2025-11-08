#!/bin/bash
#SBATCH --job-name=maker_resume
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=4-00:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err
#SBATCH --mail-user=yuwei.liu1@students.unibe.ch
#SBATCH --mail-type=END,FAIL

set -euo pipefail

WORKDIR="/data/users/yliu2/Organization_and_annotation"
ANNDIR="$WORKDIR/gene_annotation"
LOGDIR="$WORKDIR/scripts/logs"
mkdir -p "$LOGDIR"

cd "$ANNDIR"

SIF="/data/courses/assembly-annotation-course/CDS_annotation/containers/MAKER_3.01.03.sif"

if command -v apptainer >/dev/null 2>&1; then
  RUN="apptainer exec --bind /data/users,/data/courses $SIF maker"
  EXEC="apptainer exec --bind /data/users,/data/courses $SIF"
elif command -v maker >/dev/null 2>&1; then
  RUN="maker"
  EXEC=""
else
  echo "ERROR: neither apptainer nor maker found in PATH." >&2
  exit 2
fi

echo "Using MAKER command: $RUN"
echo "CPUs per task: ${SLURM_CPUS_PER_TASK:-12}"
echo "Resuming from existing MAKER output..."

# ========================
# Check unfinished contigs before resuming
# ========================
INDEX_FILE="assembly_p_ctg.maker.output/assembly_p_ctg_master_datastore_index.log"

if [[ -f "$INDEX_FILE" ]]; then
  total=$(awk 'NR>1 {count++} END {print count+0}' "$INDEX_FILE")
  done=$(grep -c 'FINISHED' "$INDEX_FILE" || true)
  failed=$(grep -c 'FAILED' "$INDEX_FILE" || true)
  running=$(grep -c 'RUNNING' "$INDEX_FILE" || true)
  pending=$((total - done))
  echo "-------------------------------------------"
  echo "MAKER progress before resuming:"
  echo "Total contigs:   $total"
  echo "Finished:        $done"
  echo "Failed:          $failed"
  echo "Running:         $running"
  echo "Pending (to run):$pending"
  echo "-------------------------------------------"
else
  echo "Warning: master_datastore_index.log not found. Assuming full rerun."
fi

# ========================
# Resume MAKER
# ========================

$RUN -base assembly_p_ctg -fix_nucleotides -cpus ${SLURM_CPUS_PER_TASK:-12} \
  maker_opts.ctl maker_bopts.ctl maker_exe.ctl | tee -a "$LOGDIR/maker_resume.log"

# ========================
# Merge results after run
# ========================

if [ -n "${EXEC}" ]; then
  ${EXEC} gff3_merge -d assembly_p_ctg.maker.output/assembly_p_ctg_master_datastore_index.log \
    > assembly_p_ctg.all.maker.gff
  ${EXEC} fasta_merge -d assembly_p_ctg.maker.output/assembly_p_ctg_master_datastore_index.log
else
  gff3_merge -d assembly_p_ctg.maker.output/assembly_p_ctg_master_datastore_index.log \
    > assembly_p_ctg.all.maker.gff
  fasta_merge -d assembly_p_ctg.maker.output/assembly_p_ctg_master_datastore_index.log
fi

echo "=============================="
echo "MAKER resumed annotation completed."
echo "GFF3:        $ANNDIR/assembly_p_ctg.all.maker.gff"
echo "Transcripts: $ANNDIR/assembly_p_ctg.all.maker.transcripts.fasta"
echo "Proteins:    $ANNDIR/assembly_p_ctg.all.maker.proteins.fasta"
echo "=============================="