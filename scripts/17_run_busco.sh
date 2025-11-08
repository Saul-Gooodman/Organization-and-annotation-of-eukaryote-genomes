#!/bin/bash
#SBATCH --job-name=busco_prot
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuwei.liu1@students.unibe.ch

set -euo pipefail
set -x

# Paths
WRK="/data/users/yliu2/Organization_and_annotation/gene_annotation/final"
LOGDIR="/data/users/yliu2/Organization_and_annotation/logs"
mkdir -p "$LOGDIR" "$WRK"
cd "$WRK"

# Inputs
PROT="$WRK/assembly_p_ctg.all.maker.proteins.renamed.filtered.fasta"
FALLBACK="$WRK/assembly_p_ctg.all.maker.proteins.renamed.fasta"
if [ ! -s "$PROT" ]; then
  echo "[WARN] filtered proteins missing or empty, fallback to renamed proteins"
  PROT="$FALLBACK"
fi
if [ ! -s "$PROT" ]; then
  echo "[FATAL] No usable protein FASTA: $PROT" >&2
  exit 2
fi

# Environment
module load BUSCO/5.4.2-foss-2021a
echo "[DEBUG] node=$(hostname)"
echo "[DEBUG] SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-unset}"
echo "[DEBUG] PROT=$PROT"
ls -lh "$PROT"
echo "[DEBUG] first header:"; grep -m1 '^>' "$PROT" || true
echo "[DEBUG] n_records=$(grep -c '^>' "$PROT" || echo 0)"
which busco || { echo "[FATAL] busco not found in PATH"; exit 3; }

# Make BUSCO DB dir and temp dir
DBDIR="$WRK/busco_db"
mkdir -p "$DBDIR"
export TMPDIR="$WRK/tmp_busco_${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR"

# Use absolute OUT dir
LINEAGE="brassicales_odb10"
OUTTAG="busco_proteins_${SLURM_JOB_ID:-manual}"
OUTDIR="$WRK/$OUTTAG"

# Unbuffer stdout/stderr for live logging
export PYTHONUNBUFFERED=1
if command -v stdbuf >/dev/null 2>&1; then
  STDBUF="stdbuf -oL -eL"
else
  STDBUF=""
fi

# Run
$STDBUF srun --cpu-bind=cores busco \
  -i "$PROT" \
  -l "$LINEAGE" \
  -o "$OUTDIR" \
  -m proteins \
  --cpu ${SLURM_CPUS_PER_TASK:-12} \
  --download_path "$DBDIR"

echo "[OK] BUSCO short summary:"
grep -H . "$OUTDIR"/run_*/short_summary*.txt || true

echo "[OK] BUSCO result dir: $OUTDIR"