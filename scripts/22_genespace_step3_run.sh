#!/usr/bin/env bash
#SBATCH --job-name=genespace_run
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuwei.liu1@students.unibe.ch

set -euo pipefail

echo "[INFO] Node=$(hostname)"
echo "[INFO] Start=$(date -Iseconds)"

BASE="/data/users/yliu2/Organization_and_annotation"
WD="$BASE/genespace"
SCR="$BASE/tmp/${SLURM_JOB_ID}"
COURSEDIR="/data/courses/assembly-annotation-course"
LOGDIR="$BASE/logs"
SCRIPTDIR="$BASE/scripts"
mkdir -p "$SCR" "$LOGDIR" "$SCRIPTDIR"

RS="$SCRIPTDIR/genespace.R"
cat > "$RS" << 'RSCRIPT'
suppressPackageStartupMessages({
  library(GENESPACE)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Please provide working directory (wd) as arg[1].")
wd <- args[1]
pepDir <- file.path(wd, "peptide")
bedDir <- file.path(wd, "bed")

pepFiles <- list.files(pepDir, pattern="\\.fa(sta)?$", full.names=TRUE)
if (length(pepFiles) < 2) stop("Need >=2 genomes (peptide FASTA) to run GENESPACE.")

genomeIDs <- sub("\\.fa(sta)?$", "", basename(pepFiles))
df <- data.frame(
  genome = genomeIDs,
  pep    = file.path(pepDir, paste0(genomeIDs, ".fa")),
  bed    = file.path(bedDir,  paste0(genomeIDs, ".bed")),
  stringsAsFactors = FALSE
)

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = "/usr/local/bin/MCScanX"
)

out <- run_genespace(
  gpar,
  overwrite = TRUE
)

refGenome <- if ("TAIR10" %in% df$genome) "TAIR10" else df$genome[1]

pan <- query_pangenes(
  out,
  bed = NULL,
  refGenome = refGenome,
  transform = TRUE,
  showArrayMem = TRUE,
  showNSOrtho = TRUE,
  maxMem2Show = Inf
)
saveRDS(pan, file = file.path(wd, "pangenome_matrix.rds"))
cat("[OK] GENESPACE finished.\n")
cat("[OK] Pangenome matrix: ", file.path(wd, "pangenome_matrix.rds"), "\n", sep="")
RSCRIPT

SIF="$COURSEDIR/CDS_annotation/containers/genespace_latest.sif"
if [ ! -s "$SIF" ]; then
  echo "[FATAL] genespace_latest.sif not found at $SIF" >&2
  exit 2
fi

export LC_ALL=C.UTF-8
export LANG=C.UTF-8

echo "[INFO] Checking required BED/FASTA:"
for f in "$WD/bed/TAIR10.bed" "$WD/peptide/TAIR10.fa" \
         "$WD/bed/Yuwei.bed"  "$WD/peptide/Yuwei.fa"; do
  if [ ! -s "$f" ]; then
    echo "[FATAL] Missing file: $f" >&2
    exit 3
  fi
  echo "  OK  $f"
done

echo "[INFO] Running GENESPACE via Apptainer..."
apptainer exec \
  --bind "$COURSEDIR" \
  --bind "$BASE" \
  --bind "$SCR:/temp" \
  "$SIF" \
  Rscript "$RS" "$WD"

echo "[INFO] End=$(date -Iseconds)"
echo "[DONE] GENESPACE run completed successfully."