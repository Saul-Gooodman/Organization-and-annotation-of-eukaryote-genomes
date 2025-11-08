#!/usr/bin/env bash
#SBATCH --job-name=genespace_make_rds
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err

set -euo pipefail

BASE="/data/users/yliu2/Organization_and_annotation"
WD="$BASE/genespace"
SIF="/data/courses/assembly-annotation-course/CDS_annotation/containers/genespace_latest.sif"

echo "[INFO] Node=$(hostname)"
echo "[INFO] Start=$(date -Iseconds)"
echo "[INFO] WD=$WD"
mkdir -p "$WD/pangenes"

# Run a tiny R snippet inside the container to extract and save the pangenome matrix.
apptainer exec \
  --bind "$BASE":/work \
  "$SIF" bash -lc '
    Rscript - << "RS"
    suppressPackageStartupMessages(library(GENESPACE))

    wd <- "/work/genespace"
    # init_genespace with correct MCScanX path (per instructor update)
    gpar <- init_genespace(wd = wd, path2mcscanx = "/usr/local/bin")

    # This will re-use existing OrthoFinder/DIAMOND results; no heavy recompute.
    out <- run_genespace(gpar, overwrite = FALSE)

    # choose reference genome
    gids <- out$params$genomeIDs
    refGenome <- if ("TAIR10" %in% gids) "TAIR10" else gids[1]

    # build pangene matrix and save
    pan <- query_pangenes(
      out,
      bed = NULL,
      refGenome = refGenome,
      transform = TRUE,
      showArrayMem = TRUE,
      showNSOrtho = TRUE,
      maxMem2Show = Inf
    )

    saveRDS(pan, file = file.path(wd, "pangenes", "pangenome_matrix.rds"))
    cat("[OK] Saved:", file.path(wd, "pangenes", "pangenome_matrix.rds"), "\n")
RS
  '

echo "[INFO] End=$(date -Iseconds)"
echo "[DONE] Pangenome RDS created."