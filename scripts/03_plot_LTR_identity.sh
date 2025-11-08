#!/bin/bash
#SBATCH --job-name=Plot_LTR_identity
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err
set -euo pipefail


module load R/4.3.2-foss-2021a
export R_LIBS_USER="/data/users/yliu2/Rlibs/4.3"
mkdir -p "$R_LIBS_USER"

Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE);
            repos <- c(CRAN="https://cloud.r-project.org");
            pkgs <- c("tidyverse","data.table","cowplot","optparse");
            inst <- rownames(installed.packages(lib.loc=Sys.getenv("R_LIBS_USER")));
            need <- setdiff(pkgs, inst);
            if(length(need)) install.packages(need, repos=repos, lib=Sys.getenv("R_LIBS_USER"), Ncpus=2)'


WORKDIR="/data/users/yliu2/Organization_and_annotation"
OUTDIR="$WORKDIR/results/EDTA_annotation"
RAW="$OUTDIR/assembly.p_ctg.fa.mod.EDTA.raw"
SCRIPT="/data/courses/assembly-annotation-course/CDS_annotation/scripts/02-full_length_LTRs_identity.R"


cd "$OUTDIR"
ln -sf "$RAW/assembly.p_ctg.fa.mod.LTR.intact.raw.gff3" genomic.fna.mod.LTR.intact.raw.gff3
ln -sf "assembly.p_ctg.fa.mod.LTR.raw.fa.rexdb-plant.cls.tsv" genomic.fna.mod.LTR.raw.fa.rexdb-plant.cls.tsv


mkdir -p plots


Rscript "$SCRIPT"


echo "Plot generated files:"
ls -lh "$OUTDIR"/plots/*pdf 2>/dev/null || ls -lh "$OUTDIR"/*pdf || true