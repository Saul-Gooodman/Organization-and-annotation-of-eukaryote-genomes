#!/bin/bash
#SBATCH --job-name=clade_summary
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.out
#SBATCH --error=/data/users/yliu2/Organization_and_annotation/logs/%x_%j.err
set -euo pipefail

module load R/4.3.2-foss-2021a
export R_LIBS_USER="/data/users/yliu2/Rlibs/4.3"
mkdir -p "$R_LIBS_USER"

Rscript -e 'repos=c(CRAN="https://cloud.r-project.org");
            pkgs=c("data.table","dplyr","stringr","ggplot2","tidytext");
            inst=rownames(installed.packages(lib.loc=Sys.getenv("R_LIBS_USER")));
            need=setdiff(pkgs, inst);
            if(length(need)) install.packages(need, repos=repos, lib=Sys.getenv("R_LIBS_USER"), Ncpus=2)'

WORKDIR="/data/users/yliu2/Organization_and_annotation"
OUTDIR="$WORKDIR/results/EDTA_annotation"
PLOTDIR="$OUTDIR/summary"
mkdir -p "$PLOTDIR"

CLS1="$OUTDIR/assembly.p_ctg.fa.mod.LTR.raw.fa.rexdb-plant.cls.tsv"
CLS_C="$OUTDIR/Copia_sequences.fa.rexdb-plant.cls.tsv"
CLS_G="$OUTDIR/Gypsy_sequences.fa.rexdb-plant.cls.tsv"

CLS_LIST=()
[[ -f "$CLS1" ]] && CLS_LIST+=("$CLS1")
[[ -f "$CLS_C" ]] && CLS_LIST+=("$CLS_C")
[[ -f "$CLS_G" ]] && CLS_LIST+=("$CLS_G")

if [[ ${#CLS_LIST[@]} -eq 0 ]]; then
  echo "ERROR: No TEsorter classification *.cls.tsv found under $OUTDIR" >&2
  exit 2
fi

echo "Using cls files:"
printf '  %s\n' "${CLS_LIST[@]}"

Rscript - "$PLOTDIR" "${CLS_LIST[@]}" <<'RSCRIPT'
args <- commandArgs(TRUE)
outdir <- args[1]
files  <- args[-1]

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(stringr); library(ggplot2); library(tidytext)
})

read_cls <- function(f){
  dt <- data.table::fread(f)
  data.table::setnames(dt, tolower(names(dt)))
  if(!"superfamily" %in% names(dt)) stop("No 'superfamily' column in: ", f)
  if(!"order" %in% names(dt)) dt[, order := NA_character_]
  if(!"clade" %in% names(dt)) dt[, clade := "UNCLASSIFIED"]
  dt[, source_file := basename(f)]
  dt
}

all <- data.table::rbindlist(lapply(files, read_cls), use.names = TRUE, fill = TRUE)
all$superfamily <- toupper(all$superfamily)

sub <- all %>% 
  dplyr::filter(grepl("^LTR", toupper(order)) | superfamily %in% c("COPIA","GYPSY")) %>%
  dplyr::mutate(
    superfamily = dplyr::case_when(
      grepl("COPIA", superfamily, ignore.case = TRUE) ~ "COPIA",
      grepl("GYPSY", superfamily, ignore.case = TRUE) ~ "GYPSY",
      TRUE ~ superfamily
    ),
    clade = dplyr::if_else(is.na(clade) | clade=="", "UNCLASSIFIED", clade)
  )

sum_by_clade <- sub %>% 
  dplyr::group_by(superfamily, clade) %>% 
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(superfamily, dplyr::desc(n))

sum_by_sf <- sub %>% dplyr::count(superfamily, name = "total")
sum_by_clade <- sum_by_clade %>% 
  dplyr::left_join(sum_by_sf, by="superfamily") %>%
  dplyr::mutate(pct = 100 * n / total)

f_csv  <- file.path(outdir, "TEsorter_Copia_Gypsy_clade_summary.csv")
f_top  <- file.path(outdir, "TEsorter_Copia_Gypsy_clade_top12.csv")
f_plot1 <- file.path(outdir, "TEsorter_Copia_Gypsy_clade_bar.pdf")
f_plot2 <- file.path(outdir, "TEsorter_Copia_Gypsy_clade_bar.png")

data.table::fwrite(sum_by_clade, f_csv)

top12 <- sum_by_clade %>% 
  dplyr::group_by(superfamily) %>% 
  dplyr::slice_max(n, n = 12, with_ties = FALSE) %>%
  dplyr::ungroup()
data.table::fwrite(top12, f_top)

p <- ggplot2::ggplot(top12, ggplot2::aes(
        x = n,
        y = tidytext::reorder_within(clade, n, superfamily),
        fill = superfamily)) +
  ggplot2::geom_col() +
  tidytext::scale_y_reordered() +
  ggplot2::facet_wrap(~ superfamily, scales = "free_y", ncol = 1) +
  ggplot2::labs(x = "Count (classified sequences)", y = "Clade",
       title = "Copia & Gypsy clade counts (TEsorter)") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "none")

ggplot2::ggsave(f_plot1, p, width = 8, height = 10)
ggplot2::ggsave(f_plot2, p, width = 8, height = 10, dpi = 200)

message("Saved: ", f_csv)
message("Saved: ", f_top)
message("Saved: ", f_plot1)
message("Saved: ", f_plot2)
RSCRIPT

echo "Done."