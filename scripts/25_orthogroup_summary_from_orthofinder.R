suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})


ORTHO_FILE <- "/data/users/yliu2/Organization_and_annotation/genespace/orthofinder/Results_Nov08/Orthogroups/Orthogroups.tsv"
OUT_SUMMARY <- "/data/users/yliu2/Organization_and_annotation/genespace/orthogroup_summary_from_orthofinder.csv"

cat("[INFO] Reading Orthogroups from:\n  ", ORTHO_FILE, "\n")

og <- fread(ORTHO_FILE, sep = "\t", header = TRUE)

cat("[INFO] Dimensions of Orthogroups table: ", paste(dim(og), collapse = " x "), "\n")
cat("[INFO] Column names:\n")
print(colnames(og))


if (!all(c("Orthogroup", "TAIR10", "Yuwei") %in% colnames(og))) {
  stop("Cannot find columns 'Orthogroup', 'TAIR10', 'Yuwei' in Orthogroups.tsv")
}

og2 <- og[, .(Orthogroup, TAIR10, Yuwei)]




present_TAIR10 <- og2$TAIR10 != "" & !is.na(og2$TAIR10)
present_Yuwei  <- og2$Yuwei  != "" & !is.na(og2$Yuwei)

core_raw      <- og2$Orthogroup[present_TAIR10 & present_Yuwei]
tair_only_raw <- og2$Orthogroup[present_TAIR10 & !present_Yuwei]
yuwei_only_raw<- og2$Orthogroup[!present_TAIR10 & present_Yuwei]

cat("\n[RAW ORTHOGROUP COUNTS]\n")
cat("Core orthogroups (TAIR10 & Yuwei): ", length(core_raw), "\n")
cat("TAIR10-only orthogroups:          ", length(tair_only_raw), "\n")
cat("Yuwei-only orthogroups:           ", length(yuwei_only_raw), "\n")



split_genes <- function(x) {
  if (is.na(x) || x == "") return(character(0))

  y <- unlist(strsplit(x, ","))
  y <- str_trim(y)
  y[y != ""]
}


expand_genes <- function(df) {
  res_list <- list()
  k <- 1
  for (i in seq_len(nrow(df))) {
    og_id <- df$Orthogroup[i]
    # TAIR10 genes
    g_t <- split_genes(df$TAIR10[i])
    if (length(g_t) > 0) {
      res_list[[k]] <- data.table(Orthogroup = og_id,
                                  Genome = "TAIR10",
                                  Gene   = g_t)
      k <- k + 1
    }

    g_y <- split_genes(df$Yuwei[i])
    if (length(g_y) > 0) {
      res_list[[k]] <- data.table(Orthogroup = og_id,
                                  Genome = "Yuwei",
                                  Gene   = g_y)
      k <- k + 1
    }
  }
  if (length(res_list) == 0) {
    return(data.table(Orthogroup=character(0),
                      Genome=character(0),
                      Gene=character(0)))
  }
  rbindlist(res_list)
}

og_gene_long <- expand_genes(og2)


core_genes_raw <- og_gene_long[Orthogroup %in% core_raw]
tair_genes_raw <- og_gene_long[Genome == "TAIR10" & Orthogroup %in% tair_only_raw]
yuwei_genes_raw<- og_gene_long[Genome == "Yuwei"  & Orthogroup %in% yuwei_only_raw]

cat("\n[RAW GENE COUNTS]\n")
cat("Core genes (both genomes): ", nrow(core_genes_raw), "\n")
cat("TAIR10-only genes:         ", nrow(tair_genes_raw), "\n")
cat("Yuwei-only genes:          ", nrow(yuwei_genes_raw), "\n")


gene_counts <- og_gene_long[, .N, by = .(Orthogroup, Genome)]
gene_counts_wide <- dcast(gene_counts, Orthogroup ~ Genome, value.var = "N", fill = 0)


flt_ogs <- gene_counts_wide[
  TAIR10 >= 1 & TAIR10 <= 5 &
    Yuwei  >= 1 & Yuwei  <= 5,
  Orthogroup
]

cat("\n[FILTER] Keep orthogroups with 1–5 genes in TAIR10 and 1–5 genes in Yuwei\n")
cat("[FILTER] Kept ", length(flt_ogs), " orthogroups out of ", nrow(gene_counts_wide), " total (",
    round(100 * length(flt_ogs)/nrow(gene_counts_wide), 2), "%)\n", sep = "")


present_TAIR10_f <- present_TAIR10 & og2$Orthogroup %in% flt_ogs
present_Yuwei_f  <- present_Yuwei  & og2$Orthogroup %in% flt_ogs

core_f      <- og2$Orthogroup[present_TAIR10_f & present_Yuwei_f]
tair_only_f <- og2$Orthogroup[present_TAIR10_f & !present_Yuwei_f]
yuwei_only_f<- og2$Orthogroup[!present_TAIR10_f & present_Yuwei_f]

cat("\n[FILTERED ORTHOGROUP COUNTS]\n")
cat("Core orthogroups (TAIR10 & Yuwei): ", length(core_f), "\n")
cat("TAIR10-only orthogroups:          ", length(tair_only_f), "\n")
cat("Yuwei-only orthogroups:           ", length(yuwei_only_f), "\n")


core_genes_f <- og_gene_long[Orthogroup %in% core_f]
tair_genes_f <- og_gene_long[Genome == "TAIR10" & Orthogroup %in% tair_only_f]
yuwei_genes_f<- og_gene_long[Genome == "Yuwei"  & Orthogroup %in% yuwei_only_f]

cat("\n[FILTERED GENE COUNTS]\n")
cat("Core genes (both genomes): ", nrow(core_genes_f), "\n")
cat("TAIR10-only genes:         ", nrow(tair_genes_f), "\n")
cat("Yuwei-only genes:          ", nrow(yuwei_genes_f), "\n")



summary_df <- data.frame(
  Category = c("Core", "TAIR10_unique", "Kas1_unique"),
  Orthogroups_raw      = c(length(core_raw), length(tair_only_raw), length(yuwei_only_raw)),
  Genes_raw_TAIR10     = c(
    nrow(core_genes_raw[Genome == "TAIR10"]),
    nrow(tair_genes_raw),
    0
  ),
  Genes_raw_Kas1       = c(
    nrow(core_genes_raw[Genome == "Yuwei"]),
    0,
    nrow(yuwei_genes_raw)
  ),
  Orthogroups_filtered = c(length(core_f), length(tair_only_f), length(yuwei_only_f)),
  Genes_filtered_TAIR10= c(
    nrow(core_genes_f[Genome == "TAIR10"]),
    nrow(tair_genes_f),
    0
  ),
  Genes_filtered_Kas1  = c(
    nrow(core_genes_f[Genome == "Yuwei"]),
    0,
    nrow(yuwei_genes_f)
  ),
  stringsAsFactors = FALSE
)

cat("\n[SUMMARY TABLE]\n")
print(summary_df)

fwrite(summary_df, OUT_SUMMARY)
cat("\n[OK] Summary written to:\n  ", OUT_SUMMARY, "\n")