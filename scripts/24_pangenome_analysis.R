suppressPackageStartupMessages({
  library(GENESPACE)
})

# Load GENESPACE output RDS
pan <- readRDS("/data/users/yliu2/Organization_and_annotation/genespace/pangenome_matrix.rds")

# Check structure
cat("[INFO] Loaded pangenome object:\n")
print(names(pan))

# Extract presence/absence matrix
pam <- pan$pangeneMat
cat("[INFO] Presence/absence matrix dimensions: ", dim(pam), "\n")

# Define core and unique orthogroups (TAIR10 vs Yuwei)
core_ogs      <- rownames(pam)[pam[,"TAIR10"]==1 & pam[,"Yuwei"]==1]
tair_only_ogs <- rownames(pam)[pam[,"TAIR10"]==1 & pam[,"Yuwei"]==0]
yuwei_only_ogs<- rownames(pam)[pam[,"TAIR10"]==0 & pam[,"Yuwei"]==1]

cat("Core orthogroups:", length(core_ogs), "\n")
cat("TAIR10-unique orthogroups:", length(tair_only_ogs), "\n")
cat("Yuwei-unique orthogroups:", length(yuwei_only_ogs), "\n")

# Map orthogroups back to genes
og2genes <- pan$og2Gene
core_genes       <- og2genes[og2genes$og %in% core_ogs, ]
tair_only_genes  <- og2genes[og2genes$og %in% tair_only_ogs & og2genes$genome=="TAIR10", ]
yuwei_only_genes <- og2genes[og2genes$og %in% yuwei_only_ogs & og2genes$genome=="Yuwei", ]

# Output summary table
summary_df <- data.frame(
  Category = c("Core", "TAIR10_unique", "Yuwei_unique"),
  Orthogroups = c(length(core_ogs), length(tair_only_ogs), length(yuwei_only_ogs)),
  Genes = c(nrow(core_genes), nrow(tair_only_genes), nrow(yuwei_only_genes))
)
print(summary_df)

# Save results
write.csv(summary_df,
  file="/data/users/yliu2/Organization_and_annotation/genespace/pangenome_summary.csv",
  row.names=FALSE
)
cat("[OK] Summary table saved: pangenome_summary.csv\n")
