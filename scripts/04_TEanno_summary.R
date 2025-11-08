#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

option_list <- list(
  make_option("--sum", type="character", help="EDTA .TEanno.sum file"),
  make_option("--outdir", type="character", default=NULL, help="Output dir [default: next to input]/summary"),
  make_option("--prefix", type="character", default="TEanno", help="Output prefix [default: %default]"),
  make_option("--fai", type="character", default=NA, help="Optional .fai for genome size"),
  make_option("--topN", type="integer", default=12, help="Top-N superfamilies to plot [default: %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))
stopifnot(!is.null(opt$sum))

if (is.null(opt$outdir)) opt$outdir <- file.path(dirname(opt$sum), "summary")
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

sum_file <- opt$sum
out_prefix <- file.path(opt$outdir, opt$prefix)
message("Reading: ", sum_file)

raw_lines <- readLines(sum_file, warn = FALSE)

# genome size (optional)
genome_size <- NA_real_
gs_line <- raw_lines[grepl("Genome size", raw_lines, ignore.case = TRUE)]
if (length(gs_line) > 0) {
  gs_num <- str_extract(gs_line[1], "\\d+(?:,\\d+)*")
  if (!is.na(gs_num)) genome_size <- as.numeric(gsub(",", "", gs_num))
}
if (is.na(genome_size) && !is.na(opt$fai) && file.exists(opt$fai)) {
  fai <- read_tsv(opt$fai, col_names = FALSE, show_col_types = FALSE)
  genome_size <- sum(fai$X2, na.rm=TRUE)
}
if (is.na(genome_size)) {
  message("Warning: genome size not found; will trust percent columns from .sum")
}

# parse per-line summary
# expected lines like: "LTR/Copia   123   4567890   12.34%"
pattern <- "^\\s*([[:alnum:]_./+-]+)\\s+(\\d+)\\s+(\\d+)\\s+([0-9.]+)%"
keep <- grepl(pattern, raw_lines)
tbl_lines <- raw_lines[keep]
if (length(tbl_lines) == 0) stop("No superfamily summary lines detected.")

df <- stringr::str_match(tbl_lines, pattern) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  transmute(
    superfamily_raw = V2,
    count = as.numeric(V3),
    bp = as.numeric(V4),
    percent = as.numeric(V5)
  )

# split into class and subtype; normalize labels
df <- df %>%
  mutate(
    class = toupper(ifelse(grepl("/", superfamily_raw), sub("/.*$", "", superfamily_raw), superfamily_raw)),
    subtype = ifelse(grepl("/", superfamily_raw), sub("^.*?/", "", superfamily_raw), superfamily_raw),
    superfamily = ifelse(grepl("/", superfamily_raw),
                         paste0(toupper(sub("/.*$", "", superfamily_raw)), "/", subtype),
                         toupper(superfamily_raw))
  )

# --- CLASS-LEVEL summary (clean, small number of categories) ---
class_sum <- df %>%
  group_by(class) %>%
  summarise(count = sum(count, na.rm=TRUE),
            bp = sum(bp, na.rm=TRUE),
            percent = sum(percent, na.rm=TRUE),
            .groups = "drop") %>%
  arrange(desc(percent))

# write table
class_csv <- paste0(out_prefix, "_class_summary.csv")
write_csv(class_sum, class_csv)
message("Saved: ", class_csv)

# plot class-level composition
p_class <- ggplot(class_sum, aes(x = reorder(class, percent), y = percent)) +
  geom_col(fill = "#4E79A7", width = 0.7) +
  coord_flip() +
  labs(x = "Class", y = "Genome %", title = "TE composition by Class") +
  theme_minimal(base_size = 12)
ggsave(paste0(out_prefix, "_class_percent_bar.pdf"), p_class, width = 6, height = 4)
ggsave(paste0(out_prefix, "_class_percent_bar.png"), p_class, width = 6, height = 4, dpi = 300)

# --- SUPERFAMILY-LEVEL summary (aggregate, top N + Others) ---
sf_sum <- df %>%
  group_by(superfamily) %>%
  summarise(count = sum(count, na.rm=TRUE),
            bp = sum(bp, na.rm=TRUE),
            percent = sum(percent, na.rm=TRUE),
            .groups = "drop") %>%
  arrange(desc(percent))

sf_csv <- paste0(out_prefix, "_superfamily_summary.csv")
write_csv(sf_sum, sf_csv)
message("Saved: ", sf_csv)

topN <- min(opt$topN, nrow(sf_sum))
sf_top <- sf_sum %>% slice_head(n = topN)
others_pct <- sum(sf_sum$percent, na.rm=TRUE) - sum(sf_top$percent, na.rm=TRUE)
others_row <- if (others_pct > 0) tibble(superfamily = "Others", count = NA_real_, bp = NA_real_, percent = others_pct) else NULL
sf_plot <- bind_rows(sf_top, others_row) %>%
  mutate(superfamily = factor(superfamily, levels = rev(superfamily)))

p_sf <- ggplot(sf_plot, aes(x = superfamily, y = percent)) +
  geom_col(fill = "#59A14F", width = 0.7) +
  coord_flip() +
  labs(x = "Superfamily (Top-N)", y = "Genome %", title = sprintf("Top-%d TE superfamilies", topN)) +
  theme_minimal(base_size = 12)
ggsave(paste0(out_prefix, "_superfamily_top", topN, "_percent_bar.pdf"), p_sf, width = 7, height = 6)
ggsave(paste0(out_prefix, "_superfamily_top", topN, "_percent_bar.png"), p_sf, width = 7, height = 6, dpi = 300)

# optional totals
if (!is.na(genome_size)) {
  total_bp <- sum(df$bp, na.rm=TRUE)
  total_pct <- total_bp / genome_size * 100
  totals_file <- paste0(out_prefix, "_totals.txt")
  cat(sprintf("Genome_size_bp\t%.0f\nTotal_TE_bp\t%.0f\nTotal_TE_%%\t%.3f\n",
              genome_size, total_bp, total_pct),
      file = totals_file)
  message("Saved: ", totals_file)
}

message("\nCLASS summary (top):")
print(class_sum)

message("\nSUPERFAMILY summary (top rows):")
print(sf_sum %>% slice_head(n = 15))