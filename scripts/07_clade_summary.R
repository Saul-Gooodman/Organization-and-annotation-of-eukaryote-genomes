#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(readr)
})

option_list <- list(
  make_option("--workdir", type="character",
              default="/data/users/yliu2/Organization_and_annotation",
              help="Base workdir"),
  make_option("--edta_dir", type="character",
              default=NULL,
              help="EDTA result dir [default: <workdir>/results/EDTA_annotation]"),
  make_option("--cls_dir", type="character",
              default=NULL,
              help="Directory containing *.rexdb-plant.cls.tsv [default: search under <edta_dir>]"),
  make_option("--gff", type="character",
              default=NULL,
              help="EDTA TEanno.gff3 [default: <edta_dir>/assembly.p_ctg.fa.mod.EDTA.TEanno.gff3]"),
  make_option("--outdir", type="character",
              default=NULL,
              help="Output dir [default: <edta_dir>/te_clades]"),
  make_option("--prefix", type="character",
              default="assembly_p_ctg",
              help="Output file prefix")
)
opt <- parse_args(OptionParser(option_list = option_list))

workdir <- opt$workdir
edta_dir <- if (is.null(opt$edta_dir)) file.path(workdir, "results/EDTA_annotation") else opt$edta_dir
gff_file <- if (is.null(opt$gff)) file.path(edta_dir, "assembly.p_ctg.fa.mod.EDTA.TEanno.gff3") else opt$gff
outdir   <- if (is.null(opt$outdir)) file.path(edta_dir, "te_clades") else opt$outdir
cls_dir  <- if (is.null(opt$cls_dir)) edta_dir else opt$cls_dir

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# find cls.tsv (recursively, robust)
cls_files <- list.files(cls_dir, pattern="\\.rexdb-plant\\.cls\\.tsv$", full.names = TRUE, recursive = TRUE)
if (length(cls_files) == 0) {
  message("No *.rexdb-plant.cls.tsv found under: ", cls_dir)
  quit(save="no", status=0)
}
message("Using cls files:")
message(paste0("  - ", paste(basename(cls_files), collapse = "\n  - ")))

# reader with column name harmonization
read_cls <- function(f){
  dt <- suppressWarnings(read_tsv(f, show_col_types = FALSE, progress = FALSE))
  nm <- tolower(names(dt))
  names(dt) <- nm

  # harmonize "query" column
  if (!"query" %in% nm) {
    cand <- c("qseqid","id","seqid","name","qname","entry","seq_name")
    has  <- intersect(cand, nm)
    if (length(has)) dt$query <- dt[[has[1]]] else dt$query <- NA_character_
  }
  # make sure we have superfamily and clade
  if (!"superfamily" %in% nm) dt$superfamily <- NA_character_
  if (!"clade"       %in% nm) dt$clade <- NA_character_

  dt$source_file <- basename(f)
  dt
}

cls_all <- bind_rows(lapply(cls_files, read_cls)) %>%
  mutate(class = case_when(
           str_detect(superfamily, regex("copia", ignore_case = TRUE)) ~ "Copia",
           str_detect(superfamily, regex("gypsy", ignore_case = TRUE)) ~ "Gypsy",
           str_detect(source_file, regex("copia", ignore_case = TRUE)) ~ "Copia",
           str_detect(source_file, regex("gypsy", ignore_case = TRUE)) ~ "Gypsy",
           TRUE ~ "Unknown"
         ),
         clade = if_else(is.na(clade) | clade=="", "UNCLASSIFIED", clade))

# library entries per clade
lib_summary <- cls_all %>%
  count(class, clade, name = "library_entries") %>%
  arrange(class, desc(library_entries))

write_csv(lib_summary, file.path(outdir, paste0(opt$prefix, "_clade_library_entries.csv")))

# quick barplots (top N per class)
plot_top <- function(df, clazz, topn=12){
  sub <- df %>% filter(class==clazz) %>% arrange(desc(library_entries)) %>% slice_head(n=topn)
  if (nrow(sub) == 0) return(NULL)
  ggplot(sub, aes(x=reorder(clade, library_entries), y=library_entries)) +
    geom_col() +
    coord_flip() +
    labs(title=paste(clazz, "clades (library entries, Top", topn, ")"),
         x="Clade", y="Count") +
    theme_bw()
}
p1 <- plot_top(lib_summary, "Copia", 12)
p2 <- plot_top(lib_summary, "Gypsy", 12)

pdf(file.path(outdir, paste0(opt$prefix, "_clade_library_entries_top12.pdf")), width=7, height=8)
if (!is.null(p1)) print(p1)
if (!is.null(p2)) print(p2)
dev.off()

# genome copy counts by joining GFF if available
if (file.exists(gff_file)) {
  # fread with fill to tolerate header/comment lines
  gff <- tryCatch(
    fread(gff_file, sep="\t", header=FALSE, quote="", fill=TRUE,
          col.names=c("seqid","source","type","start","end","score","strand","phase","attr")),
    error = function(e) {
      # fallback: read lines, drop comment-lines, then fread text
      x <- readLines(gff_file, warn = FALSE)
      x <- x[!startsWith(x, "#")]
      fread(text = x, sep = "\t", header = FALSE, quote = "", fill = TRUE,
            col.names=c("seqid","source","type","start","end","score","strand","phase","attr"))
    }
  )

  gff <- gff %>%
    filter(type %in% c("transposable_element","repeat_region","dispersed_repeat","mobile_genetic_element") |
             grepl("TE", type, ignore.case=TRUE) |
             grepl("repeat", type, ignore.case=TRUE))

  # extract a library ID-like field from attributes
  get_name_like <- function(s){
    nm <- str_match(s, "Name=([^;]+)")[,2]; if (!is.na(nm)) return(nm)
    tg <- str_match(s, "Target=([^;\\s]+)")[,2]; if (!is.na(tg)) return(tg)
    id <- str_match(s, "ID=([^;]+)")[,2];     if (!is.na(id)) return(id)
    cl <- str_match(s, "Classification=([^;]+)")[,2]; if (!is.na(cl)) return(cl)
    return(NA_character_)
  }
  gff$lib_id <- vapply(gff$attr, get_name_like, FUN.VALUE=character(1))
  gff$lib_id <- gsub("\\s.*$", "", gff$lib_id)
  gff$lib_id <- gsub(",.*$", "", gff$lib_id)

  # clean query to improve join
  cls_all$query_clean <- gsub("\\s.*$", "", cls_all$query)
  cls_all$query_clean <- gsub(",.*$", "", cls_all$query_clean)

  gff_join <- gff %>%
    left_join(cls_all %>% select(query_clean, class, clade) %>% distinct(),
              by = c("lib_id" = "query_clean"))

  copies_summary <- gff_join %>%
    mutate(class = if_else(is.na(class), "Unknown", class),
           clade = if_else(is.na(clade), "UNCLASSIFIED_OR_NOHIT", clade)) %>%
    count(class, clade, name = "genome_copies") %>%
    arrange(class, desc(genome_copies))

  write_csv(copies_summary, file.path(outdir, paste0(opt$prefix, "_clade_genome_copies.csv")))

  plot_top_copies <- function(df, clazz, topn=12){
    sub <- df %>% filter(class==clazz) %>% arrange(desc(genome_copies)) %>% slice_head(n=topn)
    if (nrow(sub) == 0) return(NULL)
    ggplot(sub, aes(x=reorder(clade, genome_copies), y=genome_copies)) +
      geom_col() +
      coord_flip() +
      labs(title=paste(clazz, "clades (genome copies, Top", topn,")"),
           x="Clade", y="Copies") +
      theme_bw()
  }
  q1 <- plot_top_copies(copies_summary, "Copia", 12)
  q2 <- plot_top_copies(copies_summary, "Gypsy", 12)
  pdf(file.path(outdir, paste0(opt$prefix, "_clade_genome_copies_top12.pdf")), width=7, height=8)
  if (!is.null(q1)) print(q1)
  if (!is.null(q2)) print(q2)
  dev.off()
} else {
  message("GFF not found: ", gff_file, " -> skip genome copies summary.")
}

message("Done. Outputs in: ", outdir)