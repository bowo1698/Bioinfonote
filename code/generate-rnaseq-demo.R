#!/usr/bin/env Rscript
# Demo RNA-seq dataset for the Bioinfonote RNA-Seq tutorial.
# Produces a small reference genome (FASTA) + annotation (GTF), paired-end
# FASTQ for 10 samples (5 kontrol, 5 infeksi) with a built-in differential-
# expression signal, and sample metadata. The HISAT2 -> StringTie -> prepDE
# pipeline run on these inputs yields gene_count_matrix.csv, which the live
# DESeq2 analysis on the tutorial page consumes. Deterministic.
# Run from repo root:  Rscript code/generate-rnaseq-demo.R

set.seed(20260604)

out_dir   <- "docs/data/rnaseq-demo"
reads_dir <- file.path(out_dir, "raw_reads")
ref_dir   <- file.path(out_dir, "ref")
dir.create(reads_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ref_dir,   recursive = TRUE, showWarnings = FALSE)

bases <- c("A", "C", "G", "T")
revcomp <- function(x) chartr("ACGT", "TGCA",
                              paste(rev(strsplit(x, "")[[1]]), collapse = ""))

## 1. Gen sintetik = genom referensi (tiap gen satu contig) -----------------
n_genes  <- 150
gene_ids <- sprintf("gene%03d", seq_len(n_genes))
gene_len <- sample(600:1400, n_genes, replace = TRUE)
gene_seq <- vapply(gene_len,
                   function(L) paste(sample(bases, L, replace = TRUE), collapse = ""),
                   character(1))
names(gene_seq) <- gene_ids

# genom referensi FASTA (sekuens dibungkus 70 char/baris)
fa <- file(file.path(ref_dir, "demo_genome.fa"), "w")
for (i in seq_len(n_genes)) {
  writeLines(paste0(">", gene_ids[i]), fa)
  s <- gene_seq[i]
  st <- seq(1, nchar(s), by = 70)
  writeLines(substring(s, st, pmin(st + 69L, nchar(s))), fa)
}
close(fa)

# anotasi GTF (tiap gen = satu transkrip satu exon membentang penuh contig)
gtf <- file(file.path(ref_dir, "demo_annotation.gtf"), "w")
for (i in seq_len(n_genes)) {
  attr <- sprintf('gene_id "%s"; transcript_id "%s.1";', gene_ids[i], gene_ids[i])
  writeLines(sprintf("%s\tdemo\ttranscript\t1\t%d\t.\t+\t.\t%s",
                     gene_ids[i], gene_len[i], attr), gtf)
  writeLines(sprintf('%s\tdemo\texon\t1\t%d\t.\t+\t.\t%s exon_number "1";',
                     gene_ids[i], gene_len[i], attr), gtf)
}
close(gtf)

## 2. Desain & ekspresi dasar ----------------------------------------------
samples <- c(sprintf("kontrol%d", 1:5), sprintf("infeksi%d", 1:5))
groups  <- c(rep("kontrol", 5), rep("infeksi", 5))
write.csv(data.frame(sample = samples, group = groups),
          file.path(out_dir, "metadata.csv"), row.names = FALSE, quote = FALSE)

base_mu <- round(exp(rnorm(n_genes, log(40), 0.7)))
base_mu[base_mu < 6] <- 6
n_de   <- 40
de_idx <- sample(n_genes, n_de)
lfc    <- numeric(n_genes)
lfc[de_idx] <- sample(c(-1, 1), n_de, replace = TRUE) * runif(n_de, 1.5, 3.2)
size_factor <- setNames(runif(length(samples), 0.85, 1.15), samples)

## 3. Helpers reads ---------------------------------------------------------
rl <- 100
frag <- 250
mutate_seq <- function(s, err = 0.003) {
  v <- strsplit(s, "")[[1]]
  idx <- which(runif(length(v)) < err)
  for (j in idx) v[j] <- sample(setdiff(bases, v[j]), 1)
  paste(v, collapse = "")
}
qual <- strrep("I", rl)

## 4. Tulis reads -----------------------------------------------------------
for (k in seq_along(samples)) {
  smp <- samples[k]
  grp <- groups[k]
  con1 <- gzfile(file.path(reads_dir, paste0(smp, "_R1.fastq.gz")), "w")
  con2 <- gzfile(file.path(reads_dir, paste0(smp, "_R2.fastq.gz")), "w")
  buf1 <- character(0)
  buf2 <- character(0)
  rid <- 0
  for (g in seq_len(n_genes)) {
    mu <- base_mu[g] * size_factor[smp] * if (grp == "infeksi") 2^lfc[g] else 1
    n_pairs <- rnbinom(1, mu = mu, size = 8)
    if (n_pairs == 0) next
    L <- gene_len[g]
    fr <- min(frag, L)
    for (r in seq_len(n_pairs)) {
      rid <- rid + 1
      start <- if (L > fr) sample(seq_len(L - fr + 1), 1) else 1
      r1 <- mutate_seq(substr(gene_seq[g], start, start + rl - 1))
      r2 <- mutate_seq(revcomp(substr(gene_seq[g], start + fr - rl, start + fr - 1)))
      id <- sprintf("@%s.%d", smp, rid)
      buf1 <- c(buf1, paste(id, "1:N:0:1"), r1, "+", qual)
      buf2 <- c(buf2, paste(id, "2:N:0:1"), r2, "+", qual)
    }
  }
  writeLines(buf1, con1)
  writeLines(buf2, con2)
  close(con1)
  close(con2)
  cat("wrote", smp, "(", rid, "read pairs )\n")
}

## 5. Bundel reads untuk tombol unduh ---------------------------------------
local({
  wd <- setwd(out_dir)
  on.exit(setwd(wd))
  utils::tar("raw_reads.tar.gz", "raw_reads", compression = "gzip")
})

cat("DONE. genome+gtf:", ref_dir, " reads:", reads_dir, "\n")
