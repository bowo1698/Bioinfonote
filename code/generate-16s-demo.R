#!/usr/bin/env Rscript
# Synthetic 16S rRNA V3-V4 demo dataset for the Bioinfonote metagenomics tutorial.
# Produces paired-end FASTQ for 20 samples (10 sehat, 10 AHPND), a small DADA2
# taxonomy training set, and metadata. Deterministic (seeded).
# Run from repo root:  Rscript code/generate-16s-dummy.R

set.seed(20260604)

out_dir   <- "docs/data/16s-demo"
reads_dir <- file.path(out_dir, "raw_reads")
ref_dir   <- file.path(out_dir, "ref")
dir.create(reads_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ref_dir,   recursive = TRUE, showWarnings = FALSE)

## 1. Taxa table (genus + 6-level lineage) ---------------------------------
# Semua genus adalah bakteri laut/payau yang lazim pada lingkungan tambak dan
# saluran pencernaan udang vanamei (tidak ada genus khas air tawar).
tax <- read.csv(text = "
genus,kingdom,phylum,class,order,family
Vibrio,Bacteria,Proteobacteria,Gammaproteobacteria,Vibrionales,Vibrionaceae
Photobacterium,Bacteria,Proteobacteria,Gammaproteobacteria,Vibrionales,Vibrionaceae
Aliivibrio,Bacteria,Proteobacteria,Gammaproteobacteria,Vibrionales,Vibrionaceae
Shewanella,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Shewanellaceae
Pseudoalteromonas,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Pseudoalteromonadaceae
Alteromonas,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae
Colwellia,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Colwelliaceae
Marinobacter,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Marinobacteraceae
Pseudomonas,Bacteria,Proteobacteria,Gammaproteobacteria,Pseudomonadales,Pseudomonadaceae
Psychrobacter,Bacteria,Proteobacteria,Gammaproteobacteria,Pseudomonadales,Moraxellaceae
Halomonas,Bacteria,Proteobacteria,Gammaproteobacteria,Oceanospirillales,Halomonadaceae
Marinomonas,Bacteria,Proteobacteria,Gammaproteobacteria,Oceanospirillales,Oceanospirillaceae
Tenacibaculum,Bacteria,Bacteroidota,Flavobacteriia,Flavobacteriales,Flavobacteriaceae
Polaribacter,Bacteria,Bacteroidota,Flavobacteriia,Flavobacteriales,Flavobacteriaceae
Muricauda,Bacteria,Bacteroidota,Flavobacteriia,Flavobacteriales,Flavobacteriaceae
Ruegeria,Bacteria,Proteobacteria,Alphaproteobacteria,Rhodobacterales,Rhodobacteraceae
Phaeobacter,Bacteria,Proteobacteria,Alphaproteobacteria,Rhodobacterales,Rhodobacteraceae
Sulfitobacter,Bacteria,Proteobacteria,Alphaproteobacteria,Rhodobacterales,Rhodobacteraceae
Roseobacter,Bacteria,Proteobacteria,Alphaproteobacteria,Rhodobacterales,Rhodobacteraceae
Bacillus,Bacteria,Firmicutes,Bacilli,Bacillales,Bacillaceae
", stringsAsFactors = FALSE, strip.white = TRUE)
tax <- tax[tax$genus != "", ]
n_taxa <- nrow(tax)

## 2. Reference amplicon sequences (synthetic V3-V4, 400 bp) ----------------
La    <- 400
bases <- c("A","C","G","T")
backbone <- sample(bases, La, replace = TRUE)   # conserved across taxa
n_var  <- 150
var_pos <- sort(sample(seq_len(La), n_var))     # variable positions
make_seq <- function(i) {
  s <- backbone
  s[var_pos] <- sample(bases, n_var, replace = TRUE)  # taxon-specific variants
  paste(s, collapse = "")
}
ref_seqs <- character(n_taxa)
for (i in seq_len(n_taxa)) { set.seed(1000 + i); ref_seqs[i] <- make_seq(i) }
names(ref_seqs) <- tax$genus
set.seed(20260604)

## 3. DADA2 training fasta --------------------------------------------------
hdr <- sprintf(">%s;%s;%s;%s;%s;%s;",
               tax$kingdom, tax$phylum, tax$class, tax$order, tax$family, tax$genus)
ref_fa <- file.path(ref_dir, "train_set_demo.fa.gz")
con <- gzfile(ref_fa, "w")
writeLines(as.vector(rbind(hdr, unname(ref_seqs))), con)
close(con)

## 4. Per-group composition (Dirichlet) -------------------------------------
rdirichlet <- function(alpha) { g <- rgamma(length(alpha), alpha, 1); g / sum(g) }
alpha_sehat <- rep(0.7, n_taxa); names(alpha_sehat) <- tax$genus
alpha_sehat[c("Bacillus","Phaeobacter","Ruegeria","Shewanella","Pseudoalteromonas")] <- 2.5
alpha_ahpnd <- rep(0.35, n_taxa); names(alpha_ahpnd) <- tax$genus
alpha_ahpnd[c("Vibrio","Photobacterium")] <- c(9, 5)
alpha_ahpnd[c("Aliivibrio","Tenacibaculum")] <- 1.5

## 5. Primers + helpers -----------------------------------------------------
primerF <- "CCTACGGGAGGCAGCAG"        # 341F, 17 nt
primerR <- "GACTACCAGGGTATCTAATCC"    # 805R, 21 nt
Lb <- 215                              # biological read length (post primer trim)
revcomp <- function(x) chartr("ACGT","TGCA", paste(rev(strsplit(x,"")[[1]]), collapse=""))
mutate_seq <- function(s, err = 0.005) {
  v <- strsplit(s,"")[[1]]
  idx <- which(runif(length(v)) < err)
  for (j in idx) v[j] <- sample(setdiff(bases, v[j]), 1)
  paste(v, collapse = "")
}
qual_str <- function(n) {                     # ~12% reads berkualitas rendah (gagal filter maxEE)
  # dua profil tetap (entropi rendah, kompresi gzip baik); ~12% read "buruk" gagal filter
  if (runif(1) < 0.12) intToUtf8(33 + round(seq(30, 8,  length.out = n)))   # bad read
  else                 intToUtf8(33 + round(seq(38, 28, length.out = n)))   # good read
}

## 6. Samples + reads -------------------------------------------------------
samples <- c(sprintf("SH%02d", 1:10), sprintf("AH%02d", 1:10))
groups  <- c(rep("sehat", 10), rep("AHPND", 10))
write.csv(data.frame(sample_id = samples, group = groups),
          file.path(out_dir, "metadata.csv"), row.names = FALSE)

for (k in seq_along(samples)) {
  smp <- samples[k]; grp <- groups[k]
  alpha  <- if (grp == "sehat") alpha_sehat else alpha_ahpnd
  p      <- rdirichlet(alpha)
  N      <- rpois(1, 2000)
  counts <- as.vector(rmultinom(1, N, p))
  con1 <- gzfile(file.path(reads_dir, paste0(smp, "_1.fastq.gz")), "w")
  con2 <- gzfile(file.path(reads_dir, paste0(smp, "_2.fastq.gz")), "w")
  buf1 <- character(0); buf2 <- character(0); rid <- 0
  for (t in seq_len(n_taxa)) {
    if (counts[t] == 0) next
    A       <- ref_seqs[t]
    fwd_bio <- substr(A, 1, Lb)
    rev_bio <- revcomp(substr(A, La - Lb + 1, La))
    for (r in seq_len(counts[t])) {
      rid <- rid + 1
      id  <- sprintf("@%s.%d", smp, rid)
      r1  <- mutate_seq(paste0(primerF, fwd_bio))
      r2  <- mutate_seq(paste0(primerR, rev_bio))
      buf1 <- c(buf1, paste(id, "1:N:0:1"), r1, "+", qual_str(nchar(r1)))
      buf2 <- c(buf2, paste(id, "2:N:0:1"), r2, "+", qual_str(nchar(r2)))
    }
  }
  writeLines(buf1, con1); writeLines(buf2, con2)
  close(con1); close(con2)
  cat("wrote", smp, "(", sum(counts), "reads )\n")
}

## 7. Bundel reads untuk tombol unduh ---------------------------------------
local({
  wd <- setwd(out_dir); on.exit(setwd(wd))
  utils::tar("raw_reads.tar.gz", "raw_reads", compression = "gzip")
})

cat("DONE. reads:", reads_dir, " ref:", ref_fa,
    " bundle:", file.path(out_dir, "raw_reads.tar.gz"), "\n")
