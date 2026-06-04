#!/usr/bin/env bash
# Jalankan pipeline HISAT2 -> StringTie -> prepDE pada data demo RNA-Seq untuk
# menghasilkan docs/data/rnaseq-demo/gene_count_matrix.csv. Pipeline ini memakai
# tools yang sama persis seperti yang ditampilkan pada tutorial. Dijalankan
# sekali; hasil count matrix di-commit sehingga build situs tidak perlu tools ini.
# Prasyarat: hisat2, stringtie, samtools, python3 tersedia di PATH.
# Jalankan dari root repo:  bash code/run-rnaseq-pipeline.sh
set -euo pipefail

DD=docs/data/rnaseq-demo
REF=$DD/ref
WORK=$DD/_pipe
rm -rf "$WORK"; mkdir -p "$WORK/bam" "$WORK/ballgown"

# prepDE.py3 (skrip resmi StringTie) diunduh bila belum ada
if [ ! -f "$WORK/prepDE.py3" ]; then
  curl -fsSL https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3 -o "$WORK/prepDE.py3"
fi

# 1. Index genom referensi
hisat2-build -q "$REF/demo_genome.fa" "$WORK/idx"

# daftar sampel dari metadata
samples=$(tail -n +2 "$DD/metadata.csv" | cut -d, -f1 | tr -d '"')

: > "$WORK/sample_list.txt"
for s in $samples; do
  # 2. Alignment (HISAT2) + sorting (samtools)
  hisat2 -p 4 --no-spliced-alignment -x "$WORK/idx" \
    -1 "$DD/raw_reads/${s}_R1.fastq.gz" -2 "$DD/raw_reads/${s}_R2.fastq.gz" \
    -S "$WORK/bam/${s}.sam" 2> "$WORK/bam/${s}.hisat2.log"
  samtools sort -@ 4 -o "$WORK/bam/${s}.bam" "$WORK/bam/${s}.sam"
  rm -f "$WORK/bam/${s}.sam"
  samtools index "$WORK/bam/${s}.bam"

  # 3. Kuantifikasi (StringTie -e -B mengikuti anotasi referensi)
  mkdir -p "$WORK/ballgown/${s}"
  stringtie -e -B -p 4 -G "$REF/demo_annotation.gtf" \
    -o "$WORK/ballgown/${s}/${s}.gtf" "$WORK/bam/${s}.bam"
  printf '%s\t%s\n' "$s" "$WORK/ballgown/${s}/${s}.gtf" >> "$WORK/sample_list.txt"
done

# 4. Matriks hitungan gen (prepDE.py)
python3 "$WORK/prepDE.py3" -i "$WORK/sample_list.txt" \
  -g "$DD/gene_count_matrix.csv" -t "$WORK/transcript_count_matrix.csv"

echo "DONE -> $DD/gene_count_matrix.csv"
