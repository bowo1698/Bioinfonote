#!/bin/bash

# ===================================
# FASTQ Cleaning Pipeline for Paired-End Data
# ===================================
# Script untuk memproses multiple file FASTQ paired-end
# - Mengatasi kualitas basis menurun
# - Menghilangkan adapter
# - Standardisasi panjang read
# - Menghasilkan laporan kualitas 
# - Jalankan dengan perintah: bash process_paired_fastq.sh -i direktori_fastq_input -o direktori_output
# - Lihat dokumentasi: bash process_paired_fastq.sh -h
# ===================================

echo "==================================="
echo "FASTQ Cleaning Pipeline for Paired-End Data"
echo "==================================="
echo "Author: Agus Wibowo"
echo "Version: 1.0"
echo "Last updated: 2025-05-20"
echo "==================================="
echo "Lihat dokumentasi: bash process_paired_fastq.sh -h"
echo "==================================="

# Default values
INPUT_DIR=${INPUT_DIR%/}
OUTPUT_DIR=${OUTPUT_DIR%/}
THREADS=4
QUALITY=20
MIN_LENGTH=36
ADAPTER="AGATCGGAAGAGC"

# Fungsi untuk menampilkan bar progres
show_progress() {
    local percent=$1
    local width=50
    local num_filled=$(( width * percent / 100 ))
    local num_empty=$(( width - num_filled ))
    
    # Buat bar progres dengan karakter blok
    local bar=""
    for ((i=0; i<num_filled; i++)); do
        bar="${bar}█"
    done
    for ((i=0; i<num_empty; i++)); do
        bar="${bar}░"
    done
    
    # Tampilkan bar progres dengan persentase
    printf "\r[%s] %3d%%" "$bar" "$percent"
}

# Fungsi untuk menampilkan bantuan penggunaan
show_usage() {
    echo "Penggunaan: $0 [opsi]"
    echo ""
    echo "Opsi:"
    echo "  -i, --input DIR      Direktori input yang berisi file FASTQ (default: ./raw_data)"
    echo "  -o, --output DIR     Direktori output untuk file hasil (default: ./cleaned_data)"
    echo "  -t, --threads NUM    Jumlah thread yang digunakan (default: 4)"
    echo "  -q, --quality NUM    Nilai kualitas minimum (default: 20)"
    echo "  -l, --length NUM     Panjang minimum read (default: 36)"
    echo "  -a, --adapter SEQ    Sekuens adapter (default: AGATCGGAAGAGC)"
    echo "  -h, --help           Tampilkan bantuan ini"
    echo ""
    echo "Contoh:"
    echo "  $0 -i ./fastq_raw -o ./fastq_clean -t 8 -q 30"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -q|--quality)
            QUALITY="$2"
            shift 2
            ;;
        -l|--length)
            MIN_LENGTH="$2"
            shift 2
            ;;
        -a|--adapter)
            ADAPTER="$2"
            shift 2
            ;;
        -h|--help)
            show_usage
            ;;
        *)
            echo "Opsi tidak dikenal: $1"
            show_usage
            ;;
    esac
done

# Validasi direktori input
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Direktori input '$INPUT_DIR' tidak ditemukan."
    exit 1
fi

# Buat direktori output jika belum ada
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/reports"
mkdir -p "${OUTPUT_DIR}/logs"

# Nama file log
DATE_STAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${OUTPUT_DIR}/logs/processing_log_${DATE_STAMP}.txt"

# Header log
echo "===================================" | tee $LOG_FILE
echo "FASTQ Processing Pipeline" | tee -a $LOG_FILE
echo "Started at: $(date)" | tee -a $LOG_FILE
echo "===================================" | tee -a $LOG_FILE
echo "Input directory: $INPUT_DIR" | tee -a $LOG_FILE
echo "Output directory: $OUTPUT_DIR" | tee -a $LOG_FILE
echo "Threads: $THREADS" | tee -a $LOG_FILE
echo "Quality threshold: $QUALITY" | tee -a $LOG_FILE
echo "Minimum read length: $MIN_LENGTH" | tee -a $LOG_FILE
echo "Adapter sequence: $ADAPTER" | tee -a $LOG_FILE
echo "===================================" | tee -a $LOG_FILE

# Hitung jumlah file R1 untuk laporan
TOTAL_PAIRS=$(ls $INPUT_DIR/*_R1*.fastq.gz 2>/dev/null | wc -l)
if [ $TOTAL_PAIRS -eq 0 ]; then
  echo "PERINGATAN: Tidak ada file *_R1*.fastq.gz yang ditemukan di $INPUT_DIR" | tee -a $LOG_FILE
  echo "Periksa apakah direktori input benar dan berisi file FASTQ" | tee -a $LOG_FILE
  exit 1
fi

echo "Total pasangan file yang akan diproses: $TOTAL_PAIRS" | tee -a $LOG_FILE
echo "===================================" | tee -a $LOG_FILE

# Variabel untuk menghitung progres
CURRENT_PAIR=0
TOTAL_STEPS=$(( TOTAL_PAIRS * 2 )) # Jumlah total langkah (fastp + fastqc untuk setiap pasangan)
CURRENT_STEP=0

# Header untuk progress bar
echo "Progress:"

# Cari semua file R1
for R1_FILE in $INPUT_DIR/*_R1*.fastq.gz
do
  # Tambahkan counter
  ((CURRENT_PAIR++))
  
  # Nama file tanpa path
  R1_FILENAME=$(basename "$R1_FILE")
  
  # Konstruksi nama R2 yang sesuai
  R2_FILE="${R1_FILE/_R1/_R2}"
  R2_FILENAME=$(basename "$R2_FILE")
  
  # Verifikasi R2 ada
  if [ ! -f "$R2_FILE" ]; then
    echo -e "\n[ERROR] Pasangan file tidak ditemukan untuk $R1_FILENAME" | tee -a $LOG_FILE
    echo "        Harap periksa apakah $R2_FILENAME ada" | tee -a $LOG_FILE
    echo "        Melanjutkan ke file berikutnya..." | tee -a $LOG_FILE
    echo "" | tee -a $LOG_FILE
    # Update progress bar saat melewati file
    CURRENT_STEP=$(( CURRENT_STEP + 2 ))
    PERCENT=$(( CURRENT_STEP * 100 / TOTAL_STEPS ))
    show_progress $PERCENT
    continue
  fi
  
  # Buat nama untuk file output
  SAMPLE_NAME=$(basename "$R1_FILE" | sed 's/_R1.*//g')
  OUT_R1="${OUTPUT_DIR}/${SAMPLE_NAME}_R1_cleaned.fastq.gz"
  OUT_R2="${OUTPUT_DIR}/${SAMPLE_NAME}_R2_cleaned.fastq.gz"
  HTML_REPORT="${OUTPUT_DIR}/reports/${SAMPLE_NAME}_report.html"
  JSON_REPORT="${OUTPUT_DIR}/reports/${SAMPLE_NAME}_report.json"
  
  # Log info tanpa mengganggu progress bar
  echo -e "\nMemproses sampel: $SAMPLE_NAME (${CURRENT_PAIR}/${TOTAL_PAIRS})" >> $LOG_FILE
  echo "  Input  R1: $R1_FILENAME" >> $LOG_FILE
  echo "  Input  R2: $R2_FILENAME" >> $LOG_FILE
  echo "  Output R1: $(basename "$OUT_R1")" >> $LOG_FILE
  echo "  Output R2: $(basename "$OUT_R2")" >> $LOG_FILE
  
  # Update status di layar (di atas progress bar)
  echo -ne "\rMemproses: $SAMPLE_NAME                                                          "
  
  # Jalankan fastp untuk paired-end dengan parameter yang dioptimalkan
  echo "  Menjalankan fastp..." >> $LOG_FILE
  fastp \
    -i "$R1_FILE" -I "$R2_FILE" \
    -o "$OUT_R1" -O "$OUT_R2" \
    -q $QUALITY \
    --cut_right \
    --cut_right_mean_quality $QUALITY \
    --adapter_sequence=$ADAPTER \
    --adapter_sequence_r2=$ADAPTER \
    -l $MIN_LENGTH \
    --compression 6 \
    --html "$HTML_REPORT" \
    --json "$JSON_REPORT" \
    --report_title "${SAMPLE_NAME} Quality Report" \
    --thread $THREADS \
    >> $LOG_FILE 2>&1
  
  # Update progress setelah fastp
  ((CURRENT_STEP++))
  PERCENT=$(( CURRENT_STEP * 100 / TOTAL_STEPS ))
  show_progress $PERCENT
  
  # Periksa apakah fastp berhasil
  if [ $? -eq 0 ]; then
    echo "  ✓ Fastp berhasil dijalankan" >> $LOG_FILE
  else
    echo -e "\n  ✗ ERROR menjalankan fastp untuk $SAMPLE_NAME" | tee -a $LOG_FILE
    echo "    Periksa file log untuk detail error" | tee -a $LOG_FILE
    # Tetap update progress
    ((CURRENT_STEP++))
    PERCENT=$(( CURRENT_STEP * 100 / TOTAL_STEPS ))
    show_progress $PERCENT
    continue
  fi
  
  # Jalankan FastQC pada file hasil untuk validasi kualitas
  echo "  Menjalankan FastQC untuk validasi..." >> $LOG_FILE
  fastqc "$OUT_R1" "$OUT_R2" -o "$OUTPUT_DIR/reports/" -t $THREADS >> $LOG_FILE 2>&1
  
  # Update progress setelah FastQC
  ((CURRENT_STEP++))
  PERCENT=$(( CURRENT_STEP * 100 / TOTAL_STEPS ))
  show_progress $PERCENT
  
  # Periksa apakah FastQC berhasil
  if [ $? -eq 0 ]; then
    echo "  ✓ FastQC berhasil dijalankan" >> $LOG_FILE
  else
    echo -e "\n  ✗ WARNING: FastQC gagal dijalankan, tetapi file output tetap dihasilkan" | tee -a $LOG_FILE
  fi
  
  # Ekstrak statistik penting dari laporan JSON fastp
  if [ -f "$JSON_REPORT" ]; then
    KEPT_READS=$(grep -o '"after_filtering":{"total_reads":[0-9]*' "$JSON_REPORT" | grep -o '[0-9]*$')
    if [ -n "$KEPT_READS" ]; then
      echo "  Statistik: $KEPT_READS reads dipertahankan setelah filtering" >> $LOG_FILE
    fi
    
    # Hitung persentase adapter yang ditemukan
    ADAPTER_RATE=$(grep -o '"adapter_trimmed_reads":[0-9]*' "$JSON_REPORT" | grep -o '[0-9]*$')
    TOTAL_READS=$(grep -o '"before_filtering":{"total_reads":[0-9]*' "$JSON_REPORT" | grep -o '[0-9]*$')
    
    if [ -n "$ADAPTER_RATE" ] && [ -n "$TOTAL_READS" ] && [ "$TOTAL_READS" -ne 0 ]; then
      PERCENT_ADAPTER=$(echo "scale=2; 100 * $ADAPTER_RATE / $TOTAL_READS" | bc)
      echo "  Statistik: $PERCENT_ADAPTER% reads mengandung adapter" >> $LOG_FILE
    fi
  else
    echo "  " | tee -a $LOG_FILE
  fi
  
  echo "  Selesai memproses: $SAMPLE_NAME" >> $LOG_FILE
  echo "-----------------------------------" >> $LOG_FILE
done

# Baris baru setelah progress bar selesai
echo -e "\n"

# Buat ringkasan MultiQC jika tersedia
if command -v multiqc &> /dev/null; then
  echo "Membuat laporan MultiQC..." | tee -a $LOG_FILE
  multiqc "$OUTPUT_DIR/reports/" -o "$OUTPUT_DIR/reports/" >> $LOG_FILE 2>&1
  if [ $? -eq 0 ]; then
    echo "✓ Laporan MultiQC berhasil dibuat" | tee -a $LOG_FILE
  else
    echo "✗ Gagal membuat laporan MultiQC" | tee -a $LOG_FILE
  fi
else
  echo "MultiQC tidak ditemukan. Melewati pembuatan laporan gabungan." | tee -a $LOG_FILE
fi

# Ringkasan akhir
echo "===================================" | tee -a $LOG_FILE
echo "Pemrosesan selesai pada: $(date)" | tee -a $LOG_FILE
echo "Total file yang berhasil diproses: $CURRENT_PAIR dari $TOTAL_PAIRS" | tee -a $LOG_FILE
echo "Laporan tersedia di: $OUTPUT_DIR/reports/" | tee -a $LOG_FILE
echo "File output tersedia di: $OUTPUT_DIR/" | tee -a $LOG_FILE
echo "Log lengkap tersedia di: $LOG_FILE" | tee -a $LOG_FILE
echo "===================================" | tee -a $LOG_FILE

# Tambahkan pesan tentang cara menjalankan analisis selanjutnya
echo ""
echo "Proses pembersihan FASTQ selesai. File yang sudah dibersihkan siap untuk analisis selanjutnya."
echo "Untuk menggunakan file tersebut dalam analisis, gunakan path: $OUTPUT_DIR/<nama_file>_cleaned.fastq.gz"