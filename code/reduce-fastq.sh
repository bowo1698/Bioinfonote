#!/bin/bash

# =========================================================================
# FASTQ-Reducer: Script untuk mengurangi ukuran file FASTQ.GZ
# Mendukung file paired-end dengan mempertahankan pasangan reads
# Version: 1.5 (Fixed Edition)
# Dependencies: seqtk, bc, parallel
# cara install:
# conda install -c bioconda seqtk
# conda install -c conda-forge bc
# conda install -c conda-forge parallel
# =========================================================================

set -o pipefail  # Menangkap error dalam pipeline

# Default values for parameters
INPUT_DIR="./raw_data"
OUTPUT_DIR="./reduced_data"
TARGET_SIZE=100
SEED=42
PREFIX="reduced_"
R1_PATTERN="_R1_"
R2_PATTERN="_R2_"
FILE_PATTERN="*.fastq.gz"
THREADS=1

# Function to display help message
show_help() {
    echo "Usage: $(basename "$0") [OPTIONS]"
    echo
    echo "Script untuk mengurangi ukuran file FASTQ.GZ paired-end menjadi ukuran target yang ditentukan."
    echo
    echo "OPTIONS:"
    echo "  -h, --help              Menampilkan pesan bantuan ini"
    echo "  -i, --input DIR         Direktori input yang berisi file FASTQ.GZ (default: $INPUT_DIR)"
    echo "  -o, --output DIR        Direktori output untuk file hasil (default: $OUTPUT_DIR)"
    echo "  -t, --target SIZE       Ukuran target dalam MB per file (default: $TARGET_SIZE)"
    echo "  -s, --seed NUMBER       Random seed untuk subsampling (default: $SEED)"
    echo "  -p, --prefix STRING     Prefiks untuk file output (default: $PREFIX)"
    echo "  -r1, --r1-pattern STR   Pattern untuk mengidentifikasi file R1 (default: $R1_PATTERN)"
    echo "  -r2, --r2-pattern STR   Pattern untuk mengidentifikasi file R2 (default: $R2_PATTERN)"
    echo "  -f, --file-pattern STR  Pattern untuk mencari file (default: $FILE_PATTERN)"
    echo "  -j, --threads NUMBER    Jumlah thread untuk digunakan (default: $THREADS)"
    echo
    echo "CONTOH:"
    echo "  $(basename "$0") --input /data/raw --output /data/reduced --target 150 --seed 123"
    echo "  $(basename "$0") -i /data/raw -o /data/reduced -t 150 -s 123 -p 'small_'"
    echo "  $(basename "$0") -i /data/raw -o /data/reduced -t 150 -s 123 -p 'small_' -r1 '_R1_' -r2 '_R2_'"
    echo
    echo "CATATAN:"
    echo "  - Script ini membutuhkan seqtk yang terinstal di sistem"
    echo "  - File akan diproses secara berpasangan (R1 dan R2) menggunakan pattern yang ditentukan"
    echo "  - Semua file hasil akan memiliki ukuran mendekati ukuran target yang ditentukan"
    echo
    exit 0
}

declare -a PARAMS=()

# Fungsi untuk mendapatkan ukuran file dalam bytes
get_file_size() {
    local file=$1
    if [ -f "$file" ]; then
        # Detect OS and use appropriate stat command
        if [[ "$OSTYPE" == "darwin"* ]]; then
            # macOS (BSD stat)
            stat -f%z "$file" 2>/dev/null || echo "0"
        else
            # Linux (GNU stat)
            stat -c%s "$file" 2>/dev/null || echo "0"
        fi
    else
        echo "0"
    fi
}

# Fungsi untuk konversi bytes ke MB
bytes_to_mb() {
    local bytes=$1
    echo $((bytes / 1024 / 1024))
}

# Fungsi untuk mendapatkan jumlah baris dalam file gzip
get_read_count() {
    local file=$1
    local count=0
    if [ -f "$file" ]; then
        # Gunakan zcat untuk files gzip dan hitung jumlah baris, dibagi 4 karena FASTQ
        count=$(zcat "$file" 2>/dev/null | head -n 40000 | wc -l)
        count=$((count / 4))
    fi
    echo "$count"
}

# Fungsi validasi file
validate_file_pair() {
    local r1_file=$1
    local r2_file=$2

    # Debug untuk ukuran file
    local r1_size_raw=$(get_file_size "$r1_file")
    local r2_size_raw=$(get_file_size "$r2_file")
    
    # Cek keberadaan file
    if [ ! -f "$r1_file" ]; then
        echo -e "\nPeringatan: File $r1_file tidak ditemukan. Melewati..."
        return 1
    fi
    
    if [ ! -f "$r2_file" ]; then
        echo -e "\nPeringatan: Pasangan file untuk $r1_file ($r2_file) tidak ditemukan. Melewati..."
        return 1
    fi
    
    # Cek ukuran file
    local r1_size=$(get_file_size "$r1_file")
    local r2_size=$(get_file_size "$r2_file")
    
    if [ "$r1_size" -eq 0 ] || [ "$r2_size" -eq 0 ]; then
        echo -e "\nPeringatan: File $r1_file atau $r2_file berukuran 0. Melewati..."
        return 1
    fi
    
    return 0
}

# Fungsi untuk memantau progress subsampling
monitor_subsampling() {
    local output_file=$1
    local input_file=$2
    local pid=$3
    local file_name=$4
    local target_size_mb=$5
    
    local delay=0.2
    local temp_file=$(mktemp -t progress_monitor.XXXXXX)
    local start_time=$(date +%s)
    
    # Dapatkan ukuran input file
    local input_file_size=$(get_file_size "$input_file")
    if [ "$input_file_size" -eq 0 ]; then
        echo -e "\r[$file_name] Error: Ukuran file input 0!"
        kill $pid 2>/dev/null
        return 1
    fi
    
    local target_bytes=$((target_size_mb * 1024 * 1024))
    local expected_ratio=$(echo "scale=4; $target_bytes / $input_file_size" | bc)
    # Format expected_ratio untuk memastikan ada 0 di depan nilai desimal
    expected_ratio=$(echo $expected_ratio | sed 's/^\./0./')
    
    # Use string comparison with bc instead of bash arithmetic
    if [ $(echo "$expected_ratio > 1" | bc -l) -eq 1 ]; then
        expected_ratio=0.99
    fi
    
    echo -ne "\r[$file_name] Memulai subsampling... [------------------------------] 0%"
    
    local prev_size=0
    local prev_time=$start_time
    local speed=0
    local eta_str="menghitung..."
    
    while kill -0 $pid 2>/dev/null; do
        sleep $delay
        [ ! -f "$output_file" ] && continue
        
        local current_size=$(get_file_size "$output_file")
        current_size=${current_size:-0}
        
        if [ "$current_size" -eq 0 ]; then
            continue
        fi

        local current_time=$(date +%s)
        local time_diff=$((current_time - prev_time))
        
        if [ $time_diff -ge 1 ]; then
            local size_diff=$((current_size - prev_size))
            speed=$((size_diff / 1024 / time_diff))
            prev_size=$current_size
            prev_time=$current_time
        fi
        
        local target_with_ratio=$(echo "$target_bytes * $expected_ratio" | bc)
        target_with_ratio=${target_with_ratio%.*}  # Ambil bagian integer
        # Jika 0 atau kosong, gunakan nilai minimal 1
        if [ -z "$target_with_ratio" ] || [ "$target_with_ratio" -eq 0 ]; then
            target_with_ratio=$target_bytes
            [ "$target_with_ratio" -eq 0 ] && target_with_ratio=1
        fi
        # Jika current_size 0, set progress ke 0, jika tidak hitung seperti biasa
        if [ "$current_size" -eq 0 ]; then
            local progress_by_size=0
        else
            local progress_by_size=$((current_size * 100 / target_with_ratio))
        fi
        [ "$progress_by_size" -gt 99 ] && progress_by_size=99
        
        # Perhitungan ETA
        if [ "$speed" -gt 0 ]; then
            local remaining_bytes=$(echo "$target_bytes * $expected_ratio - $current_size" | bc)
            remaining_bytes=${remaining_bytes%.*}  # Ambil bagian integer
            # Pastikan tidak negatif
            if [ -z "$remaining_bytes" ] || [ "$remaining_bytes" -lt 0 ]; then
                remaining_bytes=0
            fi
            local eta=$((remaining_bytes / 1024 / speed))
            if [ $eta -gt 60 ]; then
                eta_str="$((eta / 60))m $((eta % 60))s"
            else
                eta_str="${eta}s"
            fi
        fi
        
        # Tampilan progress bar
        local bar_length=30
        local completed=$((progress_by_size * bar_length / 100))
        printf -v bar "%${completed}s" | tr ' ' '#'
        printf -v empty "%$((bar_length - completed))s" | tr ' ' '-'

        # Tambahkan padding untuk membuat panjang output konsisten
        printf -v padding "%60s" " "
        echo -ne "\r[$file_name] [${bar}${empty}] ${progress_by_size}% (${speed} KB/s, ETA: ${eta_str})${padding}\r"
        echo -ne "\r[$file_name] [${bar}${empty}] ${progress_by_size}% (${speed} KB/s, ETA: ${eta_str})"
    done
    
    local elapsed=$(( $(date +%s) - start_time ))
    local elapsed_str=$( [ $elapsed -gt 60 ] && echo "$((elapsed / 60))m $((elapsed % 60))s" || echo "${elapsed}s" )
    echo -ne "\r[$file_name] [##############################] 100% (selesai dalam ${elapsed_str})      \n"
    
    rm -f "$temp_file"
}

# Fungsi utama pemrosesan file
process_pair() {
    local r1_file=$1
    local r2_file=$2
    local r1_output=$3
    local r2_output=$4
    local sampling_ratio=$5
    local index=$6
    local total=$7
    local log_file=$8
    local seed=$9
    local target_size=${10}
    
    # Hitung ukuran file
    local r1_size=$(bytes_to_mb $(get_file_size "$r1_file"))
    local r2_size=$(bytes_to_mb $(get_file_size "$r2_file"))
    local total_size=$((r1_size + r2_size))
    
    echo "[$index/$total] Memproses: $(basename "$r1_file") & $(basename "$r2_file")"
    echo "  Ukuran awal: ${total_size} MB, Rasio sampling: $(printf "%.4f" $sampling_ratio)"
    
    # Bersihkan file output lama
    rm -f "$r1_output" "$r2_output"
    
    # Proses R1
    seqtk sample -s$seed "$r1_file" $sampling_ratio | gzip > "$r1_output" &
    local r1_pid=$!
    monitor_subsampling "$r1_output" "$r1_file" $r1_pid "R1" $target_size
    
    if ! wait $r1_pid; then
        echo "Gagal memproses R1. Melewati pasangan ini..."
        rm -f "$r1_output" "$r2_output"
        return 1
    fi
    
    # Proses R2
    seqtk sample -s$seed "$r2_file" $sampling_ratio | gzip > "$r2_output" &
    local r2_pid=$!
    monitor_subsampling "$r2_output" "$r2_file" $r2_pid "R2" $target_size
    
    if ! wait $r2_pid; then
        echo "Gagal memproses R2. Melewati pasangan ini..."
        rm -f "$r1_output" "$r2_output"
        return 1
    fi
    
    # Hitung hasil
    local new_size=$(bytes_to_mb $(get_file_size "$r1_output"))
    new_size=$((new_size + $(bytes_to_mb $(get_file_size "$r2_output"))))
    local reduction_percent=$(echo "scale=2; (1 - $new_size/$total_size)*100" | bc)
    
    echo "  Selesai: Ukuran baru ${new_size} MB (${reduction_percent}% reduksi)"
    
    return 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help) show_help ;;
        -i|--input) INPUT_DIR="$2"; shift 2 ;;
        -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
        -t|--target) TARGET_SIZE="$2"; shift 2 ;;
        -s|--seed) SEED="$2"; shift 2 ;;
        -p|--prefix) PREFIX="$2"; shift 2 ;;
        -r1|--r1-pattern) R1_PATTERN="$2"; shift 2 ;;
        -r2|--r2-pattern) R2_PATTERN="$2"; shift 2 ;;
        -f|--file-pattern) FILE_PATTERN="$2"; shift 2 ;;
        -j|--threads) THREADS="$2"; shift 2 ;;
        *) echo "Parameter tidak dikenal: $1"; exit 1 ;;
    esac
done

# Validasi dependensi
for cmd in seqtk bc; do
    if ! command -v $cmd &> /dev/null; then
        echo "Error: $cmd tidak terinstal!"
        exit 1
    fi
done

# Validasi direktori
[ ! -d "$INPUT_DIR" ] && { echo "Direktori input tidak valid!"; exit 1; }
mkdir -p "$OUTPUT_DIR"

echo "Mencari file dengan pattern: *${R1_PATTERN}*${FILE_PATTERN} di $INPUT_DIR"

# Temukan file R1
R1_FILES=($(find "$INPUT_DIR" -name "*${R1_PATTERN}*${FILE_PATTERN}" | sort))
if [ ${#R1_FILES[@]} -eq 0 ]; then 
    echo "Tidak ada file yang ditemukan dengan pattern *${R1_PATTERN}*${FILE_PATTERN}!"
    echo "Mencoba pattern alternatif..."
    # Coba mencari dengan file pattern saja tanpa R1 pattern
    R1_FILES=($(find "$INPUT_DIR" -name "*${FILE_PATTERN}" | grep "${R1_PATTERN}" | sort))
    if [ ${#R1_FILES[@]} -eq 0 ]; then
        echo "Mencoba mencari file tanpa pattern spesifik..."
        R1_FILES=($(find "$INPUT_DIR" -name "*.fastq.gz" | grep -E "_1|_R1_|_R1" | sort))
        if [ ${#R1_FILES[@]} -eq 0 ]; then
            echo "Tidak dapat menemukan file FASTQ.GZ yang valid!"
            exit 1
        fi
    fi
fi

echo "Ditemukan ${#R1_FILES[@]} file yang cocok dengan pola R1"

# Debug: tampilkan semua file R1 yang ditemukan
echo "Daftar file R1 yang ditemukan:"
for R1_FILE in "${R1_FILES[@]}"; do
    echo "  - $R1_FILE"
    R2_FILE="${R1_FILE/$R1_PATTERN/$R2_PATTERN}"
    echo "    Pasangan: $R2_FILE"
    if [ -f "$R2_FILE" ]; then
        echo "    Pasangan ada"
    else
        echo "    Pasangan tidak ditemukan"
    fi
done

# Persiapan parameter
for ((i=0; i<${#R1_FILES[@]}; i++)); do
    R1_FILE="${R1_FILES[i]}"
    R2_FILE="${R1_FILE/$R1_PATTERN/$R2_PATTERN}"
    
    echo "Memeriksa pasangan: $R1_FILE dan $R2_FILE"
    
    if validate_file_pair "$R1_FILE" "$R2_FILE"; then
        # Hitung rasio sampling
        total_size_mb=$(( $(bytes_to_mb $(get_file_size "$R1_FILE")) + $(bytes_to_mb $(get_file_size "$R2_FILE")) ))
        sampling_ratio=$(echo "scale=4; (2*$TARGET_SIZE)/$total_size_mb" | bc | sed 's/^\./0./')
        [ $(echo "$sampling_ratio > 1" | bc -l) -eq 1 ] && sampling_ratio=1
        
        # Tambahkan ke parameter
        PARAMS+=("$R1_FILE" "$R2_FILE" \
            "$OUTPUT_DIR/${PREFIX}$(basename "$R1_FILE")" \
            "$OUTPUT_DIR/${PREFIX}$(basename "$R2_FILE")" \
            "$sampling_ratio" "$((i+1))" "${#R1_FILES[@]}" \
            "$OUTPUT_DIR/reduce_log.txt" "$SEED" "$TARGET_SIZE")
    else
        echo "Pasangan file tidak valid, mencoba metode alternatif pencarian pasangan..."
        # Coba cari pasangan dengan mengganti pola secara manual
        base_name=$(basename "$R1_FILE" | sed -E "s/${R1_PATTERN}.*//")
        R2_FILE=$(find "$INPUT_DIR" -name "${base_name}*${R2_PATTERN}*${FILE_PATTERN}" | head -n 1)
        
        if [ -n "$R2_FILE" ] && [ -f "$R2_FILE" ]; then
            echo "Menemukan pasangan alternatif: $R2_FILE"
            if validate_file_pair "$R1_FILE" "$R2_FILE"; then
                # Hitung rasio sampling
                total_size_mb=$(( $(bytes_to_mb $(get_file_size "$R1_FILE")) + $(bytes_to_mb $(get_file_size "$R2_FILE")) ))
                sampling_ratio=$(echo "scale=4; (2*$TARGET_SIZE)/$total_size_mb" | bc | sed 's/^\./0./')
                (( $(echo "$sampling_ratio > 1" | bc -l) )) && sampling_ratio=1
                
                # Tambahkan ke parameter
                PARAMS+=("$R1_FILE" "$R2_FILE" \
                    "$OUTPUT_DIR/${PREFIX}$(basename "$R1_FILE")" \
                    "$OUTPUT_DIR/${PREFIX}$(basename "$R2_FILE")" \
                    "$sampling_ratio" "$((i+1))" "${#R1_FILES[@]}" \
                    "$OUTPUT_DIR/reduce_log.txt" "$SEED" "$TARGET_SIZE")
            fi
        fi
    fi
done

if [ ${#PARAMS[@]} -eq 0 ]; then
    echo "Tidak ada file valid! Periksa pola pencarian R1/R2 atau file yang ada."
    exit 1
fi

echo "Ditemukan ${#PARAMS[@]} pasang file yang valid untuk diproses"

# Implementasi pemrosesan paralel langsung tanpa wrapper
process_pair_parallel() {
    local args=("$@")
    process_pair "${args[@]}"
}
export -f process_pair
export -f process_pair_parallel
export -f get_file_size
export -f bytes_to_mb
export -f monitor_subsampling
export -f get_read_count
export -f validate_file_pair

# Eksekusi pemrosesan
if (( THREADS > 1 )) && command -v parallel &>/dev/null; then
    echo "Memproses dengan $THREADS thread menggunakan parallel..."
    parallel -j $THREADS -n10 process_pair_parallel ::: "${PARAMS[@]}"
else
    echo "Memproses secara sekuensial..."
    for ((i=0; i<${#PARAMS[@]}; i+=10)); do
        r1_file="${PARAMS[$i]}"
        r2_file="${PARAMS[$i+1]}"
        r1_output="${PARAMS[$i+2]}"
        r2_output="${PARAMS[$i+3]}"
        sampling_ratio="${PARAMS[$i+4]}"
        index="${PARAMS[$i+5]}"
        total="${PARAMS[$i+6]}"
        log_file="${PARAMS[$i+7]}"
        seed="${PARAMS[$i+8]}"
        target_size="${PARAMS[$i+9]}"
        
        process_pair "$r1_file" "$r2_file" "$r1_output" "$r2_output" "$sampling_ratio" "$index" "$total" "$log_file" "$seed" "$target_size"
    done
fi

echo "Proses selesai! Hasil disimpan di: $OUTPUT_DIR"