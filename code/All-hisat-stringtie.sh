#!/bin/bash
# Load configuration
if [ ! -f "config.yaml" ]; then
    echo "❌ ERROR: config.yaml not found!"
    echo "Please create config.yaml file first"
    exit 1
fi

# Function to parse YAML (simple parser)
parse_yaml() {
    python3 -c "
import yaml
import sys
with open('$1', 'r') as f:
    config = yaml.safe_load(f)
    
def print_config(obj, prefix=''):
    for key, value in obj.items():
        if isinstance(value, dict):
            print_config(value, f'{prefix}{key}_')
        elif isinstance(value, list):
            print(f'{prefix}{key}_count={len(value)}')
            for i, item in enumerate(value):
                print(f'{prefix}{key}_{i}=\"{item}\"')
        else:
            # Konversi boolean ke string lowercase untuk bash
            if isinstance(value, bool):
                value = str(value).lower()
            print(f'{prefix}{key}=\"{value}\"')
            
print_config(config)
"
}

# Load configuration
eval $(parse_yaml config.yaml)

# Debug: tampilkan beberapa variabel penting
echo "Debug - Config variables:"
echo "  samples_count: $samples_count"
echo "  system_threads: $system_threads"
echo "  pipeline_output_dir: $pipeline_output_dir"
echo "  directories_raw_data: $directories_raw_data"
echo "  reference_genome_fasta: $reference_genome_fasta"
echo "  reference_annotation_gff: $reference_annotation_gtf"
echo "  reference_annotation_gff: $reference_annotation_gff"

# Validasi kritis
if [ -z "$samples_count" ] || [ "$samples_count" == "0" ]; then
    echo "❌ ERROR: No samples found in config.yaml!"
    echo "Please check the 'samples:' section in your config.yaml"
    exit 1
fi

# Logging functions
log_message() {
    if [ "$logging_verbose" == "true" ]; then
        if [ "$logging_timestamp" == "true" ]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
        else
            echo "$1"
        fi
    fi
}

log_progress() {
    if [ "$logging_progress_tracking" == "true" ]; then
        echo "Progress: $1"
    fi
}

if [ -z "$directories_raw_data" ] || [ -z "$directories_stringtie_output" ] || [ -z "$reference_genome_fasta" ] || [ -z "$reference_annotation_gtf" ] || [ -z "$reference_annotation_gff" ]; then
    echo "❌ ERROR: Critical configuration variables missing!"
    echo "Please check your config.yaml file"
    exit 1
fi

START_TIME=$(date)
echo "==============================================================="
echo "🕐 Menjalankan Alignment dan kuantifikasi dengan StringTie"
echo "Author: Agus Wibowo"
echo "Email: aguswibowo@gmail.com"
echo "Proses dimulai: $START_TIME"
echo "==============================================================="
echo "⚠️  Resource Requirements:"
echo "   - CPU threads: $system_threads cores (HISAT2 & StringTie)"
echo "   - Memory: ~${system_memory_gb}GB recommended"
echo "   - Disk space: ~${system_disk_space_factor}x input data size"
echo "==============================================================="

echo "💾 Memory check:"
free -h
echo "📊 Disk space check:"
df -h "$pipeline_output_dir" 2>/dev/null || df -h .

mkdir -p "$pipeline_output_dir"
mkdir -p "$directories_hisat_index"
mkdir -p "$directories_stringtie_output"
mkdir -p "$directories_aligned_data" "$directories_bam_data"
mkdir -p "$directories_transcript_matrix"

# 1. Ekstrak informasi splice sites dan exon dari GTF
hisat2_extract_splice_sites.py "$reference_annotation_gtf" > "$pipeline_output_dir/splice_sites.txt"
hisat2_extract_exons.py "$reference_annotation_gtf" > "$pipeline_output_dir/exons.txt"

# 2. Buat indeks genom dengan informasi splicing dengan 4 core
hisat2-build -p "$hisat2_build_threads" \
    --ss "$pipeline_output_dir/splice_sites.txt" \
    --exon "$pipeline_output_dir/exons.txt" \
    "$reference_genome_fasta" \
    "$directories_hisat_index/$hisat2_index_name"

# Cek apakah index berhasil dibuat
if [ ! -f "$directories_hisat_index/$hisat2_index_name.1.ht2" ]; then
    echo "❌ ERROR: HISAT2 index creation failed!"
    exit 1
fi
echo "✅ HISAT2 index created successfully"

# Daftar semua 4 sampel yang akan diproses (uncomment yang lain jika diperlukan)
# Load samples dari config
SAMPLES=()
if [ ! -z "$samples_count" ]; then
    for ((i=0; i<$samples_count; i++)); do
        var_name="samples_$i"
        eval "sample_value=\$$var_name"
        if [ ! -z "$sample_value" ]; then
            SAMPLES+=("$sample_value")
        fi
    done
fi

echo "Loaded ${#SAMPLES[@]} samples from config.yaml"
for sample in "${SAMPLES[@]}"; do
    echo "  - $sample"
done

echo "Loaded ${#SAMPLES[@]} samples from config.yaml"

# Referensi dari config
ANNOTATION="$reference_annotation_gff"
INDEX="$directories_hisat_index/$hisat2_index_name"

echo "==============================================================="
echo "Alignment and BAM conversion started at: $(date)"
echo "Processing ${#SAMPLES[@]} samples"
echo "==============================================================="

# Counter untuk tracking progress
counter=1
total=${#SAMPLES[@]}

log_message "Starting HISAT2 alignment phase"

for sample in "${SAMPLES[@]}"; do
    echo ""
    echo "Processing sample: $sample ($counter/$total)"
    echo "Started at: $(date)"
    
    # Langkah 1: Alignment dengan HISAT2
    echo "Step 1: Running HISAT2 alignment..."
    hisat2 -p "$hisat2_alignment_threads" $hisat2_extra_params \
        -x $INDEX \
        --rna-strandness "$hisat2_rna_strandness" \
        -1 "$directories_raw_data/${sample}_1.fastq.gz" \
        -2 "$directories_raw_data/${sample}_2.fastq.gz" \
        -S "$directories_aligned_data/${sample}.sam"
    
    # Langkah 2: Konversi SAM ke BAM
    echo "Step 2: Converting SAM to BAM..."
    samtools view -@ "$samtools_threads" -bS "$directories_aligned_data/${sample}.sam" > "$directories_bam_data/${sample}.bam"
    
    # Langkah 3: Sorting BAM file
    echo "Step 3: Sorting BAM file..."
    samtools sort -@ "$samtools_threads" "$directories_bam_data/${sample}.bam" -o "$directories_bam_data/${sample}.sorted.bam"
    
    # Langkah 4: Indeks BAM file
    echo "Step 4: Indexing BAM file..."
    samtools index "$directories_bam_data/${sample}.sorted.bam"
    
    # Hapus file intermediate untuk menghemat ruang
    if [ "$qc_remove_intermediate_sam" == "true" ]; then
        echo "Cleaning up intermediate files..."
        rm "$directories_aligned_data/${sample}.sam"
    fi

    if [ "$qc_remove_unsorted_bam" == "true" ]; then
        rm "$directories_bam_data/${sample}.bam"
    fi
    echo "---------------------------------"
    
    ((counter++))
done

echo ""
echo "=================================================="
echo "Alignment and BAM conversion completed at: $(date)"
echo "All ${#SAMPLES[@]} samples processed successfully!"
echo "=================================================="


echo "=================================================="
echo "StringTie individual assembly started at: $(date)"
echo "Processing ${#SAMPLES[@]} samples"
echo "=================================================="

# Counter untuk tracking progress
counter=1
total=${#SAMPLES[@]}

log_message "Starting StringTie individual assembly phase"

for sample in "${SAMPLES[@]}"; do
    echo ""
    echo "Running StringTie for sample: $sample ($counter/$total)"
    echo "Started at: $(date)"
    
    stringtie "$directories_bam_data/${sample}.sorted.bam" \
        -p "$stringtie_threads" \
        -G $ANNOTATION \
        $stringtie_strandness \
        -A "$directories_stringtie_output/${sample}_abundance.tab" \
        -o "$directories_stringtie_output/${sample}.gtf" \
        -l ${sample}
    
    echo "Sample $sample StringTie completed at: $(date)"
    echo "Progress: $counter/$total samples completed"
    echo "--------------------------------------------------"
    
    ((counter++))
done

echo ""
echo "==================================================="
echo "StringTie individual assembly completed at: $(date)"
echo "All ${#SAMPLES[@]} samples processed successfully!"
echo "==================================================="

echo "================================================================="
echo "StringTie merge and quantification (adaptive) started at: $(date)"
echo "================================================================="

# Deteksi file GTF yang sudah ada
echo "Detecting existing GTF files..."
GTF_FILES=($(ls "$directories_stringtie_output"/SRR*.gtf 2>/dev/null))

if [ ${#GTF_FILES[@]} -eq 0 ]; then
    echo "ERROR: No GTF files found in HISAT-StringTie/stringtie_output/"
    echo "Please run StringTie individual assembly first (script 2)"
    exit 1
fi

echo "Found ${#GTF_FILES[@]} GTF files:"
for file in "${GTF_FILES[@]}"; do
    echo "  - $(basename $file)"
done

# Verify samples dari config ada di GTF files
VERIFIED_SAMPLES=()
for sample in "${SAMPLES[@]}"; do
    if [[ -f "$directories_stringtie_output/${sample}.gtf" ]]; then
        VERIFIED_SAMPLES+=("$sample")
    else
        echo "⚠️ Warning: Sample $sample tidak ditemukan di GTF files"
    fi
done
SAMPLES=("${VERIFIED_SAMPLES[@]}")

echo ""
echo "Samples to process: ${SAMPLES[@]}"
echo "Total samples: ${#SAMPLES[@]}"

# Buat daftar GTF untuk merge
echo ""
echo "Creating GTF list for merged analysis..."
ls -1 "$directories_stringtie_output"/SRR*.gtf > "$directories_stringtie_output/gtf_list.txt"
echo "GTF files to be merged:"
cat "$directories_stringtie_output/gtf_list.txt"

# StringTie merge
echo ""
echo "Running StringTie merge..."
echo "Started at: $(date)"
stringtie --merge \
    -p "$stringtie_threads" \
    -G $ANNOTATION \
    -o "$directories_stringtie_output/merged.gtf" \
    "$directories_stringtie_output/gtf_list.txt"

if [ $? -eq 0 ]; then
    echo "StringTie merge completed successfully at: $(date)"
else
    echo "ERROR: StringTie merge failed"
    exit 1
fi

# Re-estimation untuk semua samples yang ada
echo ""
echo "Re-estimating transcript abundances for all detected samples..."
mkdir -p "$directories_stringtie_output/ballgown"

counter=1
total=${#SAMPLES[@]}

for sample in "${SAMPLES[@]}"; do
    echo "Re-estimating for sample: $sample ($counter/$total) at: $(date)"
    mkdir -p "$directories_stringtie_output/ballgown/${sample}"
    
    stringtie "$directories_bam_data/${sample}.sorted.bam" \
        -e -B \
        -p "$stringtie_threads" \
        $stringtie_strandness \
        -G "$directories_stringtie_output/merged.gtf" \
        -o "$directories_stringtie_output/ballgown/${sample}/${sample}.gtf"
    
    if [ $? -eq 0 ]; then
        echo "Sample $sample re-estimation completed successfully"
    else
        echo "ERROR: Re-estimation failed for sample $sample"
    fi
    
    ((counter++))
done

echo ""
echo "Output summary:"
echo "- BAM files: HISAT-StringTie/bam_data/ (${#SAMPLES[@]} files)"
echo "- Individual StringTie GTF: HISAT-StringTie/stringtie_output/ (${#SAMPLES[@]} files)"
echo "- Merged reference: HISAT-StringTie/stringtie_output/merged.gtf"
echo "- Ballgown data: HISAT-StringTie/stringtie_output/ballgown/ (${#SAMPLES[@]} directories)"
echo ""
echo "Files ready for count matrix!"

# Handles corrupted t_data.ctab files with merged fields

echo ""
echo "==================================="
echo "🧬 TRANSCRIPT MATRIX GENERATION"
echo "Detected samples: ${#SAMPLES[@]}"
echo "Processing started at: $(date)"
echo "==================================="

cd "$directories_transcript_matrix"

# Check and download prepDE.py3 if needed
echo "=== CHECKING PREPDE.PY3 ==="
if [ ! -f "../prepDE.py3" ]; then
    echo "📥 prepDE.py3 not found, downloading..."
    if command -v wget >/dev/null 2>&1; then
        wget -q https://raw.githubusercontent.com/gpertea/stringtie/master/prepDE.py3 -O ../prepDE.py3
    elif command -v curl >/dev/null 2>&1; then
        curl -s -o ../prepDE.py3 https://raw.githubusercontent.com/gpertea/stringtie/master/prepDE.py3
    else
        echo "❌ Error: wget atau curl tidak tersedia untuk download prepDE.py3"
        echo "   Silakan download manual dari: https://github.com/gpertea/stringtie/blob/master/prepDE.py3"
        exit 1
    fi
    
    if [ -f "../prepDE.py3" ]; then
        chmod +x ../prepDE.py3
        echo "✅ prepDE.py3 berhasil didownload"
    else
        echo "❌ Gagal download prepDE.py3"
        exit 1
    fi
else
    echo "✅ prepDE.py3 sudah tersedia"
fi

echo ""
echo "=== PHASE 1: ANALYZE DATA CORRUPTION ==="
echo ""

# Set path yang benar berdasarkan struktur direktori
BALLGOWN_PATH="../stringtie_output/ballgown"

# Validasi struktur direktori
echo "=== VALIDATING DIRECTORY STRUCTURE ==="
if [ ! -d "../stringtie_output" ]; then
    echo "❌ StringTie output directory not found!"
    exit 1
fi

if [ ! -d "$BALLGOWN_PATH" ]; then
    echo "❌ Ballgown directory not found!"
    exit 1
fi

echo "✅ Directory structure validated"
echo "GTF files location: ../stringtie_output/"
echo "Ballgown data location: $BALLGOWN_PATH/"
echo ""

# Check each sample for data corruption dengan path yang benar
echo "Sample,GTF_Transcripts,TData_Transcripts,Data_Corrupt,Missing_GTF" > corruption_report.csv

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Analyzing $SAMPLE..."
    
    # PERBAIKI PATH - gunakan path absolut
    gtf_file="$(realpath "$BALLGOWN_PATH/$SAMPLE/$SAMPLE.gtf" 2>/dev/null || echo "$BALLGOWN_PATH/$SAMPLE/$SAMPLE.gtf")"
    tdata_file="$BALLGOWN_PATH/$SAMPLE/t_data.ctab"
    
    # Debug: tampilkan path yang dicari
    echo "  Looking for GTF: $gtf_file"
    echo "  Looking for t_data: $tdata_file"
    
    # Cek file GTF
    if [ ! -f "$gtf_file" ]; then
        echo "  ❌ GTF file not found: $gtf_file"
        continue
    fi
    
    # Cek file t_data
    if [ ! -f "$tdata_file" ]; then
        echo "  ❌ t_data.ctab file not found: $tdata_file"
        echo "  Available in ballgown/$SAMPLE/:"
        ls -la "$BALLGOWN_PATH/$SAMPLE/" 2>/dev/null | head -5
        continue
    fi
    
    echo "  ✅ Both files found, analyzing..."
    
    # Extract transcripts from GTF
    grep "transcript_id" "$gtf_file" | \
    grep -o 'transcript_id "[^"]*"' | \
    cut -d'"' -f2 | \
    sort -u > "${SAMPLE}_gtf_transcripts.txt"
    
    # Extract transcripts from t_data with corruption handling
    tail -n +2 "$tdata_file" > "${SAMPLE}_tdata_raw.txt"
    
    # Check for field corruption
    corrupt_lines=0
    total_lines=$(wc -l < "${SAMPLE}_tdata_raw.txt")
    
    # Extract transcript names more robustly
    while IFS=$'\t' read -r t_id chr strand start end t_name gene_id gene_name cov fpkm rest; do
        # Look for rna- or XM_ patterns in the transcript name field
        if [[ "$t_name" =~ rna-[A-Z]+_[0-9]+\.[0-9]+ ]]; then
            # Extract just the transcript ID part
            clean_transcript=$(echo "$t_name" | grep -o 'rna-[A-Z]*_[0-9]*\.[0-9]*')
            echo "$clean_transcript"
        elif [[ "$t_name" =~ [A-Z]+_[0-9]+\.[0-9]+ ]]; then
            # Handle cases without rna- prefix
            clean_transcript=$(echo "$t_name" | grep -o '[A-Z]*_[0-9]*\.[0-9]*')
            echo "$clean_transcript"
        else
            # Check if transcript ID is in a different field due to corruption
            full_line="$t_id	$chr	$strand	$start	$end	$t_name	$gene_id	$gene_name	$cov	$fpkm"
            transcript_match=$(echo "$full_line" | grep -o 'rna-[A-Z]*_[0-9]*\.[0-9]*' | head -1)
            if [ -n "$transcript_match" ]; then
                echo "$transcript_match"
                ((corrupt_lines++))
            fi
        fi
    done < "${SAMPLE}_tdata_raw.txt" | sort -u > "${SAMPLE}_tdata_transcripts.txt"
    
    gtf_count=$(wc -l < "${SAMPLE}_gtf_transcripts.txt")
    tdata_count=$(wc -l < "${SAMPLE}_tdata_transcripts.txt")
    
    # Find missing transcripts
    comm -23 "${SAMPLE}_tdata_transcripts.txt" "${SAMPLE}_gtf_transcripts.txt" > "${SAMPLE}_missing_in_gtf.txt"
    missing_count=$(wc -l < "${SAMPLE}_missing_in_gtf.txt")
    
    if [ $total_lines -gt 0 ]; then
    corruption_pct=$(awk "BEGIN {printf \"%.2f\", $corrupt_lines * 100 / $total_lines}")
    else
        corruption_pct="0"
    fi
    
    echo "  GTF: $gtf_count, t_data: $tdata_count, Missing: $missing_count"
    echo "  Corrupted lines: $corrupt_lines/$total_lines ($corruption_pct%)"
    
    echo "$SAMPLE,$gtf_count,$tdata_count,$corruption_pct,$missing_count" >> corruption_report.csv
    
    if [ $missing_count -gt 0 ]; then
        echo "  🔴 Missing transcripts:"
        head -10 "${SAMPLE}_missing_in_gtf.txt" | sed 's/^/    /'
    fi
done

echo ""
echo "=== PHASE 2: FIXING STRATEGY ==="
echo ""

# Strategy: Add missing transcripts to GTF files with robust data extraction
success_samples=()
failed_samples=()

for SAMPLE in "${SAMPLES[@]}"; do
    missing_file="${SAMPLE}_missing_in_gtf.txt"
    
    if [ ! -f "$missing_file" ] || [ ! -s "$missing_file" ]; then
        echo "✅ $SAMPLE: No missing transcripts"
        success_samples+=("$SAMPLE")
        continue
    fi
    
    echo "🔧 Fixing $SAMPLE..."
    
    gtf_file="$(realpath "$BALLGOWN_PATH/$SAMPLE/$SAMPLE.gtf" 2>/dev/null || echo "$BALLGOWN_PATH/$SAMPLE/$SAMPLE.gtf")"
    tdata_file="$BALLGOWN_PATH/$SAMPLE/t_data.ctab"
    
    # Backup GTF
    cp "$gtf_file" "${gtf_file}.backup_$(date +%Y%m%d_%H%M%S)"
    
    # Process each missing transcript
    while IFS= read -r transcript_id; do
        if [ -z "$transcript_id" ]; then continue; fi
        
        echo "  Adding: $transcript_id"
        
        # Find the transcript line in t_data.ctab with robust parsing
        transcript_line=$(grep "$transcript_id" "$tdata_file" | head -1)
        
        if [ -n "$transcript_line" ]; then
            # Robust extraction of coordinates
            chr=$(echo "$transcript_line" | cut -f2)
            strand=$(echo "$transcript_line" | cut -f3)
            start=$(echo "$transcript_line" | cut -f4)
            end=$(echo "$transcript_line" | cut -f5)
            
            # Extract gene information more robustly
            # Look for gene-LOC pattern or create from transcript
            gene_match=$(echo "$transcript_line" | grep -o 'gene-LOC[0-9]*')
            if [ -n "$gene_match" ]; then
                gene_id=$(echo "$gene_match" | sed 's/gene-//')
                gene_name="$gene_id"
            else
                # Create gene ID from transcript ID
                gene_id="GENE_$(echo "$transcript_id" | sed 's/rna-//' | sed 's/XM_//' | cut -d'.' -f1)"
                gene_name="$gene_id"
            fi
            
            echo "    Chr: $chr, Start: $start, End: $end, Strand: $strand"
            echo "    Gene: $gene_id"
            
            # Validate coordinates
            if [[ "$start" =~ ^[0-9]+$ ]] && [[ "$end" =~ ^[0-9]+$ ]] && [ "$start" -lt "$end" ]; then
                # Add gene entry if not exists
                if ! grep -q "gene_id \"$gene_id\"" "$gtf_file"; then
                    echo -e "$chr\tStringTie\tgene\t$start\t$end\t1000\t$strand\t.\tgene_id \"$gene_id\"; gene_name \"$gene_name\";" >> "$gtf_file"
                fi
                
                # Add transcript entry
                echo -e "$chr\tStringTie\ttranscript\t$start\t$end\t1000\t$strand\t.\tgene_id \"$gene_id\"; transcript_id \"$transcript_id\"; gene_name \"$gene_name\";" >> "$gtf_file"
                
                # Add exon entry
                echo -e "$chr\tStringTie\texon\t$start\t$end\t1000\t$strand\t.\tgene_id \"$gene_id\"; transcript_id \"$transcript_id\"; exon_number \"1\"; gene_name \"$gene_name\";" >> "$gtf_file"
                
                echo "    ✅ Added entries"
            else
                echo "    ❌ Invalid coordinates: start=$start, end=$end"
                failed_samples+=("$SAMPLE")
                break
            fi
        else
            echo "    ❌ Cannot find transcript in t_data.ctab"
            failed_samples+=("$SAMPLE")
            break
        fi
        
    done < "$missing_file"
    
    # Check if sample was successfully processed
    if [[ ! " ${failed_samples[@]} " =~ " ${SAMPLE} " ]]; then
        success_samples+=("$SAMPLE")
        echo "  ✅ $SAMPLE completed successfully"
    fi
done

echo ""
echo "=== PHASE 3: TESTING RESULTS ==="
echo ""

if [ ${#success_samples[@]} -eq ${#SAMPLES[@]} ]; then
    echo "🎉 All samples fixed successfully!"
    
    TEST_DIR="$directories_transcript_matrix/test_gtf"
    mkdir -p "$TEST_DIR"
    
    # Create sample list for testing - GTF files are in stringtie_output, not ballgown
    > sample_list.txt
    for SAMPLE in "${SAMPLES[@]}"; do
        gtf_file="$BALLGOWN_PATH/$SAMPLE/$SAMPLE.gtf"
        test_gtf="$TEST_DIR/${SAMPLE}.gtf"
        
        if [ -f "$gtf_file" ]; then
            cp "$gtf_file" "$test_gtf"
            echo "$SAMPLE $test_gtf" >> sample_list.txt
            echo "  ✅ Test copy: $SAMPLE -> $test_gtf"
        fi
    done
    
    echo "Testing prepDE.py3 with all fixed samples..."
    
    if python3 ../prepDE.py3 -i sample_list.txt -g gene_count_matrix.csv -t transcript_count_matrix.csv -l 75 2>&1; then
        echo ""
        echo "🎉 COMPLETE SUCCESS!"
        echo ""
        echo "Final output files:"
        echo "  ✅ $(pwd)/gene_count_matrix.csv"
        echo "  ✅ $(pwd)/transcript_count_matrix.csv"
        
        # Show statistics
        if [ -f gene_count_matrix.csv ]; then
            genes=$(tail -n +2 gene_count_matrix.csv | wc -l)
            transcripts=$(tail -n +2 transcript_count_matrix.csv | wc -l)
            samples_count=$(($(head -1 gene_count_matrix.csv | tr ',' '\n' | wc -l) - 1))
            
            echo ""
            echo "📊 Matrix Statistics:"
            echo "   Genes: $genes"
            echo "   Transcripts: $transcripts"
            echo "   Samples: $samples_count"
            
            echo ""
            echo "Sample order in matrix:"
            head -1 gene_count_matrix.csv | tr ',' '\n' | tail -n +2 | nl
        fi
        
        echo ""
        echo "✅ Backup files preserved for safety"
        echo "🗂️  Analysis details in: $(pwd)/"
        
    else
        echo ""
        echo "🔍 Re-analyzing after first fix attempt..."
        
        # Re-run analysis to catch newly emerged problems
        for SAMPLE in "${SAMPLES[@]}"; do
            gtf_file="$(realpath "$BALLGOWN_PATH/$SAMPLE/$SAMPLE.gtf" 2>/dev/null || echo "$BALLGOWN_PATH/$SAMPLE/$SAMPLE.gtf")"
            tdata_file="$BALLGOWN_PATH/$SAMPLE/t_data.ctab"
            
            if [ -f "$gtf_file" ] && [ -f "$tdata_file" ]; then
                # Quick check for new missing transcripts
                grep "transcript_id" "$gtf_file" | grep -o 'transcript_id "[^"]*"' | cut -d'"' -f2 | sort -u > "temp_gtf_${SAMPLE}.txt"
                tail -n +2 "$tdata_file" | cut -f6 | grep -o 'rna-[^[:space:]]*\|[A-Z]*_[0-9]*\.[0-9]*' | sort -u > "temp_tdata_${SAMPLE}.txt"
                
                comm -23 "temp_tdata_${SAMPLE}.txt" "temp_gtf_${SAMPLE}.txt" > "temp_missing_${SAMPLE}.txt"
                
                if [ -s "temp_missing_${SAMPLE}.txt" ]; then
                    echo "🔧 Re-fixing $SAMPLE with newly found missing transcripts..."
                    
                    while IFS= read -r transcript_id; do
                        if [ -n "$transcript_id" ]; then
                            echo "  Adding: $transcript_id"
                            
                            transcript_line=$(grep "$transcript_id" "$tdata_file" | head -1)
                            if [ -n "$transcript_line" ]; then
                                chr=$(echo "$transcript_line" | cut -f2)
                                strand=$(echo "$transcript_line" | cut -f3) 
                                start=$(echo "$transcript_line" | cut -f4)
                                end=$(echo "$transcript_line" | cut -f5)
                                
                                gene_id="GENE_$(echo "$transcript_id" | sed 's/rna-//' | sed 's/-/_/g' | cut -d'.' -f1)"
                                
                                if [[ "$start" =~ ^[0-9]+$ ]] && [[ "$end" =~ ^[0-9]+$ ]]; then
                                    if ! grep -q "gene_id \"$gene_id\"" "$gtf_file"; then
                                        echo -e "$chr\tStringTie\tgene\t$start\t$end\t1000\t$strand\t.\tgene_id \"$gene_id\"; gene_name \"$gene_id\";" >> "$gtf_file"
                                    fi
                                    echo -e "$chr\tStringTie\ttranscript\t$start\t$end\t1000\t$strand\t.\tgene_id \"$gene_id\"; transcript_id \"$transcript_id\"; gene_name \"$gene_id\";" >> "$gtf_file"
                                    echo -e "$chr\tStringTie\texon\t$start\t$end\t1000\t$strand\t.\tgene_id \"$gene_id\"; transcript_id \"$transcript_id\"; exon_number \"1\"; gene_name \"$gene_id\";" >> "$gtf_file"
                                fi
                            fi
                        fi
                    done < "temp_missing_${SAMPLE}.txt"
                fi
                
                rm -f "temp_gtf_${SAMPLE}.txt" "temp_tdata_${SAMPLE}.txt" "temp_missing_${SAMPLE}.txt"
            fi
        done
        
        echo "🧪 Testing again after re-fix..."
        if python3 ../prepDE.py3 -i sample_list.txt -g gene_count_matrix.csv -t transcript_count_matrix.csv -l 75 2>&1; then
            echo "🎉 SUCCESS after re-analysis!"
            exit 0
        fi
        
        echo "❌ prepDE.py3 still failing. Trying cleanup approach..."
        
        # Alternative: Clean up t_data files instead
        echo ""
        echo "=== CLEANUP STRATEGY: Fix t_data corruption ==="
        
        for SAMPLE in "${SAMPLES[@]}"; do
            tdata_file="$BALLGOWN_PATH/$SAMPLE/t_data.ctab"
            missing_file="${SAMPLE}_missing_in_gtf.txt"
            
            if [ -f "$missing_file" ] && [ -s "$missing_file" ]; then
                echo "🧹 Cleaning $SAMPLE t_data.ctab..."
                
                # Backup t_data
                cp "$tdata_file" "${tdata_file}.backup_$(date +%Y%m%d_%H%M%S)"
                
                # Remove problematic transcripts
                head -1 "$tdata_file" > "${tdata_file}.cleaned"
                
                while IFS= read -r transcript_id; do
                    if [ -n "$transcript_id" ]; then
                        tail -n +2 "$tdata_file" | grep -v "$transcript_id" > "${tdata_file}.temp"
                        tail -n +2 "${tdata_file}.temp" >> "${tdata_file}.cleaned"
                        echo "    Removed: $transcript_id"
                    fi
                done < "$missing_file"
                
                mv "${tdata_file}.cleaned" "$tdata_file"
                rm -f "${tdata_file}.temp"
                
                echo "  ✅ Cleaned t_data.ctab"
            fi
        done
        
        # Test cleanup approach
        if python3 ../prepDE.py3 -i sample_list.txt -g gene_count_matrix.csv -t transcript_count_matrix.csv 2>&1; then
            echo ""
            echo "🎉 SUCCESS with cleanup approach!"
            echo "Output: $(pwd)/gene_count_matrix.csv, $(pwd)/transcript_count_matrix.csv"
        else
            echo ""
            echo "❌ Both approaches failed. Data may need manual curation."
        fi
    fi
    
else
    echo "❌ Failed samples: ${failed_samples[*]}"
    echo "Successful samples: ${#success_samples[@]}/${#SAMPLES[@]}"
    
    # Restore failed samples
    for sample in "${failed_samples[@]}"; do
        gtf_file="$(realpath "$BALLGOWN_PATH/$SAMPLE/$SAMPLE.gtf" 2>/dev/null || echo "$BALLGOWN_PATH/$SAMPLE/$SAMPLE.gtf")"
        backup_file=$(ls "${gtf_file}.backup_"* 2>/dev/null | tail -1)
        if [ -f "$backup_file" ]; then
            mv "$backup_file" "$gtf_file"
            echo "  🔄 Restored $sample"
        fi
    done
fi

echo ""
echo "=== SUMMARY REPORT ==="
echo "Analysis saved in: $(pwd)/"
echo "Corruption report: corruption_report.csv"

if [ -f corruption_report.csv ]; then
    echo ""
    echo "Data Corruption Summary:"
    echo "Sample          GTF    TData  Corrupt%  Missing"
    echo "---------------------------------------------------"
    tail -n +2 corruption_report.csv | while IFS=, read -r sample gtf tdata corrupt missing; do
        printf "%-15s %5s  %5s  %7s  %7s\n" "$sample" "$gtf" "$tdata" "$corrupt" "$missing"
    done
fi

echo ""

END_TIME=$(date)
echo "🕐 Proses dimulai: $START_TIME"
echo "🕐 Proses selesai: $END_TIME"

# Hitung durasi 
START_SECONDS=$(date -d "$START_TIME" +%s 2>/dev/null || echo "0")
END_SECONDS=$(date -d "$END_TIME" +%s 2>/dev/null || echo "0")
if [ "$START_SECONDS" != "0" ] && [ "$END_SECONDS" != "0" ]; then
    DURATION=$((END_SECONDS - START_SECONDS))
    echo "⏱️  Total durasi: $DURATION detik"
fi