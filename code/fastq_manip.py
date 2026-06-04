"""
Code untuk memanipulasi file fastq untuk edukasi
Author: Agus Wibowo

Contoh penggunaan:
python fastq_manip.py input_R1.fastq.gz output_R1.fastq.gz 100

input_R1.fastq.gz: file input R1
output_R1.fastq.gz: file output R1
"""


import random
import gzip
import sys
import os
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def degrade_quality(phred_str, start_pos=None, intensity=20):
    """Menurunkan nilai kualitas mulai dari posisi tertentu"""
    if start_pos is None:
        start_pos = len(phred_str) // 2  # Default: degradasi dari tengah read
    
    phred_list = list(phred_str)
    for i in range(start_pos, len(phred_list)):
        # Turunkan nilai kualitas (dalam ASCII)
        ascii_val = ord(phred_list[i])
        # Pastikan tidak turun di bawah '!' (ASCII 33 - nilai terendah dalam FASTQ)
        new_val = max(33, ascii_val - intensity)
        phred_list[i] = chr(new_val)
    
    return ''.join(phred_list)

def add_adapter(seq, qual, adapter="AGATCGGAAGAGC"):
    """Menambahkan sequence adapter pada akhir read"""
    adapter_qual = ''.join([chr(random.randint(33, 40)) for _ in range(len(adapter))])
    return seq + adapter, qual + adapter_qual

def add_n_bases(seq, qual, count=5):
    """Menambahkan N (base tidak teridentifikasi) secara acak"""
    seq_list = list(seq)
    qual_list = list(qual)
    
    for _ in range(count):
        pos = random.randint(0, len(seq_list) - 1)
        seq_list[pos] = 'N'
        qual_list[pos] = '!'  # Kualitas terendah untuk N
    
    return ''.join(seq_list), ''.join(qual_list)

def count_fastq_records(input_file):
    """Menghitung jumlah record dalam file FASTQ untuk progress bar"""
    print("Menghitung jumlah total read dalam file (mungkin memerlukan waktu)...")
    
    is_gzip = input_file.endswith(".gz")
    count = 0
    
    try:
        # Metode cepat untuk file non-gzip
        if not is_gzip:
            with open(input_file, 'r') as f:
                for i, line in enumerate(f):
                    pass
            count = (i + 1) // 4  # Setiap 4 baris = 1 read dalam FASTQ
        else:
            # Untuk file gzip, gunakan BioPython (lebih lambat)
            with gzip.open(input_file, 'rt') as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    count += 1
                    # Tampilkan progres selama penghitungan
                    if count % 100000 == 0:
                        sys.stdout.write(f"\rMenghitung read: {count:,}")
                        sys.stdout.flush()
        
        print(f"\rTotal read: {count:,}")
        return count
    except Exception as e:
        print(f"Error saat menghitung record: {e}")
        print("Melanjutkan tanpa progress bar...")
        return None

def estimate_subsample_rate(input_file, target_size_mb=50):
    """Estimasi rate subsample berdasarkan ukuran target"""
    # Mendapatkan ukuran file asli dalam MB
    original_size = os.path.getsize(input_file) / (1024 * 1024)
    
    # Menghitung proporsi yang dibutuhkan
    subsample_rate = target_size_mb / original_size
    
    # Memastikan rate tidak lebih dari 1.0 (100%)
    return min(subsample_rate, 1.0)

def print_progress_bar(iteration, total, prefix='Progress:', suffix='Complete', length=50, fill='█'):
    """Menampilkan progress bar di terminal"""
    percent = ("{0:.1f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    sys.stdout.write(f'\r{prefix} |{bar}| {percent}% {suffix}')
    sys.stdout.flush()
    # Print New Line on Complete
    if iteration == total:
        print()

def subsample_and_degrade_fastq(input_file, output_file, target_size_mb=50, 
                               degrade_percent=80, adapter_percent=25, n_base_percent=15):
    """Subsample dan manipulasi file FASTQ untuk menurunkan kualitasnya dan ukurannya"""
    # Menampilkan informasi awal
    print(f"File input: {input_file}")
    print(f"File output: {output_file}")
    print(f"Target ukuran: {target_size_mb} MB")
    
    # Menghitung rate subsample
    subsample_rate = estimate_subsample_rate(input_file, target_size_mb)
    print(f"Subsample rate: {subsample_rate:.4f} (target: {target_size_mb} MB)")
    
    # Mendapatkan jumlah total record untuk progress bar
    total_records = count_fastq_records(input_file)
    use_progress_bar = total_records is not None
    
    # Mendeteksi apakah file gzip
    is_gzip = input_file.endswith(".gz")
    in_handle = gzip.open(input_file, "rt") if is_gzip else open(input_file, "r")
    out_handle = gzip.open(output_file, "wt") if output_file.endswith(".gz") else open(output_file, "w")
    
    # Menghitung jumlah read yang diproses dan disimpan
    read_count = 0
    saved_count = 0
    start_time = time.time()
    last_update = start_time
    
    # Inisialisasi progress bar jika available
    if use_progress_bar:
        print_progress_bar(0, total_records, prefix='Processing:', suffix='Complete')
    
    for record in SeqIO.parse(in_handle, "fastq"):
        read_count += 1
        
        # Update progress bar setiap 1000 read atau setiap 1 detik
        current_time = time.time()
        if use_progress_bar and (read_count % 1000 == 0 or current_time - last_update >= 1):
            print_progress_bar(read_count, total_records, prefix='Processing:', 
                              suffix=f'Complete ({saved_count:,} saved)')
            last_update = current_time
        # Jika tidak menggunakan progress bar, tampilkan update setiap 100000 read
        elif not use_progress_bar and read_count % 100000 == 0:
            elapsed = current_time - start_time
            rate = read_count / elapsed if elapsed > 0 else 0
            sys.stdout.write(f"\rProcessed: {read_count:,} reads | Saved: {saved_count:,} | Rate: {rate:.1f} reads/sec")
            sys.stdout.flush()
        
        # Lakukan subsample berdasarkan rate yang ditentukan
        if random.random() > subsample_rate:
            continue
        
        saved_count += 1
        seq = str(record.seq)
        qual = record.letter_annotations["phred_quality"]
        qual_str = ''.join([chr(q + 33) for q in qual])  # Convert to ASCII
        
        # 1. Degradasi kualitas di akhir read
        if random.randint(1, 100) <= degrade_percent:
            degrade_pos = int(len(seq) * random.uniform(0.5, 0.8))
            qual_str = degrade_quality(qual_str, degrade_pos)
        
        # 2. Tambahkan adapter pada beberapa read
        if random.randint(1, 100) <= adapter_percent:
            seq, qual_str = add_adapter(seq, qual_str)
        
        # 3. Tambahkan N bases pada beberapa read
        if random.randint(1, 100) <= n_base_percent:
            seq, qual_str = add_n_bases(seq, qual_str)
        
        # Buat record baru dengan kualitas yang dimanipulasi
        new_record = SeqRecord(
            Seq(seq),
            id=record.id,
            name=record.name,
            description=record.description,
            letter_annotations={"phred_quality": [ord(c) - 33 for c in qual_str]}
        )
        
        SeqIO.write(new_record, out_handle, "fastq")
    
    # Finalisasi progress bar
    if use_progress_bar:
        print_progress_bar(total_records, total_records, prefix='Processing:', 
                         suffix=f'Complete ({saved_count:,} saved)')
    
    in_handle.close()
    out_handle.close()
    
    # Hitung waktu total
    total_time = time.time() - start_time
    minutes, seconds = divmod(total_time, 60)
    
    # Tampilkan statistik
    print(f"\nSelesai dalam {int(minutes)} menit {int(seconds)} detik")
    print(f"Total reads diproses: {read_count:,}")
    print(f"Reads tersimpan: {saved_count:,} ({saved_count/read_count*100:.2f}%)")
    print(f"File yang sudah dimanipulasi disimpan sebagai {output_file}")
    
    # Tampilkan ukuran akhir
    final_size = os.path.getsize(output_file) / (1024 * 1024)
    original_size = os.path.getsize(input_file) / (1024 * 1024)
    print(f"Ukuran file awal: {original_size:.2f} MB")
    print(f"Ukuran file akhir: {final_size:.2f} MB ({final_size/original_size*100:.2f}% dari asli)")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Penggunaan: python fastq_manip.py input.fastq[.gz] output.fastq[.gz] [target_size_mb]")
        print("Contoh: python fastq_manip.py input_R1.fastq.gz output_R1.fastq.gz 100")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Target ukuran file (default: 50 MB)
    target_size_mb = 50
    if len(sys.argv) > 3:
        target_size_mb = float(sys.argv[3])
    
    subsample_and_degrade_fastq(input_file, output_file, target_size_mb)