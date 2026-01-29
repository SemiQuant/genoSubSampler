#!/bin/bash
#
# genoSubSampler - FASTQ subsampling tool for file size management
# Author: SemiQuant (JasonLimberis@ucsf.edu)
# GitHub: https://github.com/SemiQuant/genoSubSampler
#
# Description: Subsample paired-end FASTQ files to meet file size limits.
#              Optionally filter for M.tb aligned sequences using bowtie2.
#              Uses 7zip for compression (better compression than pigz/gzip).
#

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

# Determine script directory for relative paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

VERSION="1.0.0"
SEED=1987
DEFAULT_MAX_SIZE_MB=98
DEFAULT_READS=2000000
DEFAULT_THREADS=4
DEFAULT_REF="${SCRIPT_DIR}/Mycobacterium_tuberculosis_H37Rv_genome_v4"

# =============================================================================
# Functions
# =============================================================================

usage() {
    cat << EOF
genoSubSampler v${VERSION}

Subsample paired-end FASTQ files to meet file size constraints.

USAGE:
    $(basename "$0") -d <directory> [OPTIONS]

REQUIRED:
    -d, --dir <path>        Directory containing FASTQ files (*_R1_001.fastq.gz)

OPTIONS:
    -m, --max-size <MB>     Maximum file size in MB (default: ${DEFAULT_MAX_SIZE_MB})
    -r, --reads <int>       Initial number of reads for subsampling (default: ${DEFAULT_READS})
    -t, --threads <int>     Number of threads for bowtie2 (default: ${DEFAULT_THREADS})
    -f, --filter            Filter for M.tb sequences using bowtie2
    --ref <path>            Bowtie2 reference index (default: ${DEFAULT_REF})
    -h, --help              Show this help message
    -v, --version           Show version

EXAMPLES:
    # Basic subsampling
    $(basename "$0") -d /path/to/fastq

    # Subsample with custom size limit
    $(basename "$0") -d /path/to/fastq -m 50

    # Filter for M.tb sequences and subsample
    $(basename "$0") -d /path/to/fastq -f -t 8

OUTPUT:
    - Subsampled FASTQ files (compressed with 7zip)
    - subsampling_summary.txt with read counts and retention statistics

NOTE:
    Uses 7zip instead of pigz/gzip for slightly better compression ratios.

EOF
    exit 0
}

show_version() {
    echo "genoSubSampler v${VERSION}"
    exit 0
}

log_info() {
    echo "[INFO] $*"
}

log_error() {
    echo "[ERROR] $*" >&2
}

log_warn() {
    echo "[WARN] $*"
}

check_dependencies() {
    local missing=()
    
    for cmd in seqtk 7z bc; do
        if ! command -v "$cmd" &> /dev/null; then
            missing+=("$cmd")
        fi
    done
    
    if [[ "$FILTER_MTB" == true ]]; then
        for cmd in bowtie2 samtools bedtools bgzip; do
            if ! command -v "$cmd" &> /dev/null; then
                missing+=("$cmd")
            fi
        done
    fi
    
    if [[ ${#missing[@]} -gt 0 ]]; then
        log_error "Missing dependencies: ${missing[*]}"
        log_error "Run 'install.sh' to install required tools"
        exit 1
    fi
}

count_reads() {
    local file="$1"
    if [[ "$file" == *.gz ]]; then
        zcat "$file" | wc -l | awk '{print int($1/4)}'
    else
        wc -l < "$file" | awk '{print int($1/4)}'
    fi
}

get_file_size_mb() {
    du -m "$1" | cut -f 1
}

compress_with_7z() {
    local input="$1"
    local output="${input}.gz"
    7z a -tgzip -mx=9 -mmt="$THREADS" "$output" "$input" > /dev/null
    rm -f "$input"
}

subsample_fastq() {
    local r1_in="$1"
    local r2_in="$2"
    local num_reads="$3"
    local r1_out="$4"
    local r2_out="$5"
    
    seqtk sample -s "$SEED" "$r1_in" "$num_reads" > "$r1_out"
    7z a -tgzip -mx=9 -mmt="$THREADS" "${r1_out}.gz" "$r1_out" > /dev/null
    rm -f "$r1_out"
    
    seqtk sample -s "$SEED" "$r2_in" "$num_reads" > "$r2_out"
    7z a -tgzip -mx=9 -mmt="$THREADS" "${r2_out}.gz" "$r2_out" > /dev/null
    rm -f "$r2_out"
}

recompress_with_7z() {
    local file="$1"
    local uncompressed="${file%.gz}"
    
    gunzip -k "$file"
    rm -f "$file"
    7z a -tgzip -mx=9 -mmt="$THREADS" "$file" "$uncompressed" > /dev/null
    rm -f "$uncompressed"
}

filter_mtb_sequences() {
    local r1="$1"
    local r2="$2"
    local threads="$3"
    local ref="$4"
    
    log_info "Filtering for M.tb sequences using bowtie2"
    log_info "Reference: $ref"
    
    bowtie2 \
        --dovetail \
        --local \
        --minins 0 \
        --no-unal \
        --time \
        -p "$threads" \
        -x "$ref" \
        -1 "$r1" \
        -2 "$r2" \
        -S "out.sam" \
        --un-conc-gz "unaligned_%.fastq.gz" 2>bowtie2.log

    log_info "Alignment complete, sorting BAM files"
    
    samtools sort "out.sam" -o "coord_sorted.bam"
    samtools sort -n "coord_sorted.bam" -o "sorted.bam"
    
    if [[ ! -f "sorted.bam" ]]; then
        log_error "BAM sorting failed"
        return 1
    fi
    
    local r1_mtb="${r1/.fastq*/_mtb_only.fastq}"
    local r2_mtb="${r2/.fastq*/_mtb_only.fastq}"
    
    bedtools bamtofastq -i "sorted.bam" -fq "$r1_mtb" -fq2 "$r2_mtb" || {
        log_error "bedtools bamtofastq failed"
        return 1
    }
    
    # Cleanup intermediate files
    rm -f "sorted.bam" "coord_sorted.bam" "out.sam" "*.bai"
    
    # Compress output files
    if [[ -f "$r1_mtb" ]] && [[ -f "$r2_mtb" ]]; then
        bgzip "$r1_mtb"
        bgzip "$r2_mtb"
        bgzip -f "unaligned_1.fastq" 2>/dev/null || true
        bgzip -f "unaligned_2.fastq" 2>/dev/null || true
        
        echo "${r1_mtb}.gz|${r2_mtb}.gz"
    else
        log_error "FASTQ output files not created"
        return 1
    fi
}

# =============================================================================
# Argument Parsing
# =============================================================================

FILE_DIR=""
MAX_SIZE_MB="$DEFAULT_MAX_SIZE_MB"
READS="$DEFAULT_READS"
THREADS="$DEFAULT_THREADS"
FILTER_MTB=false
REF="$DEFAULT_REF"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dir)
            FILE_DIR="$2"
            shift 2
            ;;
        -m|--max-size)
            MAX_SIZE_MB="$2"
            shift 2
            ;;
        -r|--reads)
            READS="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -f|--filter)
            FILTER_MTB=true
            shift
            ;;
        --ref)
            REF="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        -v|--version)
            show_version
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            ;;
    esac
done

# =============================================================================
# Validation
# =============================================================================

if [[ -z "$FILE_DIR" ]]; then
    log_error "Directory is required. Use -d or --dir to specify."
    usage
fi

if [[ ! -d "$FILE_DIR" ]]; then
    log_error "Directory does not exist: $FILE_DIR"
    exit 1
fi

check_dependencies

# =============================================================================
# Main Processing
# =============================================================================

cd "$FILE_DIR"

SUMMARY_FILE="subsampling_summary.txt"
echo "Sample_ID,Original_R1_Reads,Original_R2_Reads,Subsampled_R1_Reads,Subsampled_R2_Reads,R1_Retention_%,R2_Retention_%,File_Size_R1_MB,File_Size_R2_MB" > "$SUMMARY_FILE"

log_info "Starting genoSubSampler v${VERSION}"
log_info "Directory: $FILE_DIR"
log_info "Max file size: ${MAX_SIZE_MB} MB"
log_info "Initial reads: $READS"
[[ "$FILTER_MTB" == true ]] && log_info "M.tb filtering: enabled"

# Count total files to process
total_files=$(find "$FILE_DIR" -maxdepth 1 -name "*_R1_001.fastq.gz" | wc -l)
current_file=0

for r1 in "$FILE_DIR"/*_R1_001.fastq.gz; do
    [[ ! -f "$r1" ]] && continue
    
    ((current_file++))
    r2="${r1/_R1_001/_R2_001}"
    sample_id=$(basename "$r1" | sed 's/_R1_001.fastq.gz//')
    
    log_info "Processing [$current_file/$total_files]: $sample_id"
    
    if [[ ! -f "$r2" ]]; then
        log_warn "R2 file not found for $sample_id, skipping"
        continue
    fi
    
    # Count original reads
    log_info "Counting reads in original files..."
    original_r1_reads=$(count_reads "$r1")
    original_r2_reads=$(count_reads "$r2")
    log_info "Original reads - R1: $original_r1_reads, R2: $original_r2_reads"
    
    # Optional M.tb filtering
    if [[ "$FILTER_MTB" == true ]]; then
        result=$(filter_mtb_sequences "$r1" "$r2" "$THREADS" "$REF")
        if [[ $? -eq 0 ]]; then
            r1=$(echo "$result" | cut -d'|' -f1)
            r2=$(echo "$result" | cut -d'|' -f2)
        else
            log_error "M.tb filtering failed for $sample_id"
            continue
        fi
    fi
    
    # Recompress with 7zip for better compression
    log_info "Recompressing with 7zip for optimal compression..."
    recompress_with_7z "$r1"
    recompress_with_7z "$r2"
    
    r1_size=$(get_file_size_mb "$r1")
    r2_size=$(get_file_size_mb "$r2")
    
    # Check if subsampling is needed
    if [[ "$r1_size" -le "$MAX_SIZE_MB" ]] && [[ "$r2_size" -le "$MAX_SIZE_MB" ]]; then
        log_info "Files already within size limit (R1: ${r1_size}MB, R2: ${r2_size}MB), skipping subsampling"
        
        final_r1_reads=$(count_reads "$r1")
        final_r2_reads=$(count_reads "$r2")
        
        echo "$sample_id,$original_r1_reads,$original_r2_reads,$final_r1_reads,$final_r2_reads,100.00,100.00,$r1_size,$r2_size" >> "$SUMMARY_FILE"
        continue
    fi
    
    # Subsample to target size
    log_info "Subsampling to approximately $READS reads..."
    r1_sub="${r1%.fastq.gz}_subsampled.fastq"
    r2_sub="${r2%.fastq.gz}_subsampled.fastq"
    
    reads_sub="$READS"
    subsample_fastq "$r1" "$r2" "$reads_sub" "$r1_sub" "$r2_sub"
    
    # Iteratively reduce reads if still too large
    while [[ "$(get_file_size_mb "${r1_sub}.gz")" -gt "$MAX_SIZE_MB" ]] || \
          [[ "$(get_file_size_mb "${r2_sub}.gz")" -gt "$MAX_SIZE_MB" ]]; do
        reads_sub=$((reads_sub - 50000))
        log_info "Files too large, trying $reads_sub reads..."
        rm -f "${r1_sub}.gz" "${r2_sub}.gz"
        subsample_fastq "$r1" "$r2" "$reads_sub" "$r1_sub" "$r2_sub"
    done
    
    final_r1_size=$(get_file_size_mb "${r1_sub}.gz")
    final_r2_size=$(get_file_size_mb "${r2_sub}.gz")
    log_info "Final sizes - R1: ${final_r1_size}MB, R2: ${final_r2_size}MB"
    
    # Calculate statistics
    subsampled_r1_reads=$(count_reads "${r1_sub}.gz")
    subsampled_r2_reads=$(count_reads "${r2_sub}.gz")
    r1_retention=$(echo "scale=2; $subsampled_r1_reads * 100 / $original_r1_reads" | bc)
    r2_retention=$(echo "scale=2; $subsampled_r2_reads * 100 / $original_r2_reads" | bc)
    
    echo "$sample_id,$original_r1_reads,$original_r2_reads,$subsampled_r1_reads,$subsampled_r2_reads,$r1_retention,$r2_retention,$final_r1_size,$final_r2_size" >> "$SUMMARY_FILE"
    
    # Cleanup original files
    rm -f "$r1" "$r2"
    
    log_info "Completed $sample_id"
done

# =============================================================================
# Summary Output
# =============================================================================

echo ""
echo "=============================================="
echo "           SUBSAMPLING SUMMARY"
echo "=============================================="

if [[ -f "$SUMMARY_FILE" ]] && [[ $(wc -l < "$SUMMARY_FILE") -gt 1 ]]; then
    echo "Summary saved to: $SUMMARY_FILE"
    echo ""
    column -t -s',' "$SUMMARY_FILE"
else
    echo "No FASTQ files were processed"
fi

log_info "genoSubSampler complete"
