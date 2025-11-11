#!/bin/bash

# -----------------------------------------------------------------------------
# ExcludonFinder - Main Pipeline Script
# -----------------------------------------------------------------------------
# This script performs the core analysis: alignment, coverage calculation,
# and TU annotation. It is typically called by the ExcludonFinder wrapper.
# -----------------------------------------------------------------------------

# Stop the script immediately if any command fails.
set -e -o pipefail

# Automatically detect the path to the scripts directory.
if [ -z "$SCRIPTS_PATH" ]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    if [ -f "$SCRIPT_DIR/TUs_annotation.R" ]; then
        SCRIPTS_PATH="$SCRIPT_DIR"
        echo "Running from source code, using: $SCRIPTS_PATH" >&2
    else
        # Fallback for package installation.
        SCRIPT_DIR_PACKAGE="$(dirname "$(readlink -f "$(which ExcludonFinder 2>/dev/null || echo "$0")")")"
        PACKAGE_DIR="$(dirname "$SCRIPT_DIR_PACKAGE")"
        SCRIPTS_PATH="$PACKAGE_DIR/share/excludonfinder/scripts"
        echo "Trying package installation path: $SCRIPTS_PATH" >&2
    fi
    export SCRIPTS_PATH
fi

# Verify that the scripts directory and required R script exist.
if [ ! -d "$SCRIPTS_PATH" ] || [ ! -f "$SCRIPTS_PATH/TUs_annotation.R" ]; then
    echo "ERROR: Scripts directory or required R script not found at $SCRIPTS_PATH" >&2
    exit 1
fi

# Function to check if a command exists.
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to display the help message.
show_help() {
    cat << EOF
ExcludonFinder main processing script
Usage:
    main.sh -f <reference.fasta> -1 <reads_R1.fastq> [-2 <reads_R2.fastq>] -g <annotation.gff> [options]

Required arguments:
    -f <file>    Reference genome in FASTA format
    -1 <file>    Input FASTQ file (Read 1 for paired-end data)
    -g <file>    Annotation file in GFF format

Optional arguments:
    -2 <file>    Input FASTQ file (Read 2 for paired-end data)
    -t <float>   Coverage threshold (default: 0.5)
    -j <int>     Number of threads (default: 8)
    -l           Use long-read mode (uses minimap2 instead of bwa-mem2)
    -o <dir>     Supply a custom output dir (default: ./output)
    -k           Keep intermediate files (default: remove)
    -C           Run quality control checks
    -h, --help   Show this help message
EOF
}

# Handle the help flag first.
if [[ "$1" == "--help" ]] || [[ "$1" == "-h" ]]; then
    show_help
    exit 0
fi

# Check for required software dependencies.
for cmd in bwa-mem2 samtools minimap2 reformat.sh; do
  if ! command_exists "$cmd"; then
    echo "ERROR: Required command '$cmd' not found. Please ensure dependencies are installed." >&2
    exit 1
  fi
done

# Set default values for options.
threshold=0.5
N_threads=8
Paired="FALSE"
use_minimap2=false
output_dir="output"
keep_intermediate=false
run_control_script=false

# Parse command-line arguments.
while getopts "f:1:2:g:t:j:o:lkC" opt; do
  case ${opt} in
    f ) fasta_input=$(realpath "${OPTARG}");;
    1 ) fastq_input1=$(realpath "${OPTARG}");;
    2 ) fastq_input2=$(realpath "${OPTARG}");;
    g ) gff_input=$(realpath "${OPTARG}");;
    t ) threshold=$OPTARG;;
    j ) N_threads=$OPTARG;;
    o ) output_dir=$(realpath "${OPTARG}");;
    l ) use_minimap2=true;;
    k ) keep_intermediate=true;;
    C ) run_control_script=true;;
    \? ) echo "ERROR: Invalid option: -$OPTARG" 1>&2; exit 1;;
    : ) echo "ERROR: Option -$OPTARG requires an argument." 1>&2; exit 1;;
  esac
done

# Validate that mandatory files were provided.
if [ -z "$fasta_input" ] || [ -z "$gff_input" ] || [ -z "$fastq_input1" ]; then
  echo "ERROR: Missing required arguments. Use -h for help." >&2
  exit 1
fi

# Create output directories.
mkdir -p "$output_dir"/{alignment,coverage_data,index,Excludons}

# ✅ RESCUED LOGIC: Robustly detect and handle paired-end data.
# First, check if two files were provided explicitly.
if [ -n "$fastq_input2" ]; then
    Paired="TRUE"
else
    # If not, check if the single input file is interleaved.
    first_header=$(head -n 1 "$fastq_input1" 2>/dev/null)
    fifth_header=$(head -n 5 "$fastq_input1" 2>/dev/null | tail -n 1)
    if [[ -n "$first_header" && -n "$fifth_header" ]]; then
        id1=$(echo "$first_header" | cut -d' ' -f1)
        id2=$(echo "$fifth_header" | cut -d' ' -f1)
        if [[ "$id1" == "$id2" ]]; then
            echo "Interleaved paired-end data detected. Splitting..." >&2
            Paired="TRUE"
            mkdir -p "$output_dir"/temp_split
            filename=$(basename "$fastq_input1")
            sample=${filename%.fastq*}

            # Split the interleaved file into R1 and R2.
            reformat.sh in="$fastq_input1" \
                out1="$output_dir/temp_split/${sample}_R1.fastq" \
                out2="$output_dir/temp_split/${sample}_R2.fastq" \
                ow=t >/dev/null 2>&1

            # Update variables to point to the new split files.
            fastq_input2="$output_dir/temp_split/${sample}_R2.fastq"
            fastq_input1="$output_dir/temp_split/${sample}_R1.fastq"
            echo "Split completed: R1 and R2 files created" >&2
        fi
    fi
fi

# ✅ FIX: Prepare GFF file - create a local copy and normalize attributes
cp "$gff_input" "$output_dir"/$(basename "$gff_input")
gff_input="$output_dir"/$(basename "$gff_input")

# 1. Fix GFF header if corrupted (remove extra spaces in version directive)
if [[ "$OSTYPE" == "darwin"* ]]; then
  sed -i '' '1s/##gff-version 3.*/##gff-version 3/' "$gff_input"
else
  sed -i '1s/##gff-version 3.*/##gff-version 3/' "$gff_input"
fi

# 2. Normalize attributes for R script compatibility
# - Change ID= to gene_id= (only word boundaries to avoid false matches)
# - Change name= to gene= (both after semicolon and at start of attributes)
if [[ "$OSTYPE" == "darwin"* ]]; then
  sed -i '' -e 's/\bID=/gene_id=/g' \
            -e 's/;name=/;gene=/g' \
            -e 's/\tname=/\tgene=/g' "$gff_input"
else
  sed -i -e 's/\bID=/gene_id=/g' \
         -e 's/;name=/;gene=/g' \
         -e 's/\tname=/\tgene=/g' "$gff_input"
fi

# Define sample name and input files for the aligner.
if [ "$Paired" = "FALSE" ]; then
  filename=$(basename "$fastq_input1")
  sample=${filename%.fastq*}
  input_fastq_files=("$fastq_input1")
else
  filename1=$(basename "$fastq_input1")
  # Handle different naming conventions like .fastq, .fq, _R1.fastq, _1.fastq etc.
  sample=$(basename "$fastq_input1" | sed -E 's/(_R1|_1)?\.(fastq|fq)(\.gz)?$//')
  input_fastq_files=("$fastq_input1" "$fastq_input2")
fi

# Perform alignment.
echo "STARTING_ALIGNMENT" >&2
bam_output_path="$output_dir/alignment/${sample}_sorted.bam"
if $use_minimap2; then
    echo "Data type: Long reads detected" >&2
    minimap2 -t "$N_threads" -a -x map-ont "$fasta_input" "$fastq_input1" 2>/dev/null | \
      samtools view -@ "$N_threads" -bh - 2>/dev/null | \
      samtools sort -@ "$N_threads" -o "$bam_output_path" 2>/dev/null
else
    index_prefix="$output_dir/index/$(basename "$fasta_input")"
    bwa-mem2 index -p "$index_prefix" "$fasta_input" >/dev/null 2>&1
    bwa-mem2 mem -t "$N_threads" "$index_prefix" "${input_fastq_files[@]}" 2>/dev/null | \
      samtools view -@ "$N_threads" -bh - 2>/dev/null | \
      samtools sort -@ "$N_threads" -o "$bam_output_path" 2>/dev/null
fi

samtools index "$bam_output_path"

# Get alignment statistics.
mapped_percent=$(samtools flagstat "$bam_output_path" | grep "mapped" | head -n1 | awk -F'[(%]' '{print $2}')
paired_percent=$(samtools flagstat "$bam_output_path" | grep "properly paired" | awk -F'[(%]' '{print $2}')
echo "MAPPED_READS:$mapped_percent" >&2
echo "PAIRED_READS:$paired_percent" >&2

# Calculate coverage depth per strand.
echo "CALCULATING_DEPTH" >&2
if [ "$Paired" = "TRUE" ]; then
    # Paired-end logic for strand separation.
    samtools view -h -f 0x40 -F 0x10 "$bam_output_path" 2>/dev/null | samtools sort -@ "$N_threads" -o "$output_dir/alignment/plus_1.bam" 2>/dev/null
    samtools view -h -f 0x80 -f 0x10 "$bam_output_path" 2>/dev/null | samtools sort -@ "$N_threads" -o "$output_dir/alignment/plus_2.bam" 2>/dev/null
    samtools merge -f "$output_dir/alignment/plus_merged.bam" "$output_dir/alignment/plus_1.bam" "$output_dir/alignment/plus_2.bam" >/dev/null 2>&1

    samtools view -h -f 0x40 -f 0x10 "$bam_output_path" 2>/dev/null | samtools sort -@ "$N_threads" -o "$output_dir/alignment/minus_1.bam" 2>/dev/null
    samtools view -h -f 0x80 -F 0x10 "$bam_output_path" 2>/dev/null | samtools sort -@ "$N_threads" -o "$output_dir/alignment/minus_2.bam" 2>/dev/null
    samtools merge -f "$output_dir/alignment/minus_merged.bam" "$output_dir/alignment/minus_1.bam" "$output_dir/alignment/minus_2.bam" >/dev/null 2>&1
else
    # Single-end logic.
    samtools view -bh -F 0x14 "$bam_output_path" > "$output_dir/alignment/plus_merged.bam"
    samtools view -bh -f 0x10 -F 0x4 "$bam_output_path" > "$output_dir/alignment/minus_merged.bam"
fi

for strand in plus minus; do
    samtools sort -@ "$N_threads" "$output_dir/alignment/${strand}_merged.bam" -o "$output_dir/alignment/${strand}_sorted.bam" 2>/dev/null
    samtools index "$output_dir/alignment/${strand}_sorted.bam"
    samtools depth -a "$output_dir/alignment/${strand}_sorted.bam" > "$output_dir/coverage_data/${sample}_${strand}_depth.txt" 2>/dev/null
done

echo "DEPTH COMPLETED" >&2
echo "STARTING_ANNOTATION" >&2

# Run the R script for TU annotation and excludon finding.
Rscript "$SCRIPTS_PATH/TUs_annotation.R" "$fasta_input" "$gff_input" "$output_dir" "$keep_intermediate" "$sample" "$threshold" "$N_threads"

# Clean up intermediate files unless the -k flag is used.
if [ "$keep_intermediate" = false ]; then
    rm -r "$output_dir/alignment/"
    rm -r "$output_dir/coverage_data"
    rm -r "$output_dir/index"
    rm -f "$gff_input"
    rm -rf "$output_dir/temp_split"
fi

echo "### The analysis is finished ###"
